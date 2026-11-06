"""Opt-in warmed CPUNEP/GPU force throughput benchmark (no FC fitting or caches)."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess
import time
from unittest.mock import patch

from ase.io import read
from calorine.calculators import CPUNEP
import numpy as np

from nepkappa.nep_batch import predict_forces, resolve_nep_executable
from nepkappa.provenance import file_identity, installed_versions


def make_structures(primitive, dimension, count, seed):
    """Reproducible small random displacements, identical for both backends."""
    ideal = primitive.repeat((dimension,) * 3)
    rng = np.random.default_rng(seed)
    structures = []
    for _ in range(count):
        atoms = ideal.copy()
        atoms.positions += rng.normal(0.0, 0.01, size=(len(atoms), 3))
        structures.append(atoms)
    return structures


def cpu_forces(structures, calculator):
    calculator.reset()
    forces = []
    for atoms in structures:
        atoms.calc = calculator
        forces.append(atoms.get_forces())
    return np.asarray(forces)


def gpu_forces(structures, model, command, batch_size):
    """Time the exact public batch adapter, including its file I/O and checks."""
    run = subprocess.run
    process_seconds = []

    def timed_run(*args, **kwargs):
        start = time.perf_counter()
        result = run(*args, **kwargs)
        process_seconds.append(time.perf_counter() - start)
        return result

    forces = []
    with patch.object(subprocess, "run", timed_run):
        for first in range(0, len(structures), batch_size):
            forces.extend(predict_forces(
                structures[first:first + batch_size], model=model, command=command
            ))
    return np.asarray(forces), process_seconds


def compare(cpu, gpu):
    if cpu.shape != gpu.shape or not np.isfinite(cpu).all() or not np.isfinite(gpu).all():
        raise RuntimeError("Invalid force shapes or non-finite CPU/GPU results.")
    difference = gpu - cpu
    return {
        "max_abs_eV_per_A": float(np.max(np.abs(difference))),
        "rms_eV_per_A": float(np.sqrt(np.mean(difference**2))),
        "relative_l2": float(np.linalg.norm(difference) / np.linalg.norm(cpu)),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structure", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--nep-command", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--seed", type=int, default=20261106)
    parser.add_argument("--force-atol", type=float, default=1e-4)
    args = parser.parse_args()
    if args.repeats < 1 or args.force_atol <= 0:
        parser.error("repeats and force-atol must be positive")
    root = Path(args.output).resolve()
    root.mkdir(parents=True, exist_ok=False)
    model = str(Path(args.model).resolve())
    command = resolve_nep_executable(args.nep_command)
    primitive = read(args.structure)
    calculator = CPUNEP(model)
    report = {
        "finished": False,
        "versions": installed_versions(),
        "structure": file_identity(args.structure),
        "model": file_identity(model),
        "executable": file_identity(command),
        "adapter": file_identity(predict_forces.__code__.co_filename),
        "script": file_identity(__file__),
        "seed": args.seed,
        "displacement_sigma_A": 0.01,
        "force_atol_eV_per_A": args.force_atol,
        "repeats": args.repeats,
        "environment": {key: os.environ.get(key) for key in (
            "SLURM_JOB_ID", "SLURM_JOB_PARTITION", "SLURM_JOB_NODELIST",
            "SLURM_CPUS_PER_TASK", "CUDA_VISIBLE_DEVICES", "OMP_NUM_THREADS",
            "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "RAYON_NUM_THREADS", "TMPDIR",
        )},
        "timing_scope": (
            "Warmed single-process CPUNEP versus one GPU; model creation and structure "
            "generation excluded. GPU total includes Python preparation/parsing, "
            "temporary files, nep process startup, CUDA setup and evaluation. "
            "nep_process_seconds also includes nep-side file I/O and CUDA startup: "
            "it is NOT pure GPU kernel time. No force cache or FC fitting."
        ),
        "cases": [],
    }
    path = root / "benchmark.json"
    cases = [(2, 128, 128), (4, 128, 128), (6, 128, 128),
             (8, 128, 128), (8, 512, 128), (8, 512, 512)]
    for dimension, count, batch_size in cases:
        structures = make_structures(primitive, dimension, count, args.seed)
        digest = hashlib.sha256()
        for atoms in structures:
            digest.update(atoms.positions.tobytes())
        case = {
            "dimension": [dimension] * 3, "atoms": len(structures[0]),
            "structures": count, "batch_size": batch_size,
            "positions_sha256": digest.hexdigest(), "runs": [],
        }
        print(f"Warming {case['atoms']} atoms x {count} structures, batch={batch_size}", flush=True)
        cpu_forces(structures[:1], calculator)
        gpu_forces(structures[:batch_size], model, command, batch_size)
        for repeat in range(args.repeats):
            measured = {"repeat": repeat, "order": ["cpu", "gpu"] if repeat % 2 == 0 else ["gpu", "cpu"]}
            for backend in measured["order"]:
                started = time.perf_counter()
                if backend == "cpu":
                    cpu = cpu_forces(structures, calculator)
                else:
                    gpu, process_seconds = gpu_forces(structures, model, command, batch_size)
                measured[backend + "_seconds"] = time.perf_counter() - started
            measured["nep_process_seconds"] = process_seconds
            measured["gpu_python_io_seconds"] = measured["gpu_seconds"] - sum(process_seconds)
            measured["errors"] = compare(cpu, gpu)
            case["runs"].append(measured)
            print(json.dumps(measured), flush=True)
        case["cpu_median_seconds"] = statistics.median(x["cpu_seconds"] for x in case["runs"])
        case["gpu_median_seconds"] = statistics.median(x["gpu_seconds"] for x in case["runs"])
        case["speedup"] = case["cpu_median_seconds"] / case["gpu_median_seconds"]
        case["accuracy_passed"] = all(
            x["errors"]["max_abs_eV_per_A"] <= args.force_atol for x in case["runs"]
        )
        report["cases"].append(case)
        path.write_text(json.dumps(report, indent=2) + "\n")
        print(f"Case complete: speedup={case['speedup']:.3f}, accuracy={case['accuracy_passed']}", flush=True)
    report["finished"] = True
    report["accuracy_passed"] = all(case["accuracy_passed"] for case in report["cases"])
    path.write_text(json.dumps(report, indent=2) + "\n")
    return 0 if report["accuracy_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
