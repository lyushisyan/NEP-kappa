"""Opt-in real GPU/CPU FC2/FC3 comparison; run inside a GPU allocation.

This is not collected by pytest. It needs a physical NEP model and GPUMD's
nep executable. It leaves both input files, force caches and a JSON report.
"""

import argparse
import json
from pathlib import Path
import time

import h5py
import numpy as np
import yaml

from nepkappa.application import execute_workflow_command
from nepkappa.config import parse_workflow_args
from nepkappa.provenance import installed_versions


def differences(first, second):
    if first.shape != second.shape or not np.isfinite(first).all() or not np.isfinite(second).all():
        raise RuntimeError("CPU/GPU comparison has invalid shapes or non-finite values.")
    delta = second - first
    return {
        "shape": list(first.shape),
        "max_abs": float(np.max(np.abs(delta))),
        "rms": float(np.sqrt(np.mean(delta**2))),
        "cpu_max_abs": float(np.max(np.abs(first))),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structure", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--nep-command", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--dim", type=int, nargs=3, default=[2, 2, 2])
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--force-atol", type=float, default=1e-4, help="eV/Angstrom")
    parser.add_argument("--fc2-atol", type=float, default=5e-3, help="eV/Angstrom^2")
    parser.add_argument("--fc3-atol", type=float, default=0.2, help="eV/Angstrom^3")
    args = parser.parse_args()
    root = Path(args.output).resolve()
    root.mkdir(parents=True, exist_ok=False)
    report = {"versions": installed_versions(), "seconds": {}, "comparisons": {}}
    for backend in ("cpu", "gpu-batch"):
        data = {
            "structure": {"poscar": str(Path(args.structure).resolve()), "relaxation": False},
            "calculator": {
                "name": "nep", "model": str(Path(args.model).resolve()),
                "nep-force-backend": backend, "nep-batch-size": args.batch_size,
                "nep-command": args.nep_command,
            },
            "force-constant": {"dim-fc2": args.dim, "dim-fc3": args.dim},
            "output": {"result-dir": str(root / backend), "progress": False},
        }
        path = root / f"{backend}.yaml"
        path.write_text(yaml.safe_dump(data))
        cfg = parse_workflow_args(path, command="fc2fc3")
        started = time.perf_counter()
        execution = execute_workflow_command("fc2fc3", cfg)
        if execution.status != "complete":
            raise RuntimeError(f"Unexpected deferred {backend} calculation.")
        report["seconds"][backend] = time.perf_counter() - started

    first = root / "cpu"
    second = root / "gpu-batch"
    paths = sorted(first.rglob("forces.npy"))
    if not paths:
        raise RuntimeError("CPU run produced no force caches.")
    if {p.relative_to(first) for p in paths} != {
        p.relative_to(second) for p in second.rglob("forces.npy")
    }:
        raise RuntimeError("CPU/GPU force-job sets differ.")
    report["force_jobs"] = len(paths)
    report["comparisons"]["forces"] = differences(
        np.concatenate([np.load(p, allow_pickle=False).ravel() for p in paths]),
        np.concatenate([np.load(second / p.relative_to(first), allow_pickle=False).ravel() for p in paths]),
    )
    for order, dataset in ((2, "force_constants"), (3, "fc3")):
        with h5py.File(first / f"fc{order}.hdf5") as cpu, h5py.File(second / f"fc{order}.hdf5") as gpu:
            report["comparisons"][f"fc{order}"] = differences(cpu[dataset][:], gpu[dataset][:])
    limits = {"forces": args.force_atol, "fc2": args.fc2_atol, "fc3": args.fc3_atol}
    report["absolute_tolerances"] = limits
    report["passed"] = all(report["comparisons"][key]["max_abs"] <= limit for key, limit in limits.items())

    # A second GPU workflow must reuse forces, even with the predictor disabled.
    import calorine.tools

    def unexpected_prediction(*args, **kwargs):
        raise RuntimeError("A completed GPU cache was unexpectedly recalculated.")

    calorine.tools.batch_predict_properties = unexpected_prediction
    execute_workflow_command("fc2fc3", cfg)
    report["cache_reuse_passed"] = True
    output = root / "comparison.json"
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(output.read_text())
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
