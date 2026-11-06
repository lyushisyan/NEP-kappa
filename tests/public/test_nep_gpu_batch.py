"""GPU dispatch, scientific ordering and cache contracts without a GPU/Slurm."""

import json
from pathlib import Path

from ase import Atoms
from ase.calculators.lj import LennardJones
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read
import calorine.tools
import h5py
import numpy as np
import pytest
import yaml

from nepkappa.application import execute_workflow_command
from nepkappa.config import parse_workflow_args
from nepkappa.force_worker import run_task_file
from nepkappa.nep_batch import predict_forces, resolve_nep_executable
from nepkappa.workflow import NEPPhononWorkflow


@pytest.fixture
def gpu_input(tmp_path, monkeypatch):
    monkeypatch.delenv("NEPKAPPA_FORCE_COLLECT", raising=False)
    atoms = Atoms("Ar2", positions=[[0, 0, 0], [4, 0, 0]], cell=[20]*3, pbc=True)
    poscar = tmp_path / "POSCAR"
    atoms.write(poscar, format="vasp")
    model = tmp_path / "nep.txt"
    model.write_text("synthetic model identity; calculator mocked\n")
    executable = tmp_path / "nep executable"
    executable.write_text("#!/bin/sh\nexit 99\n")  # Must never be launched.
    executable.chmod(0o755)
    data = {
        "structure": {"poscar": str(poscar), "relaxation": False},
        "calculator": {
            "name": "nep", "model": str(model), "nep-force-backend": "gpu-batch",
            "nep-batch-size": 2, "nep-command": str(executable),
        },
        "force-constant": {"dim-fc2": [1, 1, 1], "dim-fc3": [1, 1, 1]},
        "output": {"result-dir": str(tmp_path / "results"), "progress": False},
    }
    path = tmp_path / "input.yaml"
    path.write_text(yaml.safe_dump(data))
    calls = []

    def predict(structures, *, model, command):
        assert Path(model).is_absolute()
        assert command == str(executable)
        calls.append([s.copy() for s in structures])
        results = []
        for atoms in structures:
            result = atoms.copy()
            result.calc = LennardJones(epsilon=0.0103, sigma=3.4, rc=8.5)
            forces = result.get_forces()
            result.calc = SinglePointCalculator(result, forces=forces)
            results.append(result)
        return results

    monkeypatch.setattr(calorine.tools, "batch_predict_properties", predict, raising=False)
    monkeypatch.setattr("nepkappa.workflow.CPUNEP", lambda _: LennardJones(
        epsilon=0.0103, sigma=3.4, rc=8.5
    ))
    return path, calls


def make_jobs(cfg, count=5):
    jobs = []
    for i in range(count):
        atoms = read(cfg.poscar)
        atoms.positions[1, 0] += i * 0.01
        jobs.append((atoms, "fc2", i + 1))
    return jobs


def test_parse_gpu_options_and_cpu_defaults(gpu_input):
    path, _ = gpu_input
    cfg = parse_workflow_args(path, command="fc2fc3")
    assert cfg.sections.calculator.nep_force_backend == "gpu-batch"
    assert cfg.sections.calculator.nep_batch_size == 2
    data = yaml.safe_load(path.read_text())
    for key in ("nep-force-backend", "nep-batch-size", "nep-command"):
        data["calculator"].pop(key)
    path.write_text(yaml.safe_dump(data))
    cfg = parse_workflow_args(path, command="fc2")
    assert cfg.nep_force_backend == "cpu"
    assert cfg.nep_batch_size == 128
    assert cfg.nep_command == "nep"


@pytest.mark.parametrize("size", [0, -1, 1.5, True])
def test_reject_invalid_batch_size(gpu_input, size):
    path, _ = gpu_input
    data = yaml.safe_load(path.read_text())
    data["calculator"]["nep-batch-size"] = size
    path.write_text(yaml.safe_dump(data))
    with pytest.raises((SystemExit, ValueError)):
        parse_workflow_args(path, command="fc2")


@pytest.mark.parametrize("option,value", [
    ("qha", True), ("sscha", True), ("four-phonon", True),
    ("fc3-backend", "thirdorder"),
])
def test_reject_unimplemented_workflows_before_execution(gpu_input, option, value):
    path, _ = gpu_input
    data = yaml.safe_load(path.read_text())
    data["force-constant"][option] = value
    path.write_text(yaml.safe_dump(data))
    with pytest.raises(SystemExit, match="2"):
        parse_workflow_args(path, command="fc2fc3")
    assert not Path(data["output"]["result-dir"]).exists()


def test_batching_cache_hits_and_original_order(gpu_input):
    path, calls = gpu_input
    cfg = parse_workflow_args(path, command="fc2")
    workflow = NEPPhononWorkflow(cfg)
    jobs = make_jobs(cfg)
    first = workflow._calculate_nep_batch_cached(jobs)
    assert [len(c) for c in calls] == [2, 2, 1]
    assert not np.array_equal(first[0], first[1])
    calls.clear()
    second = NEPPhononWorkflow(cfg)._calculate_nep_batch_cached(jobs)
    np.testing.assert_array_equal(first, second)
    assert calls == []
    # Reordered input and a hole in the cache must preserve each job's identity.
    workflow._force_job_dir("fc2", 3).joinpath("job-input.sha256").unlink()
    third = NEPPhononWorkflow(cfg)._calculate_nep_batch_cached(jobs[::-1])
    assert [len(c) for c in calls] == [1]
    np.testing.assert_array_equal(first[::-1], third)


@pytest.mark.parametrize("change", ["model", "executable", "batch-size", "cpu-cache"])
def test_gpu_cache_tracks_backend_inputs(gpu_input, change):
    path, calls = gpu_input
    cfg = parse_workflow_args(path, command="fc2")
    jobs = make_jobs(cfg, 1)
    if change == "cpu-cache":
        cfg.nep_force_backend = "cpu"
        NEPPhononWorkflow(cfg)._calculate_forces_cached(*jobs[0])
        cfg.nep_force_backend = "gpu-batch"
    else:
        NEPPhononWorkflow(cfg)._calculate_nep_batch_cached(jobs)
        calls.clear()
        if change == "batch-size":
            cfg.nep_batch_size += 1
        else:
            target = Path(cfg.nep_model if change == "model" else cfg.nep_command)
            target.write_text(target.read_text() + "\n# changed\n")
    workflow = NEPPhononWorkflow(cfg)
    workflow._calculate_nep_batch_cached(jobs)
    assert len(calls) == 1
    record = json.loads((workflow._force_job_dir("fc2", 1) / "job-input.json").read_text())
    identity = record["inputs"]["calculator_inputs"]
    assert identity["force_backend"] == "gpu-batch"
    assert identity["executable"]["sha256"]


@pytest.mark.parametrize("failure", ["count", "order", "shape", "nan", "executable"])
def test_invalid_batch_never_publishes_forces(gpu_input, monkeypatch, failure):
    path, calls = gpu_input
    cfg = parse_workflow_args(path, command="fc2")
    predict = calorine.tools.batch_predict_properties

    def broken(*args, **kwargs):
        if failure == "executable":
            raise RuntimeError("nep exited unsuccessfully")
        results = predict(*args, **kwargs)
        if failure == "count":
            return results[:-1]
        if failure == "order":
            return results[::-1]
        results[-1].calc.results["forces"] = (
            np.zeros((1, 3)) if failure == "shape" else np.full((2, 3), np.nan)
        )
        return results

    monkeypatch.setattr(calorine.tools, "batch_predict_properties", broken)
    workflow = NEPPhononWorkflow(cfg)
    with pytest.raises(RuntimeError):
        workflow._calculate_nep_batch_cached(make_jobs(cfg, 2))
    assert not list(workflow.output_dir.rglob("forces.npy"))
    assert not list(workflow.output_dir.rglob("job-input.sha256"))


def test_interrupted_batch_update_recovers_without_stale_certification(gpu_input, monkeypatch):
    path, calls = gpu_input
    cfg = parse_workflow_args(path, command="fc2")
    jobs = make_jobs(cfg, 2)
    workflow = NEPPhononWorkflow(cfg)
    reference = workflow._calculate_nep_batch_cached(jobs)
    directory = workflow._force_job_dir("fc2", 1)
    marker = directory / "job-input.sha256"
    original_digest = marker.read_text()
    cfg.nep_batch_size = 3
    replace = Path.replace
    predict = calorine.tools.batch_predict_properties

    def different(*args, **kwargs):
        results = predict(*args, **kwargs)
        for atoms in results:
            atoms.calc.results["forces"] *= 2
        return results

    def interrupt(source, target):
        result = replace(source, target)
        if Path(target) == directory / "forces.npy":
            raise KeyboardInterrupt("interrupted force publication")
        return result

    with monkeypatch.context() as patch:
        patch.setattr(calorine.tools, "batch_predict_properties", different)
        patch.setattr(Path, "replace", interrupt)
        with pytest.raises(KeyboardInterrupt):
            NEPPhononWorkflow(cfg)._calculate_nep_batch_cached(jobs)
    assert not marker.exists()
    np.testing.assert_array_equal(np.load(directory / "forces.npy"), 2 * reference[0])
    cfg.nep_batch_size = 2
    calls.clear()
    recovered = NEPPhononWorkflow(cfg)._calculate_nep_batch_cached(jobs)
    assert [len(c) for c in calls] == [1]
    np.testing.assert_array_equal(reference, recovered)
    assert marker.read_text() == original_digest


def test_missing_gpu_dependencies_fail_without_cpu_fallback(gpu_input, monkeypatch):
    path, _ = gpu_input
    cfg = parse_workflow_args(path, command="fc2")
    with pytest.raises(FileNotFoundError, match="NEP GPU executable"):
        resolve_nep_executable(str(path.parent / "missing-nep"))
    monkeypatch.delattr(calorine.tools, "batch_predict_properties")
    with pytest.raises(RuntimeError, match="Calorine >= 3.5"):
        predict_forces([read(cfg.poscar)], model=cfg.nep_model, command=cfg.nep_command)


def test_real_hiphive_batch_fitting_matches_cpu(gpu_input):
    path, calls = gpu_input
    data = yaml.safe_load(path.read_text())
    data["force-constant"].update({
        "use-hiphive": True, "n-structures": 12,
        "rattle-std": 0.01, "min-dist": 2.0, "cutoffs": [5.0, 5.0],
    })
    path.write_text(yaml.safe_dump(data))
    cfg = parse_workflow_args(path, command="fc2fc3")
    execute_workflow_command("fc2fc3", cfg)
    assert [len(c) for c in calls] == [2] * 6
    reference = Path(cfg.result_dir)
    cfg.nep_force_backend = "cpu"
    cfg.result_dir = str(path.parent / "cpu-results")
    execute_workflow_command("fc2fc3", cfg)
    for filename, dataset in (("fc2.hdf5", "force_constants"), ("fc3.hdf5", "fc3")):
        with h5py.File(reference / filename) as gpu, h5py.File(Path(cfg.result_dir) / filename) as cpu:
            np.testing.assert_allclose(gpu[dataset][:], cpu[dataset][:], rtol=0, atol=1e-12)


@pytest.mark.parametrize("cutoff", [None, 3.0])
def test_real_fc2_fc3_fitting_matches_cpu_and_slurm_collection(gpu_input, monkeypatch, cutoff):
    path, calls = gpu_input
    cfg = parse_workflow_args(path, command="fc2fc3")
    cfg.cutoff_fc3 = cutoff
    execution = execute_workflow_command("fc2fc3", cfg)
    assert execution.status == "complete"
    assert calls and max(map(len, calls)) <= 2
    reference = {}
    for order in (2, 3):
        with h5py.File(Path(cfg.result_dir) / f"fc{order}.hdf5") as handle:
            reference[order] = handle[f"fc{order}" if order == 3 else "force_constants"][:]
    cfg.nep_force_backend = "cpu"
    cfg.result_dir = str(path.parent / "cpu-results")
    execute_workflow_command("fc2fc3", cfg)
    for order, expected in reference.items():
        with h5py.File(Path(cfg.result_dir) / f"fc{order}.hdf5") as handle:
            actual = handle[f"fc{order}" if order == 3 else "force_constants"][:]
        np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-12)

    data = yaml.safe_load(path.read_text())
    if cutoff:
        data["force-constant"]["cutoff-fc3"] = cutoff
    data["output"]["result-dir"] = str(path.parent / "slurm-results")
    data["parallel"] = {"force-constant": {
        "backend": "slurm", "jobs": 2, "submit": False,
        "extra-sbatch": ["--gres=gpu:1"],
    }}
    path.write_text(yaml.safe_dump(data))
    cfg = parse_workflow_args(path, command="fc2fc3")
    execution = execute_workflow_command("fc2fc3", cfg)
    assert execution.status == "prepared"
    work = Path(cfg.result_dir) / "force-slurm"
    assert "#SBATCH --gres=gpu:1" in (work / "array.sh").read_text()
    for task in sorted((work / "tasks").glob("*.json")):
        run_task_file(path, task)
    calls.clear()
    monkeypatch.setenv("NEPKAPPA_FORCE_COLLECT", "1")
    assert execute_workflow_command("fc2fc3", cfg).status == "complete"
    assert calls == []  # Collection fits cached forces, no second GPU evaluation.
    for order, expected in reference.items():
        with h5py.File(Path(cfg.result_dir) / f"fc{order}.hdf5") as handle:
            actual = handle[f"fc{order}" if order == 3 else "force_constants"][:]
        np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-12)
