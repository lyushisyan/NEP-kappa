"""Real force/stage/cache integration on synthetic atoms; no external engines."""

import json
from pathlib import Path
import shlex
import sys
from types import SimpleNamespace

from ase import Atoms
from ase.io import read
import h5py
import numpy as np
import pytest
import yaml

from nepkappa.application import execute_workflow_command
from nepkappa.cli import build_parser, main
from nepkappa.config import parse_selected_command
from nepkappa.config import parse_workflow_args
from nepkappa.force_worker import run_task_file
from nepkappa.workflow import NEPPhononWorkflow


@pytest.fixture
def force_input(tmp_path, monkeypatch):
    for name in ("NEPKAPPA_FORCE_COLLECT", "NEPKAPPA_RESUME_AFTER"):
        monkeypatch.delenv(name, raising=False)
    atoms = Atoms("Ar2", positions=[[0, 0, 0], [4, 0, 0]], cell=[20]*3, pbc=True)
    structure = tmp_path / "POSCAR"
    atoms.write(structure, format="vasp")
    path = tmp_path / "force input.yaml"
    path.write_text(yaml.safe_dump({
        "structure": {"poscar": str(structure), "dimensionality": 3, "relaxation": False},
        "calculator": {"name": "ase", "factory": "ase.calculators.lj:LennardJones",
                       "kwargs": {"epsilon": 0.0103, "sigma": 3.4, "rc": 8.5}},
        "force-constant": {"dim-fc2": [1, 1, 1], "dim-fc3": [1, 1, 1]},
        "output": {"result-dir": str(tmp_path / "batch results"), "progress": False},
    }))
    return path


@pytest.mark.parametrize("command", ["fc2", "run"])
@pytest.mark.parametrize("submit", [False, True])
def test_local_and_slurm_collected_fc2_agree(force_input, tmp_path, monkeypatch, command, submit):
    from nepkappa import scheduler
    from nepkappa.transport import Phono3pyTransportWorkflow

    local = parse_workflow_args(force_input, command="fc2")
    local.result_dir = str(tmp_path / "local results")
    execution = execute_workflow_command("fc2", local)
    assert [name for name, _ in execution.stage_timings] == ["Finite displacement"]
    with h5py.File(Path(local.result_dir) / "fc2.hdf5") as handle:
        reference = handle["force_constants"][:]
    assert reference.shape == (2, 2, 3, 3)
    assert np.isfinite(reference).all() and np.linalg.norm(reference) > 0
    np.testing.assert_allclose(reference.sum(axis=1), 0, atol=1e-12)

    data = yaml.safe_load(force_input.read_text())
    if command == "run":
        data.update({"kappa": {"mesh": [2, 2, 2], "temps": [300]}, "plot": {}})
    else:
        with pytest.raises(ValueError, match="six sections"):
            parse_selected_command(force_input, "run")
    data["force-constant"]["parallel"] = {
        "backend": "slurm", "jobs": 2, "max-concurrent": 1, "submit": submit,
        "cpus-per-task": 2, "collect-cpus-per-task": 1,
    }
    force_input.write_text(yaml.safe_dump(data))
    selected, cfg = parse_selected_command(force_input, command)
    assert selected == command
    assert cfg.workflow_steps == (["relax", "fc2fc3", "kappa"] if command == "run" else ["fc2"])
    transport_calls = []

    def transport_mock(self):
        assert (self.output_dir / "fc2.hdf5").is_file()
        assert (self.output_dir / "fc3.hdf5").is_file()
        transport_calls.append(self.output_dir)

    monkeypatch.setattr(Phono3pyTransportWorkflow, "run", transport_mock)
    submissions = []

    def submit_mock(command, **kwargs):
        assert command[:2] == ["sbatch", "--parsable"]
        submissions.append(command)
        return SimpleNamespace(returncode=0, stdout=f"{100 + len(submissions)};test\n", stderr="")

    monkeypatch.setattr("nepkappa.slurm.SlurmScheduler", lambda: scheduler.SlurmScheduler(run=submit_mock))
    assert main([command, str(force_input)]) == 0
    assert transport_calls == []
    output = Path(cfg.result_dir)
    assert not (output / "fc2.hdf5").exists()
    work = output / "force-slurm"
    manifest = yaml.safe_load((work / "submission.yaml").read_text())
    assert manifest["status"] == ("submitted" if submit else "prepared")
    assert manifest["submitted"] is submit
    provenance = yaml.safe_load((output / "provenance.yaml").read_text())
    assert provenance["status"] == manifest["status"]
    log = (output / "run.log").read_text()
    assert f"[main] {'Submitted' if submit else 'Prepared'}: {command}" in log
    assert f"[main] Completed: {command}" not in log
    assert submissions == ([
        ["sbatch", "--parsable", str(work / "array.sh")],
        ["sbatch", "--parsable", "--dependency=afterok:101", str(work / "collect.sh")],
    ] if submit else [])
    array = (work / "array.sh").read_text()
    collect = (work / "collect.sh").read_text()
    assert "#SBATCH --array=0-1%1" in array
    assert "#SBATCH --cpus-per-task=2" in array
    assert "#SBATCH --cpus-per-task=1" in collect
    python = str(Path(sys.executable).resolve())
    assert shlex.join([python, "-m", "nepkappa.force_worker", "--config", str(force_input)]) in array
    assert ' --task-file "$TASK_FILE"' in array
    assert "NEPKAPPA_FORCE_COLLECT=1 " + shlex.join(
        [python, "-m", "nepkappa", command, str(force_input)]
    ) in collect

    collection_line = next(line for line in collect.splitlines() if line.startswith("NEPKAPPA_FORCE_COLLECT="))
    tokens = shlex.split(collection_line)
    assert tokens[:4] == ["NEPKAPPA_FORCE_COLLECT=1", python, "-m", "nepkappa"]
    cli_args = build_parser().parse_args(tokens[4:])
    parsed_command, parsed_config = parse_selected_command(cli_args.config, cli_args.command)
    assert parsed_command == command
    assert parsed_config.workflow_steps == cfg.workflow_steps

    tasks = sorted((work / "tasks").glob("*.json"))
    records = [record for task in tasks for record in json.loads(task.read_text())]
    assert len(records) == manifest["force_job_count"]
    assert len({(r["label"], r["index"]) for r in records}) == len(records)
    for task in tasks:
        run_task_file(force_input, task)
    cached = {Path(r["poscar"]).parent / "forces.npy" for r in records}
    before = {path: path.read_bytes() for path in cached}

    def unexpected_force(*args, **kwargs):
        pytest.fail("Collection must reuse the worker forces")

    monkeypatch.setattr(NEPPhononWorkflow, "_calculate_forces", unexpected_force)
    monkeypatch.setenv("NEPKAPPA_FORCE_COLLECT", "1")
    assert main(tokens[4:]) == 0
    assert transport_calls == ([output] if command == "run" else [])
    provenance = yaml.safe_load((output / "provenance.yaml").read_text())
    assert provenance["status"] == "complete"
    assert "Relax structure" not in provenance["stage_timings_seconds"]
    assert yaml.safe_load((work / "submission.yaml").read_text())["status"] == "complete"
    assert {path: path.read_bytes() for path in cached} == before
    with h5py.File(output / "fc2.hdf5") as handle:
        np.testing.assert_allclose(handle["force_constants"][:], reference, rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize("change", ["geometry", "calculator", "missing-marker", "old-marker", "nonfinite-cache"])
def test_force_cache_reuses_only_matching_valid_results(force_input, monkeypatch, change):
    cfg = parse_workflow_args(force_input, command="fc2")
    atoms = read(cfg.poscar)
    calls = []
    original = NEPPhononWorkflow._calculate_forces

    def calculate(self, structure, label, index):
        calls.append((label, index))
        return original(self, structure, label, index)

    monkeypatch.setattr(NEPPhononWorkflow, "_calculate_forces", calculate)
    core = NEPPhononWorkflow(cfg)
    first = core._calculate_forces_cached(atoms.copy(), "fc2", 1)
    directory = core._force_job_dir("fc2", 1)
    marker = directory / "job-input.sha256"
    first_digest = marker.read_text()
    np.testing.assert_array_equal(
        NEPPhononWorkflow(cfg)._calculate_forces_cached(atoms.copy(), "fc2", 1), first
    )
    assert len(calls) == 1
    if change == "geometry":
        atoms.positions[1, 0] += 0.1
    elif change == "calculator":
        cfg.calculator_kwargs["epsilon"] *= 2
    elif change == "missing-marker":
        marker.unlink()
    elif change == "old-marker":
        marker.write_text("unverified-old-fingerprint\n")
    else:
        np.save(directory / "forces.npy", np.full((2, 3), np.nan))

    fresh = NEPPhononWorkflow(cfg)
    result = fresh._calculate_forces_cached(atoms.copy(), "fc2", 1)
    assert len(calls) == 2
    assert np.isfinite(result).all()
    np.testing.assert_allclose(result.sum(axis=0), 0, atol=1e-14)
    if change == "calculator":
        np.testing.assert_allclose(result, 2 * first)
    elif change == "geometry":
        assert not np.allclose(result, first)
    else:
        np.testing.assert_array_equal(result, first)
    assert (marker.read_text() != first_digest) == (change in {"geometry", "calculator"})
    record = json.loads((directory / "job-input.json").read_text())
    assert record["sha256"] == marker.read_text().strip()


def test_failed_force_recalculation_does_not_certify_stale_cache(force_input, monkeypatch):
    cfg = parse_workflow_args(force_input, command="fc2")
    core = NEPPhononWorkflow(cfg)
    atoms = read(cfg.poscar)
    core._calculate_forces_cached(atoms.copy(), "fc2", 1)
    directory = core._force_job_dir("fc2", 1)
    old = {name: (directory / name).read_bytes() for name in
           ("forces.npy", "job-input.sha256", "job-input.json")}
    cfg.calculator_kwargs["epsilon"] *= 2
    monkeypatch.setattr(NEPPhononWorkflow, "_calculate_forces", lambda *args: np.full((2, 3), np.nan))
    with pytest.raises(RuntimeError, match="Non-finite"):
        NEPPhononWorkflow(cfg)._calculate_forces_cached(atoms.copy(), "fc2", 1)
    assert {name: (directory / name).read_bytes() for name in old} == old


@pytest.mark.parametrize("interrupted_file", ["forces.npy", "job-input.json", "job-input.sha256"])
def test_interrupted_cache_update_cannot_certify_new_forces_with_old_inputs(
    force_input, monkeypatch, interrupted_file
):
    cfg = parse_workflow_args(force_input, command="fc2")
    atoms = read(cfg.poscar)
    original = NEPPhononWorkflow(cfg)
    first = original._calculate_forces_cached(atoms.copy(), "fc2", 1)
    directory = original._force_job_dir("fc2", 1)
    marker = directory / "job-input.sha256"
    old_digest = marker.read_text()
    cfg.calculator_kwargs["epsilon"] *= 2
    replace = Path.replace

    def interrupt_after_replace(source, destination):
        result = replace(source, destination)
        if Path(destination) == directory / interrupted_file:
            raise KeyboardInterrupt("simulated process interruption")
        return result

    with monkeypatch.context() as patch:
        patch.setattr(Path, "replace", interrupt_after_replace)
        with pytest.raises(KeyboardInterrupt, match="simulated"):
            NEPPhononWorkflow(cfg)._calculate_forces_cached(atoms.copy(), "fc2", 1)
    np.testing.assert_allclose(np.load(directory / "forces.npy"), 2 * first)
    if interrupted_file == "job-input.sha256":
        assert marker.read_text() != old_digest  # A fully published new identity.
    else:
        assert not marker.exists()  # No stale identity during publication.

    # Returning to the old model used to accept the new model's force array.
    cfg.calculator_kwargs["epsilon"] /= 2
    calls = []
    calculate = NEPPhononWorkflow._calculate_forces

    def counted(self, *args):
        calls.append(True)
        return calculate(self, *args)

    monkeypatch.setattr(NEPPhononWorkflow, "_calculate_forces", counted)
    recovered = NEPPhononWorkflow(cfg)._calculate_forces_cached(atoms.copy(), "fc2", 1)
    assert len(calls) == 1
    np.testing.assert_array_equal(recovered, first)
    assert marker.read_text() == old_digest
    record = json.loads((directory / "job-input.json").read_text())
    assert record["sha256"] == old_digest.strip()
    np.testing.assert_array_equal(
        NEPPhononWorkflow(cfg)._calculate_forces_cached(atoms.copy(), "fc2", 1), first
    )
    assert len(calls) == 1


@pytest.mark.parametrize("backend", ["thirdorder", "fourthorder"])
def test_interrupted_order_force_job_invalidates_old_xml_identity(force_input, monkeypatch, backend):
    from nepkappa.stages.force_constants import _calculate_order_forces

    cfg = parse_workflow_args(force_input, command="fc2")
    core = NEPPhononWorkflow(cfg)
    # These commands are identified but never executed by the XML force path.
    cfg.thirdorder_command = "synthetic-thirdorder"
    cfg.fourthorder_command = "synthetic-fourthorder"
    cfg.cutoff_fc3 = -3.0
    pattern = Path(cfg.poscar)
    paths = _calculate_order_forces(core, [pattern], backend=backend, label="test")
    marker = paths[0].parent / "job-input.sha256"
    old_digest = marker.read_text()
    old_xml = paths[0].read_bytes()
    cfg.calculator_kwargs["epsilon"] *= 2
    changed = NEPPhononWorkflow(cfg)
    write_xml = changed._run_external_order_job

    def interrupted(*args):
        write_xml(*args)
        raise KeyboardInterrupt("interrupted XML force job")

    monkeypatch.setattr(changed, "_run_external_order_job", interrupted)
    with pytest.raises(KeyboardInterrupt, match="interrupted XML"):
        _calculate_order_forces(changed, [pattern], backend=backend, label="test")
    assert not marker.exists()
    assert paths[0].read_bytes() != old_xml
    cfg.calculator_kwargs["epsilon"] /= 2
    _calculate_order_forces(NEPPhononWorkflow(cfg), [pattern], backend=backend, label="test")
    assert marker.read_text() == old_digest
    assert paths[0].read_bytes() == old_xml
