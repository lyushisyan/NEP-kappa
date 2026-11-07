"""Distributed RTA contracts and real serial/worker/collection equivalence."""

import os
from pathlib import Path
import subprocess

import h5py
import numpy as np
import pytest
import yaml

from nepkappa.application import execute_workflow_command
from nepkappa.cli import main
from nepkappa.config import parse_workflow_args
from nepkappa.qha_transport import QHA3PhWorkflow
from nepkappa.report import _read_transport
from nepkappa.transport import Phono3pyTransportWorkflow


def transport_input(tmp_path, *, wigner=False, submit=False):
    data = {
        "kappa": {"engine": "phono3py", "method": "rta", "mesh": [3, 3, 3],
                  "temps": [200, 400, 100], "isotope": True, "bfmp": 10,
                  "wigner": wigner},
        "parallel": {"kappa": {
            "backend": "slurm", "jobs": 2, "max-concurrent": 1,
            "cpus-per-task": 2, "collect-cpus-per-task": 1, "submit": submit,
        }},
        "output": {"result-dir": str(tmp_path / "result with spaces"), "progress": False},
    }
    path = tmp_path / "input.yaml"
    path.write_text(yaml.safe_dump(data))
    return path


@pytest.mark.parametrize("submit", [False, True])
@pytest.mark.parametrize("wigner", [False, True])
@pytest.mark.parametrize("temperatures", [[300], [200, 400, 100]])
def test_rta_dispatch_flags_submission_and_deferred_status(
    tmp_path, monkeypatch, submit, wigner, temperatures
):
    path = transport_input(tmp_path, submit=submit, wigner=wigner)
    data = yaml.safe_load(path.read_text())
    data["kappa"]["temps"] = temperatures
    path.write_text(yaml.safe_dump(data))
    cfg = parse_workflow_args(path, command="kappa")
    output = Path(cfg.result_dir)
    output.mkdir()
    for name in ("fc2.hdf5", "fc3.hdf5", "phono3py_disp.yaml", "BORN"):
        (output / name).write_text("synthetic input; numerical calls mocked")
    commands, submissions = [], []

    def run(command, cwd=None):
        commands.append(command)
        assert "--wgp" in command and "--br" in command
        (Path(cwd) / "ir_grid_points.yaml").write_text(yaml.safe_dump({
            "ir_grid_points": [{"grid_point": gp} for gp in [0, 3, 7, 9, 12]],
        }))
        return 0

    def submit_job(self, script, cwd, dependency=None):
        submissions.append((script.name, dependency))
        return str(100 + len(submissions))

    monkeypatch.setattr("nepkappa.workflow.NEPPhononWorkflow._run_command", staticmethod(run))
    monkeypatch.setattr("nepkappa.slurm.SlurmScheduler.submit", submit_job)
    assert main(["kappa", str(path)]) == 0
    work = output / "rta-slurm"
    manifest = yaml.safe_load((work / "submission.yaml").read_text())
    assert manifest["workflow"] == "rta"
    assert manifest["status"] == ("submitted" if submit else "prepared")
    assert manifest["array_jobs"] == 2
    assert submissions == ([("prepare.sh", None), ("array.sh", "101"), ("collect.sh", "102")] if submit else [])
    assert yaml.safe_load((output / "provenance.yaml").read_text())["status"] == manifest["status"]
    calculation = Path(manifest["calculation_directory"])
    assert calculation.parent == work
    assert (calculation / "BORN").resolve() == output / "BORN"
    assert (work / "grid-points/000.txt").read_text().strip() == "0,7,12"
    assert (work / "grid-points/001.txt").read_text().strip() == "3,9"
    worker = (work / "array.sh").read_text()
    collect = (work / "collect.sh").read_text()
    assert "--array=0-1%1" in worker
    assert "--write-gamma" in worker and '--gp "$GRID_POINTS"' in worker
    assert "--read-gamma" in collect and "--gp" not in collect
    assert "--write-pp" not in worker and "--read-pp" not in collect
    for script in (worker, collect):
        assert "--br --nu" in script
        assert ("--ts 300.0" if len(temperatures) == 1 else
                "--tmin 200.0 --tmax 400.0 --tstep 100.0") in script
        assert "--isotope" in script and "--boundary-mfp 10.0" in script
        assert ("--tt smm19" in script) == wigner
        assert "--read-phonon" in script
    assert "#SBATCH --cpus-per-task=2" in worker
    assert "export RAYON_NUM_THREADS=${SLURM_CPUS_PER_TASK:-2}" in worker
    assert "export RAYON_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}" in collect
    assert collect.index("cp ") < collect.index("nepkappa.run_state --complete")
    for name in ("prepare.sh", "array.sh", "collect.sh"):
        result = subprocess.run(["bash", "-n", str(work / name)], capture_output=True, text=True)
        assert result.returncode == 0, result.stderr

    # Preparing again creates a fresh gamma workspace, never reuses a full result.
    before = set(work.glob("run-*"))
    assert main(["kappa", str(path)]) == 0
    assert len(set(work.glob("run-*")) - before) == 1


@pytest.mark.parametrize("submit", [False, True])
@pytest.mark.parametrize("coupled", ["qha", "sscha", "qha-sscha"])
def test_nested_rta_slurm_rejected_before_temperature_loops(tmp_path, submit, coupled):
    root = Path(__file__).resolve().parents[2]
    name = "nep-qha-rta-3ph.yaml" if coupled == "qha" else "nep-qha-sscha-rta-3ph.yaml"
    data = yaml.safe_load((root / "examples" / name).read_text())
    if coupled == "sscha":
        data["force-constant"].pop("qha")
    data["parallel"] = {"kappa": {"backend": "slurm", "submit": submit}}
    data["output"]["result-dir"] = str(tmp_path / "result")
    path = tmp_path / "input.yaml"
    path.write_text(yaml.safe_dump(data))
    with pytest.raises(SystemExit):
        parse_workflow_args(path, command="run")
    assert not (tmp_path / "result").exists()


def test_programmatic_qha_rejects_nested_rta_before_reading_outputs(tmp_path):
    cfg = parse_workflow_args("examples/nep-qha-rta-3ph.yaml", command="run")
    cfg.result_dir = str(tmp_path)
    cfg.lbte_parallel = {"backend": "slurm", "submit": False}
    with pytest.raises(ValueError, match="requires synchronous transport"):
        QHA3PhWorkflow(cfg, execution=None).run()
    assert not (tmp_path / "qha-kappa").exists()


@pytest.mark.parametrize("status", ["prepared", "submitted"])
def test_qha_never_accepts_unexpected_deferred_rta(tmp_path, monkeypatch, status):
    from ase import Atoms

    class Runner:
        def __init__(self, config, execution):
            pass

        def run(self, stage):
            return {"status": status} if stage == "kappa" else None

    monkeypatch.setattr("nepkappa.application.WorkflowStageRunner", Runner)
    cfg = parse_workflow_args("examples/nep-qha-rta-3ph.yaml", command="run")
    cfg.result_dir = str(tmp_path)
    cfg.qha_relax_internal = False
    base = Atoms("Si", cell=[5, 5, 5], pbc=True)
    with pytest.raises(RuntimeError, match="cannot continue from a deferred kappa job"):
        QHA3PhWorkflow(cfg, execution=None)._run_temperature(base, 300, 125, 1)


@pytest.fixture(scope="module")
def real_force_constants(tmp_path_factory):
    """Small physical NEP model, reused by the serial/distributed comparisons."""
    root = Path(__file__).resolve().parents[2]
    directory = tmp_path_factory.mktemp("rta-fc")
    path = directory / "forces.yaml"
    path.write_text(yaml.safe_dump({
        "structure": {"poscar": str(root / "examples/structures/3C-SiC/POSCAR_primitive"),
                      "relaxation": False},
        "calculator": {"name": "nep", "model": str(root / "potentials/3C-SiC/nep_3C-SiC.txt")},
        "force-constant": {"dim-fc2": [2, 2, 2], "dim-fc3": [2, 2, 2]},
        "output": {"result-dir": str(directory / "forces"), "progress": False},
    }))
    cfg = parse_workflow_args(path, command="fc2fc3")
    assert execute_workflow_command("fc2fc3", cfg).status == "complete"
    return Path(cfg.result_dir)


@pytest.mark.parametrize("wigner", [False, True])
def test_real_rta_workers_and_collection_match_serial(tmp_path, monkeypatch, real_force_constants, wigner):
    for key in ("RAYON_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        monkeypatch.setenv(key, "1")
    path = transport_input(tmp_path, wigner=wigner)
    cfg = parse_workflow_args(path, command="kappa")
    cfg.kappa_force_constants = str(real_force_constants)
    output = Path(cfg.result_dir)
    calls = []

    def run(command, cwd=None):
        calls.append(command)
        result = subprocess.run(command, cwd=cwd, capture_output=True, text=True, timeout=60)
        assert result.returncode == 0, result.stdout + result.stderr
        return result.returncode

    # Serial result in the same output directory deliberately becomes stale.
    cfg.lbte_parallel = {}
    serial = Phono3pyTransportWorkflow(cfg, output, run_command=run)
    serial.run()
    filename = "kappa-m333.hdf5"
    with h5py.File(output / filename) as handle:
        names = ["temperature", "gamma", "gamma_N", "gamma_U", "gamma_isotope"]
        names += ["kappa", "kappa_intra", "kappa_inter"] if wigner else ["kappa"]
        reference = {name: handle[name][:] for name in names}
        assert np.isfinite(reference["gamma"]).all()
        assert np.max(reference["gamma"]) > 0
        assert np.max(reference["kappa"]) > 0
    with h5py.File(output / filename, "r+") as handle:
        handle["gamma"][...] = -999  # Must not be read by distributed collection.

    cfg.lbte_parallel = {"backend": "slurm", "jobs": 2, "cpus_per_task": 1, "submit": False}
    manifest = Phono3pyTransportWorkflow(cfg, output, run_command=run).run()
    assert manifest["status"] == "prepared"
    work = output / "rta-slurm"

    def script(name, **environment):
        env = {**os.environ, **environment}
        result = subprocess.run(["bash", str(work / name)], cwd=output, env=env,
                                capture_output=True, text=True, timeout=60)
        assert result.returncode == 0, result.stdout + result.stderr
        return result

    script("prepare.sh")
    for index in range(manifest["array_jobs"]):
        script("array.sh", SLURM_ARRAY_TASK_ID=str(index))
    collected = script("collect.sh")
    assert "Read data from kappa-m333-g" in collected.stdout
    assert yaml.safe_load((work / "submission.yaml").read_text())["status"] == "complete"
    with h5py.File(output / filename) as handle:
        for name, expected in reference.items():
            assert np.isfinite(handle[name][:]).all()
            np.testing.assert_allclose(handle[name][:], expected, rtol=1e-10, atol=1e-10)
    warnings = []
    records = _read_transport(output, warnings)
    assert not warnings
    assert len(records) == (9 if wigner else 3)
    assert {record["file"] for record in records} == {filename}

    # Missing worker gamma must fail, without publishing output or certifying completion.
    (Path(manifest["calculation_directory"]) / filename).unlink()
    next(Path(manifest["calculation_directory"]).glob("kappa-m333-g*.hdf5")).unlink()
    state = work / "submission.yaml"
    state.write_text(yaml.safe_dump({**manifest, "status": "prepared"}))
    before = (output / filename).read_bytes()
    result = subprocess.run(["bash", str(work / "collect.sh")], cwd=output,
                            capture_output=True, text=True, timeout=60)
    assert result.returncode != 0
    assert (output / filename).read_bytes() == before
    assert yaml.safe_load(state.read_text())["status"] == "prepared"
