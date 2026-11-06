"""Synchronous SSCHA policy and real temperature loops with mocked numerics."""

from pathlib import Path

from ase import Atoms
import pytest
import yaml

from nepkappa.application import CommandExecution, WorkflowStageRunner
from nepkappa.cli import main
from nepkappa.config import parse_selected_command
from nepkappa.qha_sscha import QHASSCHAWorkflow
from nepkappa.sscha import PhonopySSCHAWorkflow, temperature_directory_name
from nepkappa.transport import Phono3pyTransportWorkflow
from nepkappa.workflow import NEPPhononWorkflow


@pytest.fixture(params=[False, True], ids=["sscha", "qha-sscha"])
def sscha_input(tmp_path, request):
    structure = tmp_path / "POSCAR"
    Atoms("Ar", cell=[5]*3, pbc=True).write(structure, format="vasp")
    data = {
        "structure": {"poscar": str(structure), "dimensionality": 3, "relaxation": False},
        "calculator": {"name": "ase", "factory": "ase.calculators.lj:LennardJones"},
        "force-constant": {
            "dim-fc2": [1, 1, 1], "dim-fc3": [1, 1, 1],
            "sscha": {"enabled": True, "run-transport": True, "temps": [300, 900, 300],
                      "snapshots": 2, "iterations": 2, "transient": 1},
        },
        "kappa": {"engine": "phono3py", "method": "rta", "mesh": [2, 2, 2],
                  "temps": [300, 900, 300]},
        "plot": {},
        "output": {"result-dir": str(tmp_path / "result"), "progress": False},
    }
    if request.param:
        data["force-constant"]["qha"] = {"enabled": True, "relax-internal": False}
    path = tmp_path / "input.yaml"
    path.write_text(yaml.safe_dump(data))
    return path


@pytest.mark.parametrize("submit", [False, True])
def test_nested_slurm_rejected_by_cli_before_execution(sscha_input, monkeypatch, capsys, submit):
    data = yaml.safe_load(sscha_input.read_text())
    data["kappa"].update({"method": "lbte", "parallel": {"backend": "slurm", "submit": submit}})
    sscha_input.write_text(yaml.safe_dump(data))

    def unexpected(*args, **kwargs):
        pytest.fail("Invalid nested transport must be rejected before application execution")

    monkeypatch.setattr("nepkappa.application.execute_workflow_command", unexpected)
    with pytest.raises(SystemExit) as exc:
        main(["run", str(sscha_input)])
    assert exc.value.code == 2
    assert "SSCHA and QHA-SSCHA require synchronous transport" in capsys.readouterr().err
    assert not Path(data["output"]["result-dir"]).exists()


@pytest.mark.parametrize("submit", [False, True])
def test_direct_workflow_calls_reject_nested_slurm_before_science(sscha_input, monkeypatch, submit):
    _, cfg = parse_selected_command(sscha_input, "run")
    # Programmatic callers can bypass the parser; run() must still reject early.
    cfg.method = "lbte"
    cfg.lbte_parallel = {"backend": "slurm", "submit": submit}
    execution = CommandExecution()

    def unexpected(*args, **kwargs):
        pytest.fail("Nested Slurm was not rejected before scientific work")

    monkeypatch.setattr(NEPPhononWorkflow, "load_force_constant_structure", unexpected)
    monkeypatch.setattr("nepkappa.qha_sscha.load_qha_equilibrium_volumes", unexpected)
    monkeypatch.setattr(PhonopySSCHAWorkflow, "_run_temperature", unexpected)
    monkeypatch.setattr(Phono3pyTransportWorkflow, "run", unexpected)
    if cfg.qha_enabled:
        workflow = QHASSCHAWorkflow(cfg, execution=execution)
    else:
        workflow = PhonopySSCHAWorkflow(cfg, NEPPhononWorkflow(cfg), execution=execution)
    with pytest.raises(ValueError, match="SSCHA and QHA-SSCHA require synchronous transport"):
        workflow.run()
    assert not workflow.workdir.exists()
    assert execution.stage_timings == []


@pytest.mark.parametrize("method,status", [
    ("rta", "complete"), ("lbte", "complete"),
    ("lbte", "prepared"), ("lbte", "submitted"),
])
def test_real_sscha_loops_never_complete_deferred_transport(
    sscha_input, monkeypatch, capsys, method, status
):
    from nepkappa.qha import QHAWorkflow

    data = yaml.safe_load(sscha_input.read_text())
    data["kappa"].update({"method": method, "parallel": {"backend": "none"}})
    sscha_input.write_text(yaml.safe_dump(data))
    _, cfg = parse_selected_command(sscha_input, "run")
    temperatures, transport_calls = [], []

    def qha_fit(self):
        (self.output_dir / "volume-temperature.dat").write_text("0 125\n900 128\n")
        (self.output_dir / "qha-summary.yaml").write_text("normalization_multiplicity: 1\n")

    def force_constants(self, *, include_fc3):
        assert include_fc3
        for name in ("fc2.hdf5", "fc3.hdf5", "phono3py_disp.yaml"):
            (self.core.output_dir / name).write_text("synthetic input; numerical operations are mocked")

    def temperature_calculation(self, temperature, initial_fc2, transport_temperatures):
        assert initial_fc2.is_file()
        assert temperature in transport_temperatures
        temperatures.append(temperature)
        label = temperature_directory_name(temperature)
        directory = self.workdir / label
        directory.mkdir(parents=True)
        # Keep the real nested transport dispatch and both real temperature loops.
        self._run_transport(temperature, directory)
        return {"temperature": temperature, "directory": label, "bubble": None}

    def transport(self):
        transport_calls.append(self.config.temps[0])
        # Fault injection: even a serial backend returning deferred work must fail.
        return {"status": status, "submitted": status == "submitted"}

    monkeypatch.setattr(QHAWorkflow, "run", qha_fit)
    monkeypatch.setattr(WorkflowStageRunner, "force_constants", force_constants)
    monkeypatch.setattr(PhonopySSCHAWorkflow, "_run_temperature", temperature_calculation)
    monkeypatch.setattr(Phono3pyTransportWorkflow, "run", transport)
    result = main(["run", str(sscha_input)])
    output = Path(cfg.result_dir)
    provenance = yaml.safe_load((output / "provenance.yaml").read_text())
    terminal = capsys.readouterr().out
    if status == "complete":
        assert result == 0
        assert temperatures == transport_calls == [300, 600, 900]
        assert provenance["status"] == "complete"
        assert len(list(output.rglob("sscha-summary.yaml"))) == (3 if cfg.qha_enabled else 1)
        if cfg.qha_enabled:
            summary = yaml.safe_load((output / "qha-sscha/qha-sscha-summary.yaml").read_text())
            assert summary["status"] == "complete" and len(summary["cases"]) == 3
    else:
        assert result == 1
        assert temperatures == transport_calls == [300]
        assert provenance["status"] == "failed"
        assert "SSCHA transport must complete synchronously" in provenance["error"]
        assert not list(output.rglob("sscha-summary.yaml"))
        assert not (output / "qha-sscha/qha-sscha-summary.yaml").exists()
        assert "[Done] Phonopy SSCHA" not in terminal
        assert "[main] Completed:" not in terminal


@pytest.mark.parametrize("submit", [False, True])
def test_direct_slurm_lbte_remains_supported(tmp_path, submit):
    source = tmp_path / "kappa.yaml"
    source.write_text(yaml.safe_dump({
        "kappa": {"method": "lbte", "mesh": [2, 2, 2], "temps": [300],
                  "parallel": {"backend": "slurm", "submit": submit}},
        "output": {"result-dir": str(tmp_path / "result")},
    }))
    command, cfg = parse_selected_command(source, "kappa")
    assert command == "kappa"
    assert cfg.lbte_parallel == {"backend": "slurm", "submit": submit}
