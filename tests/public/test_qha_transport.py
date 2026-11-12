"""QHA-coupled RTA regenerates force constants at each expanded volume."""

from pathlib import Path
from types import SimpleNamespace

from ase import Atoms
import numpy as np
import pytest
import yaml

from nepkappa.config import parse_workflow_args
from nepkappa.application import WorkflowStageRunner
from nepkappa.qha_transport import QHA3PhWorkflow


def test_qha_transport_example_selects_coupled_stages():
    config = parse_workflow_args("examples/nep-qha-rta-3ph.yaml", command="run")
    assert config.workflow_steps == ["qha", "qha-kappa"]
    assert config.kappa_engine == "phono3py"
    assert config.qha_volumes is True


def test_qha_transport_uses_interpolated_volumes(monkeypatch, tmp_path):
    config = parse_workflow_args("examples/nep-qha-rta-3ph.yaml", command="run")
    config.result_dir = str(tmp_path)
    config.temps = [300, 900, 300]
    np.savetxt(tmp_path / "volume-temperature.dat", [[0, 20], [900, 23]])
    (tmp_path / "qha-summary.yaml").write_text(
        "normalization_multiplicity: 1\n", encoding="utf-8"
    )
    workflow = QHA3PhWorkflow(config, execution=SimpleNamespace())
    monkeypatch.setattr(
        workflow, "_run_temperature",
        lambda base, temperature, volume, multiplicity: {
            "temperature_K": temperature,
            "qha_primitive_volume_A3": volume,
        },
    )
    result = workflow.run()
    assert [item["qha_primitive_volume_A3"] for item in result["cases"]] == [21, 22, 23]
    assert yaml.safe_load((tmp_path / "qha-kappa/qha-kappa-summary.yaml").read_text()) == result


def test_qha_transport_runs_fc_then_kappa_at_one_temperature(monkeypatch, tmp_path):
    calls = []
    cases = []

    class FakeRunner:
        def __init__(self, config, execution):
            cases.append(config)

        def run(self, stage):
            calls.append(stage)
            return None

    monkeypatch.setattr("nepkappa.application.WorkflowStageRunner", FakeRunner)
    config = parse_workflow_args("examples/nep-qha-rta-3ph.yaml", command="run")
    config.result_dir = str(tmp_path)
    config.qha_relax_internal = False
    workflow = QHA3PhWorkflow(config, execution=SimpleNamespace(), runner_factory=lambda case: FakeRunner(case, None))
    base = Atoms("Si", cell=np.eye(3) * 5, scaled_positions=[[0, 0, 0]], pbc=True)
    result = workflow._run_temperature(base, 300.0, 150.0, 1)
    assert calls == ["fc2fc3", "kappa"]
    assert cases[0].temps == [300.0]
    assert Path(cases[0].poscar).is_file()
    assert result["qha_primitive_volume_A3"] == 150.0


def test_engine_and_switch_must_agree(tmp_path):
    data = yaml.safe_load(Path("examples/nep-rta-3ph-4ph.yaml").read_text())
    data["kappa"]["engine"] = "phono3py"
    path = tmp_path / "wrong.yaml"
    path.write_text(yaml.safe_dump(data), encoding="utf-8")
    with pytest.raises(ValueError, match="requires kappa.engine: fourphonon"):
        parse_workflow_args(path, command="run")


def test_parallel_is_separate_from_force_constant(tmp_path):
    data = yaml.safe_load(Path("examples/nep-rta-3ph-4ph.yaml").read_text())
    data["parallel"] = {
        "force-constant": {"backend": "none"},
        "kappa": {"backend": "none", "mpi-processes": 2, "omp-threads": 4},
    }
    path = tmp_path / "parallel.yaml"
    path.write_text(yaml.safe_dump(data), encoding="utf-8")
    config = parse_workflow_args(path, command="run")
    assert config.force_parallel["backend"] == "none"
    assert config.fp_parallel["backend"] == "none"
    assert config.fp_mpi_processes == 2
    assert config.fp_omp_threads == 4


def test_qha_transport_plot_targets_selected_temperature(monkeypatch, tmp_path):
    config = parse_workflow_args("examples/nep-qha-rta-3ph.yaml", command="plot")
    config.result_dir = str(tmp_path)
    config.plot_temperature = 600.0
    observed = []
    monkeypatch.setattr(
        "nepkappa.plot.plot_results",
        lambda selected: observed.append(selected.result_dir),
    )
    execution = SimpleNamespace(run_stage=lambda label, operation: operation())
    WorkflowStageRunner(config, execution).plot()
    assert observed == [str(tmp_path / "qha-kappa" / "T0600K")]
