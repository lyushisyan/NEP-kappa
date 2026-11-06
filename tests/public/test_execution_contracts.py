"""Execution outcomes and shared input helpers, without external solvers."""

from io import StringIO
from pathlib import Path

import numpy as np
import pytest
import yaml

from nepkappa.application import execute_workflow_command, execute_workflow_plan, WorkflowStageRunner
from nepkappa.config import parse_selected_command


@pytest.mark.parametrize("backend", ["phono3py", "fourphonon"])
@pytest.mark.parametrize("status", ["prepared", "submitted", "complete"])
def test_transport_command_and_plan_share_execution_status(tmp_path, monkeypatch, backend, status):
    from nepkappa.cli import main
    from nepkappa.fourphonon import FourPhononWorkflow
    from nepkappa.transport import Phono3pyTransportWorkflow

    root = Path(__file__).resolve().parents[2]
    name = "nep-rta-3ph-4ph.yaml" if backend == "fourphonon" else "nep-lbte-wigner-3ph.yaml"
    data = yaml.safe_load((root / "examples" / name).read_text())
    data["output"]["result-dir"] = str(tmp_path / "result")
    source = tmp_path / "input.yaml"
    source.write_text(yaml.safe_dump(data))
    command, cfg = parse_selected_command(source, "kappa")
    manifest = {"status": status, "submitted": status == "submitted"}
    backend_class = FourPhononWorkflow if backend == "fourphonon" else Phono3pyTransportWorkflow
    monkeypatch.setattr(backend_class, "run", lambda self: manifest)

    assert execute_workflow_command(command, cfg).status == status
    reports = []
    monkeypatch.setattr(WorkflowStageRunner, "report", lambda self: reports.append(True))
    cfg.workflow_steps = [command, "report"]
    assert execute_workflow_plan(cfg).status == status
    assert reports == ([True] if status == "complete" else [])
    assert main(["kappa", str(source)]) == 0
    output = Path(cfg.result_dir)
    assert yaml.safe_load((output / "provenance.yaml").read_text())["status"] == status
    label = {"prepared": "Prepared", "submitted": "Submitted", "complete": "Completed"}[status]
    assert f"[main] {label}: kappa" in (output / "run.log").read_text()


@pytest.mark.parametrize("text", [
    "output: {}\noutput: {}\n",
    "output:\n  result-dir: first\n  result-dir: second\n",
])
def test_strict_yaml_duplicates_rejected_by_config_and_report(tmp_path, text):
    from nepkappa.config import parse_workflow_args
    from nepkappa.config_schema import load_yaml_strict
    from nepkappa.report import resolve_result_directory

    path = tmp_path / "duplicate.yaml"
    path.write_text(text)
    with pytest.raises(ValueError, match="Duplicate YAML key"):
        load_yaml_strict(StringIO(text))
    with pytest.raises(ValueError, match="Duplicate YAML key"):
        parse_workflow_args(path, command="report")
    with pytest.raises(ValueError, match="Duplicate YAML key"):
        resolve_result_directory(path)


def test_shared_qha_helpers_keep_geometry_and_general_errors(tmp_path):
    from ase import Atoms
    from nepkappa.qha import load_qha_equilibrium_volumes, scale_atoms_to_primitive_volume, validate_temperature_range

    with pytest.raises(FileNotFoundError, match="Completed QHA outputs are required"):
        load_qha_equilibrium_volumes(tmp_path)
    (tmp_path / "volume-temperature.dat").write_text("0 20\n300 22\n600 24\n")
    (tmp_path / "qha-summary.yaml").write_text("normalization_multiplicity: 2\n")
    temperatures, volumes, multiplicity = load_qha_equilibrium_volumes(tmp_path)
    validate_temperature_range(np.array([150, 450]), temperatures)
    atoms = Atoms("Si2", scaled_positions=[[0, 0, 0], [.25, .25, .25]], cell=[4]*3, pbc=True)
    scaled, _ = scale_atoms_to_primitive_volume(atoms, volumes[1], multiplicity)
    assert scaled.get_volume() == pytest.approx(44)
    np.testing.assert_allclose(scaled.get_scaled_positions(), atoms.get_scaled_positions())
    with pytest.raises(ValueError, match="No temperatures were requested"):
        validate_temperature_range(np.array([]), temperatures)
    with pytest.raises(ValueError, match="Requested temperatures must lie inside"):
        validate_temperature_range(np.array([700]), temperatures)
