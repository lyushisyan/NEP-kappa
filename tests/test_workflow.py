from pathlib import Path
from types import SimpleNamespace
import subprocess
import sys

import h5py
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nepkappa.workflow import NEPPhononWorkflow


class FakePhono3pyForFC:
    def __init__(self):
        self.supercells_with_displacements = []
        self.phonon_supercells_with_displacements = []
        self.phonon_primitive = SimpleNamespace(p2s_map=np.array([0], dtype="int64"))
        self.primitive = SimpleNamespace(p2s_map=np.array([0], dtype="int64"))
        self.fc3_nonzero_indices = None
        self.produce_fc2_compact = None
        self.produce_fc3_compact = None
        self.generated_fc2_only = False
        self.generated_fc2fc3 = False
        self.events = []
        self.fc2 = None
        self.fc3 = None

    def generate_displacements(self):
        self.events.append("generate_fc3_displacements")
        self.generated_fc2fc3 = True

    def generate_fc2_displacements(self):
        self.events.append("generate_fc2_displacements")
        self.generated_fc2_only = True

    def save(self, filename):
        self.events.append("save")
        Path(filename).write_text("disp", encoding="utf-8")

    def produce_fc2(self, is_compact_fc=True):
        self.events.append("produce_fc2")
        self.produce_fc2_compact = is_compact_fc
        first_dim = 1 if is_compact_fc else 2
        self.fc2 = np.zeros((first_dim, 2, 3, 3), dtype="double")

    def symmetrize_fc2(self):
        self.events.append("symmetrize_fc2")
        pass

    def produce_fc3(self, is_compact_fc=True):
        self.events.append("produce_fc3")
        self.produce_fc3_compact = is_compact_fc
        first_dim = 1 if is_compact_fc else 2
        self.fc3 = np.zeros((first_dim, 2, 2, 3, 3, 3), dtype="double")

    def symmetrize_fc3(self):
        self.events.append("symmetrize_fc3")
        pass


def test_compute_kappa_raises_when_phono3py_fails(tmp_path):
    cfg = SimpleNamespace(
        result_dir=str(tmp_path),
        mesh=[1, 1, 1],
        method="rta",
        wigner=False,
        isotope=False,
        bfmp=1.0e6,
        temps=[300],
        progress=False,
    )
    workflow = NEPPhononWorkflow(cfg)
    workflow.fc2_path.write_text("fc2", encoding="utf-8")
    workflow.fc3_path.write_text("fc3", encoding="utf-8")
    workflow.disp_path.write_text("disp", encoding="utf-8")
    workflow._run_command = lambda cmd, cwd=None: 2

    with pytest.raises(RuntimeError, match="Phono3py failed with return code 2"):
        workflow.compute_kappa()


def test_finite_displacement_defaults_to_compact_fc(tmp_path):
    cfg = SimpleNamespace(
        result_dir=str(tmp_path),
        progress=False,
        compact_fc=True,
    )
    workflow = NEPPhononWorkflow(cfg)
    ph3 = FakePhono3pyForFC()
    workflow._make_phono3py = lambda: ph3

    workflow.run_finite_disp_fitting()

    assert ph3.events == [
        "generate_fc2_displacements",
        "save",
        "produce_fc2",
        "symmetrize_fc2",
        "generate_fc3_displacements",
        "save",
        "produce_fc3",
        "symmetrize_fc3",
    ]
    assert ph3.produce_fc2_compact is True
    assert ph3.produce_fc3_compact is True
    with h5py.File(workflow.fc2_path, "r") as handle:
        assert handle["force_constants"].shape == (1, 2, 3, 3)
        assert "p2s_map" in handle
    with h5py.File(workflow.fc3_path, "r") as handle:
        assert handle["fc3"].shape == (1, 2, 2, 3, 3, 3)
        assert "p2s_map" in handle


def test_finite_displacement_can_write_full_fc(tmp_path):
    cfg = SimpleNamespace(
        result_dir=str(tmp_path),
        progress=False,
        compact_fc=False,
    )
    workflow = NEPPhononWorkflow(cfg)
    ph3 = FakePhono3pyForFC()
    workflow._make_phono3py = lambda: ph3

    workflow.run_finite_disp_fitting()

    assert ph3.produce_fc2_compact is False
    assert ph3.produce_fc3_compact is False
    with h5py.File(workflow.fc2_path, "r") as handle:
        assert handle["force_constants"].shape == (2, 2, 3, 3)
        assert "p2s_map" not in handle
    with h5py.File(workflow.fc3_path, "r") as handle:
        assert handle["fc3"].shape == (2, 2, 2, 3, 3, 3)
        assert "p2s_map" not in handle


def test_finite_displacement_fc2_only_skips_fc3(tmp_path):
    cfg = SimpleNamespace(
        result_dir=str(tmp_path),
        progress=False,
        compact_fc=True,
    )
    workflow = NEPPhononWorkflow(cfg)
    ph3 = FakePhono3pyForFC()
    workflow._make_phono3py = lambda: ph3

    workflow.run_finite_disp_fitting(include_fc3=False)

    assert ph3.generated_fc2_only is True
    assert ph3.generated_fc2fc3 is False
    assert ph3.produce_fc2_compact is True
    assert ph3.produce_fc3_compact is None
    assert workflow.fc2_path.exists()
    assert not workflow.fc3_path.exists()


def test_compute_kappa_generates_slurm_lbte_pipeline(tmp_path):
    cfg = SimpleNamespace(
        result_dir=str(tmp_path),
        mesh=[3, 3, 3],
        method="lbte",
        wigner=False,
        isotope=True,
        bfmp=100.0,
        temps=[100, 500, 100],
        progress=False,
        lbte_parallel={
            "backend": "slurm",
            "jobs": 3,
            "cpus_per_task": 4,
            "submit": False,
        },
    )
    workflow = NEPPhononWorkflow(cfg)
    workflow.fc2_path.write_text("fc2", encoding="utf-8")
    workflow.fc3_path.write_text("fc3", encoding="utf-8")
    workflow.disp_path.write_text("disp", encoding="utf-8")

    commands = []

    def fake_run(command, cwd=None):
        commands.append(command)
        if "--wgp" in command:
            (tmp_path / "ir_grid_points.yaml").write_text(
                """ir_grid_points:
- grid_point: 0
- grid_point: 3
- grid_point: 7
- grid_point: 9
- grid_point: 12
""",
                encoding="utf-8",
            )
        return 0

    workflow._run_command = fake_run
    workflow._phono3py_command = lambda: ["phono3py"]
    workflow._phono3py_needs_fc_flags = lambda: False

    manifest = workflow.compute_kappa()

    assert manifest["submitted"] is False
    assert manifest["grid_point_count"] == 5
    assert manifest["array_jobs"] == 3
    assert "--wgp" in commands[0]
    work_dir = tmp_path / "lbte-slurm"
    assert (work_dir / "grid-points" / "000.txt").read_text().strip() == "0,9"
    assert (work_dir / "grid-points" / "001.txt").read_text().strip() == "3,12"
    assert (work_dir / "grid-points" / "002.txt").read_text().strip() == "7"
    array_script = (work_dir / "array.sh").read_text(encoding="utf-8")
    collect_script = (work_dir / "collect.sh").read_text(encoding="utf-8")
    assert "#SBATCH --array=0-2" in array_script
    assert "--write-pp" in array_script
    assert '--gp "$GRID_POINTS"' in array_script
    assert "--read-pp" in collect_script
    assert "--isotope" in collect_script
    assert "--boundary-mfp 100.0" in collect_script
    assert "--tmin 100 --tmax 500 --tstep 100" in collect_script
    for script_name in ("prepare.sh", "array.sh", "collect.sh"):
        result = subprocess.run(
            ["bash", "-n", str(work_dir / script_name)],
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode == 0, result.stderr
