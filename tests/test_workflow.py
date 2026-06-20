from pathlib import Path
from types import SimpleNamespace
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nepkappa.workflow import NEPPhononWorkflow


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
