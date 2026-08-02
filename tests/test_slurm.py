from pathlib import Path
from types import SimpleNamespace
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import nepkappa.slurm as slurm


def test_minimal_slurm_config_preserves_scheduler_resource_defaults():
    settings = slurm.slurm_settings(
        SimpleNamespace(lbte_parallel={"backend": "slurm"})
    )

    assert settings["jobs"] == 32
    assert settings["time"] is None
    assert settings["memory"] is None
    assert settings["cpus_per_task"] is None
    assert settings["submit"] is True


def test_submit_slurm_pipeline_chains_afterok_dependencies(tmp_path, monkeypatch):
    scripts = {
        "prepare": tmp_path / "prepare.sh",
        "array": tmp_path / "array.sh",
        "collect": tmp_path / "collect.sh",
    }
    calls = []

    def fake_submit(path, cwd, dependency=None):
        calls.append((path.name, dependency))
        return {"prepare.sh": "101", "array.sh": "102", "collect.sh": "103"}[path.name]

    monkeypatch.setattr(slurm, "submit_script", fake_submit)

    job_ids = slurm.submit_slurm_pipeline(scripts, tmp_path)

    assert calls == [
        ("prepare.sh", None),
        ("array.sh", "101"),
        ("collect.sh", "102"),
    ]
    assert job_ids == {"prepare": "101", "array": "102", "collect": "103"}
