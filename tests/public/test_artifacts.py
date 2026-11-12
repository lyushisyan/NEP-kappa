"""Unit tests for persistent calculation artifacts."""

import hashlib

from nepkappa.identity import data_sha256

from ase import Atoms
import numpy as np
import pytest

from nepkappa.artifacts import ForceArtifactStore


def make_store(tmp_path):
    return ForceArtifactStore(
        vasp_root=tmp_path / "vasp-runs",
        generic_root=tmp_path / "force-jobs",
    )


def test_force_job_directory_depends_on_calculator(tmp_path):
    store = make_store(tmp_path)

    assert store.job_dir("FC2", 3, calculator="vasp") == (
        tmp_path / "vasp-runs" / "fc2" / "00003"
    )
    assert store.job_dir("FC2", 3, calculator="nep") == (
        tmp_path / "force-jobs" / "fc2" / "00003"
    )


def test_input_markers_are_auditable_and_queryable(tmp_path):
    store = make_store(tmp_path)
    job_dir = tmp_path / "job"
    payload = {"calculator": "nep", "model": "example.nep"}

    job_dir.mkdir()
    np.save(job_dir / "forces.npy", np.zeros((1, 3)))
    digest = data_sha256(payload)
    store.record_inputs(job_dir, digest, payload)

    assert store.matches_inputs(job_dir, digest)
    assert not store.matches_inputs(job_dir, "digest-b")
    assert '"calculator": "nep"' in (job_dir / "job-input.json").read_text(
        encoding="utf-8"
    )


def test_calculate_or_reuse_keeps_a_matching_force_result(tmp_path):
    store = make_store(tmp_path)
    atoms = Atoms("Si", positions=[[0.0, 0.0, 0.0]], cell=[5.0] * 3, pbc=True)
    calls = []

    def fingerprint(path, backend):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        return digest, {"backend": backend, "digest": digest}

    def calculate(structure, label, index):
        calls.append((label, index))
        return np.full((len(structure), 3), 2.0)

    first = store.calculate_or_reuse(
        atoms,
        "fc3",
        1,
        calculator="nep",
        fingerprint=fingerprint,
        calculate=calculate,
    )
    second = store.calculate_or_reuse(
        atoms,
        "fc3",
        1,
        calculator="nep",
        fingerprint=fingerprint,
        calculate=calculate,
    )

    assert calls == [("fc3", 1)]
    np.testing.assert_array_equal(first, second)


def test_invalid_force_shape_is_rejected(tmp_path):
    store = make_store(tmp_path)

    with pytest.raises(RuntimeError, match=r"expected \(2, 3\)"):
        store.validate_forces(np.zeros((1, 3)), 2, "test result")
