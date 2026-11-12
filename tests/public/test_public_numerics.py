"""Small conservation, shape and cache checks using synthetic local data only."""
from types import SimpleNamespace

import numpy as np
import pytest

from nepkappa.artifacts import ForceArtifactStore
from nepkappa.provenance import data_sha256, file_sha256
from nepkappa.tdbte import ShellKernel


@pytest.mark.parametrize('vertices', [2, 4])
def test_collision_conserves_energy_and_bose_equilibrium(vertices):
    k = ShellKernel(np.array([3., 1., 2.]), np.array([0]),
        np.full((1, vertices), 1), np.full((1, vertices), 2),
        np.full((1, vertices), 1/vertices), np.array([.2]), np.zeros(3, bool))
    k.validate()
    equilibrium = k.bose(300.)
    np.testing.assert_allclose(k.collision(equilibrium)[0], 0, atol=1e-12)
    excited = equilibrium * [1.1, 1., 1.]
    derivative, _ = k.collision(excited)
    assert abs(k.frequency @ derivative) < 1e-12
    assert np.log1p(1/excited) @ derivative >= -1e-12
    with pytest.raises(ValueError, match='positive'):
        k.collision(np.array([-1., 1., 2.]))


@pytest.mark.parametrize('forces', [np.zeros((2, 2)), [[np.nan, 0, 0]], [[np.inf, 0, 0]]])
def test_invalid_force_results_are_rejected(forces):
    with pytest.raises(RuntimeError):
        ForceArtifactStore.validate_forces(forces, 1, 'synthetic')


def test_force_array_preserves_components_and_atom_order():
    forces = np.array([[1., -2., 3.], [-4., 5., -6.]])
    np.testing.assert_array_equal(ForceArtifactStore.validate_forces(forces, 2, 'synthetic'), forces)


def test_force_fingerprint_tracks_model_structure_and_source(tmp_path):
    from nepkappa import workflow
    model, structure = tmp_path/'model.nep', tmp_path/'POSCAR'
    model.write_text('model-a'); structure.write_text('displacement-a')
    cfg = SimpleNamespace(result_dir=str(tmp_path/'out'), calculator='nep', nep_model=str(model))
    core = workflow.NEPPhononWorkflow(cfg)
    first, payload = core._job_input_fingerprint(structure, 'fc2')
    assert first == data_sha256(payload)
    assert payload['schema_version'] == 3
    assert payload['code']['force_implementation']['workflow.force_methods']
    assert core._job_input_fingerprint(structure, 'fc2')[0] == first
    model.write_text('model-b')
    second = core._job_input_fingerprint(structure, 'fc2')[0]
    assert second != first
    structure.write_text('displacement-b')
    assert core._job_input_fingerprint(structure, 'fc2')[0] != second


def test_provenance_keeps_inputs_timings_and_artifact_hashes(tmp_path):
    import yaml
    from nepkappa.provenance import ProvenanceRecorder
    source = tmp_path/'input.yaml'; source.write_text('output: {}\n')
    model = tmp_path/'model.nep'; model.write_text('synthetic model')
    output = tmp_path/'out'
    cfg = SimpleNamespace(result_dir=str(output), calculator='nep', nep_model=str(model), mesh=[3, 3, 3])
    recorder = ProvenanceRecorder(output, 'fc2fc3', source, cfg)
    manifest = recorder.start()
    artifact = output/'fc2.hdf5'; artifact.write_bytes(b'synthetic artifact')
    recorder.finish('complete', stage_timings=[('Finite displacement', 1.25)])
    record = yaml.safe_load(manifest.read_text())
    assert record['software']['scipy']
    assert record['status'] == 'complete'
    assert record['configuration']['mesh'] == [3, 3, 3]
    assert record['stage_timings_seconds']['Finite displacement'] == 1.25
    assert record['inputs']['nep_model']['sha256'] == file_sha256(model)
    assert record['artifacts']['fc2.hdf5']['sha256'] == file_sha256(artifact)
