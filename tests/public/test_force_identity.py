"""Persistent reuse must survive only unchanged scientific inputs and outputs."""
import json
from pathlib import Path
from types import SimpleNamespace

from ase import Atoms
import numpy as np
import pytest

from nepkappa.artifacts import ForceArtifactStore
from nepkappa.force_identity import command_identity, source_digest, provider_identity
from nepkappa.identity import data_sha256


@pytest.mark.parametrize('artifact', ['forces.npy', 'vasprun.xml'])
@pytest.mark.parametrize('change', ['finite-output', 'inputs', 'legacy', 'missing', 'opaque'])
def test_cache_rejects_unverified_records(tmp_path, artifact, change):
    payload = {'model': 'a', 'cache_reusable': change != 'opaque'}
    digest = data_sha256(payload)
    path = tmp_path / artifact
    if artifact == 'forces.npy':
        np.save(path, np.ones((1, 3)))
    else:
        path.write_text('<forces>1 2 3</forces>')
    ForceArtifactStore.record_inputs(tmp_path, digest, payload, artifact=artifact)
    assert ForceArtifactStore.matches_inputs(tmp_path, digest) == (change != 'opaque')
    record_path = tmp_path / 'job-input.json'
    record = json.loads(record_path.read_text())
    if change == 'finite-output':
        if artifact == 'forces.npy':
            np.save(path, np.full((1, 3), 7.))
        else:
            path.write_text('<forces>7 8 9</forces>')
    elif change == 'inputs':
        record['inputs']['model'] = 'b'
    elif change == 'legacy':
        record.pop('schema_version')
    elif change == 'missing':
        path.unlink()
    record_path.write_text(json.dumps(record))
    assert not ForceArtifactStore.matches_inputs(tmp_path, digest)


@pytest.mark.parametrize('change', ['positions', 'pbc', 'charges', 'magmoms'])
def test_exact_ase_state_invalidates_force_cache(tmp_path, change):
    store = ForceArtifactStore(vasp_root=tmp_path/'vasp', generic_root=tmp_path/'forces')
    atoms = Atoms('Si', cell=[5]*3, pbc=True)
    calls = []

    def calculate(atoms, *args):
        calls.append(True)
        return np.full((len(atoms), 3), len(calls), float)

    def fingerprint(*args):
        payload = {'model': 'a'}
        return data_sha256(payload), payload

    def evaluate():
        return store.calculate_or_reuse(atoms, 'fc2', 1, calculator='ase',
                                        fingerprint=fingerprint, calculate=calculate)

    np.testing.assert_array_equal(evaluate(), evaluate())
    assert len(calls) == 1
    if change == 'positions':
        atoms.positions[0, 0] += 1e-10  # Below the old structure-match tolerance.
    elif change == 'pbc':
        atoms.pbc[0] = False
    elif change == 'charges':
        atoms.set_initial_charges([.1])
    else:
        atoms.set_initial_magnetic_moments([1.])
    evaluate()
    assert len(calls) == 2


def test_command_identity_keeps_executable_and_file_arguments_distinct(tmp_path):
    engine, model = tmp_path/'engine', tmp_path/'model'
    engine.write_text('executable'); model.write_text('parameters')
    identity = command_identity([str(engine), str(model)])
    assert identity['executable']['path'] == str(engine)
    assert set(identity['files']) == {str(engine), str(model)}
    before = identity['files'][str(model)]['sha256']
    model.write_text('new parameters')
    assert command_identity([str(engine), str(model)])['files'][str(model)]['sha256'] != before


def test_force_source_identity_ignores_prose_but_tracks_force_changes(tmp_path):
    path = tmp_path/'workflow.py'
    original = '''class NEPPhononWorkflow:
    def force(self):
        """Initial prose."""
        return 1
    def plot(self):
        return 2
'''
    path.write_text(original)
    reference = source_digest(path, {'force'})
    path.write_text(original.replace('Initial prose.', 'Edited prose.').replace('return 2', 'return 3')+'\n# comment\n')
    assert source_digest(path, {'force'}) == reference
    path.write_text(original.replace('return 1', 'return 4'))
    assert source_digest(path, {'force'}) != reference


def test_same_version_provider_binary_change_invalidates_identity(tmp_path, monkeypatch):
    binary = tmp_path/'native.so'; binary.write_bytes(b'first build')
    distribution = SimpleNamespace(version='1.0', files=[Path('native.so')],
                                   locate_file=lambda entry: tmp_path/entry)
    monkeypatch.setattr('nepkappa.force_identity.metadata.distribution', lambda name: distribution)
    provider_identity.cache_clear()
    before = provider_identity('synthetic')
    binary.write_bytes(b'second build')
    provider_identity.cache_clear()  # A new process/run re-reads installed engines.
    after = provider_identity('synthetic')
    assert before['version'] == after['version']
    assert before['implementation'] != after['implementation']
    provider_identity.cache_clear()


def test_opaque_backend_reuses_only_within_one_workflow(tmp_path, monkeypatch):
    from nepkappa.workflow import NEPPhononWorkflow

    monkeypatch.setattr('nepkappa.force_identity.reusable_backend', lambda *args: False)
    config = SimpleNamespace(result_dir=str(tmp_path), calculator='ase',
                             calculator_factory='ase.calculators.lj:LennardJones',
                             calculator_kwargs={}, calculator_model_files=[])
    atoms = Atoms('Ar2', positions=[[0, 0, 0], [2, 0, 0]], cell=[10]*3, pbc=True)
    calls = []

    def calculate(*args):
        calls.append(True)
        return np.ones((2, 3))

    first = NEPPhononWorkflow(config)
    first._calculate_forces = calculate
    first._calculate_forces_cached(atoms, 'fc2', 1)
    first._calculate_forces_cached(atoms, 'fc2', 1)
    assert len(calls) == 1
    second = NEPPhononWorkflow(config)
    second._calculate_forces = calculate
    second._calculate_forces_cached(atoms, 'fc2', 1)
    assert len(calls) == 2
