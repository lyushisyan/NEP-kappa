"""Excitation type and strength have the same explicit units in both TD-BTE routes."""
import json

import numpy as np
import pytest
import yaml

from nepkappa.config import parse_workflow_args
from nepkappa.tdbte import ShellKernel, run_tdbte


def make_input(tmp_path, excitation, delta=None):
    k = ShellKernel(np.array([2., 1.]), np.array([0]), np.array([[1]*4]),
                    np.array([[1]*4]), np.full((1, 4), .25), np.array([.2]), np.zeros(2, bool))
    path = tmp_path/'kernel.npz'
    np.savez(path, frequency_THz=k.frequency, zero_modes=k.zero, parent=k.parent,
             daughter1=k.daughter1, daughter2=k.daughter2, weights=k.weights, coefficients=k.coefficients)
    path.with_suffix('.json').write_text(json.dumps(dict(
        model='synthetic excitation test', coefficient_time_unit='ps',
        mode_branch_indices=[0, 1], uniform_mode_weight=1., mesh=[2, 2, 2])))
    tdbte = dict(kernel=str(path), temperature=300, branches=[0], excitation=excitation,
                 **{'duration-ps': .001, 'max-step-ps': .0001, 'samples': 3})
    if delta is not None:
        tdbte['delta'] = delta
    source = tmp_path/'input.yaml'
    source.write_text(yaml.safe_dump({'tdbte': tdbte, 'output': {'result-dir': str(tmp_path/'result')}}))
    return k, source


@pytest.mark.parametrize('excitation,delta', [('relative', .01), ('temperature', 15.)])
def test_homogeneous_excitation_initial_occupations_and_audit(tmp_path, excitation, delta):
    k, source = make_input(tmp_path, excitation, delta)
    cfg = parse_workflow_args(source, command='tdbte')
    assert cfg.sections.tdbte.excitation == excitation
    assert cfg.sections.tdbte.delta == delta
    report = run_tdbte(cfg)
    with np.load(tmp_path/'result/tdbte/trajectories.npz') as data:
        expected = k.bose(300.)
        expected[0] = expected[0]*(1+delta) if excitation == 'relative' else k.bose(300+delta)[0]
        np.testing.assert_allclose(data['excited'][:, 0], expected, rtol=1e-14)
    assert report['numerical_checks_passed']
    assert report['excitation'] == excitation and report['delta'] == delta


@pytest.mark.parametrize('excitation,expected_delta', [('relative', .01), ('temperature', 1.), (.02, .02)])
def test_explicit_defaults_and_numeric_input_normalization(tmp_path, excitation, expected_delta):
    _, source = make_input(tmp_path, excitation)
    cfg = parse_workflow_args(source, command='tdbte')
    assert cfg.tdbte_excitation == ('relative' if isinstance(excitation, float) else excitation)
    assert cfg.tdbte_delta == expected_delta


@pytest.mark.parametrize('excitation,delta', [('relative', -1), ('temperature', float('nan')),
                                           ('invalid', .01), (.01, .02)])
def test_invalid_excitation_parameters_fail_before_running(tmp_path, excitation, delta):
    _, source = make_input(tmp_path, excitation, delta)
    with pytest.raises((SystemExit, ValueError)):
        parse_workflow_args(source, command='tdbte')


def test_spatial_relative_default_is_independent_of_homogeneous_delta(tmp_path):
    from nepkappa.tdbte_spatial import spatial_options
    _, source = make_input(tmp_path, 'relative', .2)
    data = yaml.safe_load(source.read_text())
    data['tdbte']['spatial'] = {'length-nm': 100, 'cells': 20,
        'initial': {'profile': 'ttg', 'modes': 'all', 'excitation': 'relative'}}
    source.write_text(yaml.safe_dump(data))
    cfg = parse_workflow_args(source, command='tdbte')
    assert spatial_options(cfg)['initial']['delta'] == .01
