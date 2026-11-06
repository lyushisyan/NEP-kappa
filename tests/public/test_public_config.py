"""Public YAML contract checks: no calculation or external program is launched."""
import argparse
from pathlib import Path

import pytest
import yaml

from nepkappa import config as inputs
from nepkappa.config import parse_selected_command, parse_workflow_args

ROOT = Path(__file__).resolve().parents[2]
PLANS = {
    'nep-rta-wigner-3ph.yaml': ['relax', 'fc2fc3', 'kappa'],
    'nep-lbte-wigner-3ph.yaml': ['relax', 'fc2fc3', 'kappa'],
    'vasp-rta-3ph.yaml': ['relax', 'fc2fc3', 'kappa'],
    'nep-rta-3ph-4ph.yaml': ['relax', 'fc2fc3', 'fc4', 'kappa4'],
    'nep-rta-wigner-3ph-4ph.yaml': ['relax', 'fc2fc3', 'fc4', 'kappa4'],
    'nep-lbte-3ph-rta-4ph.yaml': ['relax', 'fc2fc3', 'fc4', 'kappa4'],
    'nep-lbte-3ph-lbte-4ph.yaml': ['relax', 'fc2fc3', 'fc4', 'kappa4'],
    'nep-qha-rta-3ph.yaml': ['qha', 'qha-kappa'],
    'nep-qha-sscha-rta-3ph.yaml': ['qha', 'qha-sscha'],
    'tdbte.yaml': ['fc2fc3', 'tdbte'],
}


@pytest.mark.parametrize('name,steps', PLANS.items())
def test_public_examples_select_scientific_stages(name, steps):
    cfg = parse_workflow_args(ROOT/'examples'/name, command='run')
    assert cfg.workflow_steps == steps
    assert tuple(cfg.workflow_steps) == cfg.sections.workflow.steps
    assert tuple(cfg.dim_fc2) == cfg.sections.force_constants.dim_fc2
    assert cfg.config_path == str(ROOT/'examples'/name)


def write_input(tmp_path, data):
    path = tmp_path/'input.yaml'
    path.write_text(yaml.safe_dump(data))
    return path


def test_compact_commands_preserve_defaults_and_fresh_section_views(tmp_path):
    path = write_input(tmp_path, {'structure': {}, 'calculator': {'model': 'model.nep'}, 'force-constant': {}})
    cfg = parse_workflow_args(path, command='fc2')
    assert cfg.dim_fc2 == cfg.dim_fc3 == cfg.dim_fc4 == [4, 4, 1]
    assert cfg.mesh == [21, 21, 1]
    assert cfg.calculator == 'nep'
    cfg.result_dir = 'changed'
    assert cfg.sections.output.result_dir == 'changed'
    assert 'sections' not in vars(cfg)
    with pytest.raises(ValueError, match='six sections'):
        parse_workflow_args(path, command='run')


@pytest.mark.parametrize('name,command,selected', [
    ('nep-rta-wigner-3ph.yaml', 'fc2', 'fc2'),
    ('nep-qha-rta-3ph.yaml', 'fc2', 'qha'),
    ('nep-rta-wigner-3ph.yaml', 'kappa', 'kappa'),
    ('nep-rta-3ph-4ph.yaml', 'kappa', 'kappa4'),
    ('nep-qha-rta-3ph.yaml', 'kappa', 'qha-kappa'),
    ('nep-qha-sscha-rta-3ph.yaml', 'kappa', 'qha-sscha'),
])
def test_input_selected_command_reads_and_parses_once(monkeypatch, name, command, selected):
    counts = {'yaml': 0, 'argparse': 0}
    load, parse = inputs.load_yaml_strict, argparse.ArgumentParser.parse_args
    def counted_load(*args, **kwargs):
        counts['yaml'] += 1
        return load(*args, **kwargs)
    def counted_parse(*args, **kwargs):
        counts['argparse'] += 1
        return parse(*args, **kwargs)
    monkeypatch.setattr(inputs, 'load_yaml_strict', counted_load)
    monkeypatch.setattr(argparse.ArgumentParser, 'parse_args', counted_parse)
    actual, cfg = parse_selected_command(ROOT/'examples'/name, command)
    assert actual == selected
    assert counts == {'yaml': 1, 'argparse': 1}
    assert cfg.mesh is not None  # Plot preflight must not erase calculation defaults.


@pytest.mark.parametrize('text', [
    'plot: {}\nplot: {}\n',
    'plot: {tdbte-times: [0, 1, 2], tdbte_times: [0, 1, 2]}\n',
    'plot: {misspelled-option: true}\n',
])
def test_invalid_or_duplicate_yaml_is_rejected(tmp_path, text):
    path = tmp_path/'bad.yaml'; path.write_text(text)
    with pytest.raises(ValueError):
        parse_workflow_args(path, command='plot')


def test_plot_and_kappa_apply_different_requirements(tmp_path):
    path = write_input(tmp_path, {'plot': {}, 'calculator': {'name': 'mace'}})
    assert parse_workflow_args(path, command='plot').mesh is None
    with pytest.raises(ValueError, match='missing: kappa'):
        parse_workflow_args(path, command='kappa')
    path = write_input(tmp_path, {'kappa': {'temps': [100, 300, 70]}})
    with pytest.raises(SystemExit):
        parse_selected_command(path, 'kappa')


def test_yaml_aliases_and_nested_values(tmp_path):
    cfg = parse_workflow_args(write_input(tmp_path, {
        'structure': {'relaxation': False},
        'calculator': {'name': 'ase', 'factory': 'ase.calculators.lj:LennardJones',
                       'kwargs': {'epsilon': .01, 'sigma': 3.4}},
        'force_constant': {'dim_fc2': [2, 1, 1]},
    }), command='fc2')
    assert cfg.dim_fc2 == [2, 1, 1]
    assert cfg.calculator_kwargs == {'epsilon': .01, 'sigma': 3.4}
