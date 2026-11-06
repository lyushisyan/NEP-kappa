"""Scientific stage ordering and deferred Slurm control without a scheduler."""
from types import SimpleNamespace

import pytest

from nepkappa.application import CommandExecution, WorkflowStageRunner, execute_workflow_plan


@pytest.mark.parametrize('backend,hiphive,expected', [
    ('phono3py', False, [('finite', True)]),
    ('thirdorder', False, [('finite', False), ('thirdorder', True)]),
    ('phono3py', True, [('hiphive', True)]),
])
def test_force_constant_ordering(monkeypatch, backend, hiphive, expected):
    from nepkappa.stages import force_constants
    events = []
    core = SimpleNamespace(load_force_constant_structure=lambda: events.append('structure'))
    for name, label in [('FiniteDisplacementStage', 'finite'), ('HiPhiveStage', 'hiphive'), ('ThirdOrderStage', 'thirdorder')]:
        def factory(host, include_fc3=True, label=label):
            assert host is core
            return SimpleNamespace(run=lambda: events.append((label, include_fc3)))
        monkeypatch.setattr(force_constants, name, factory)
    cfg = SimpleNamespace(force_parallel={}, fc3_backend=backend, use_hiphive=hiphive)
    WorkflowStageRunner(cfg, CommandExecution(), core=core).force_constants(include_fc3=True)
    assert events == ['structure', *expected]


@pytest.mark.parametrize('environment,expected', [
    ({}, ['relax', 'fc2fc3']),
    ({'NEPKAPPA_FORCE_COLLECT': '1'}, ['fc2fc3', 'kappa']),
    ({'NEPKAPPA_RESUME_AFTER': 'fc2fc3'}, ['kappa']),
])
def test_plan_defers_and_resumes_without_repeating_completed_stages(monkeypatch, environment, expected):
    events = []
    for key in ['NEPKAPPA_FORCE_COLLECT', 'NEPKAPPA_RESUME_AFTER']:
        monkeypatch.delenv(key, raising=False)
    for key, value in environment.items():
        monkeypatch.setenv(key, value)
    def run(self, step):
        events.append(step)
        if step == 'fc2fc3' and not environment:
            return self.execution.run_stage('Slurm force preparation', lambda: {'submitted': True, 'status': 'submitted'})
    monkeypatch.setattr(WorkflowStageRunner, 'run', run)
    execute_workflow_plan(SimpleNamespace(workflow_steps=['relax', 'fc2fc3', 'kappa']))
    assert events == expected


def test_failed_stage_propagates_and_records_timing():
    execution = CommandExecution()
    def fail():
        raise RuntimeError('solver failed')
    with pytest.raises(RuntimeError, match='solver failed'):
        execution.run_stage('FC2', fail)
    assert execution.stage_timings[0][0] == 'FC2'
    assert execution.stage_timings[0][1] >= 0


def test_force_worker_still_uses_shared_context(tmp_path, monkeypatch):
    import json
    from nepkappa import force_worker
    from ase import Atoms
    events = []
    directory = tmp_path/'force-jobs/fc2/00001'; directory.mkdir(parents=True)
    poscar = directory/'POSCAR'
    Atoms('C', cell=[3, 3, 3], pbc=True).write(poscar, format='vasp')
    tasks = tmp_path/'tasks.json'
    tasks.write_text(json.dumps([{'label': 'fc2', 'index': 1, 'poscar': str(poscar)}]))
    monkeypatch.setattr(force_worker, 'parse_workflow_args', lambda *a, **k: SimpleNamespace(nep_force_backend='cpu'))
    context = SimpleNamespace(_force_job_dir=lambda *a: directory,
        _calculate_forces_cached=lambda atoms, label, index: events.append((len(atoms), label, index)))
    monkeypatch.setattr(force_worker, 'NEPPhononWorkflow', lambda cfg: context)
    force_worker.run_task_file(tmp_path/'input.yaml', tasks)
    assert events == [(1, 'fc2', 1)]
