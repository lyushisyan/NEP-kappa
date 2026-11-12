"""Slurm worker for cached displaced-structure force calculations."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from ase.io import read
from ase.io.jsonio import decode

from nepkappa.config import parse_workflow_args
from nepkappa.workflow import NEPPhononWorkflow
from nepkappa.identity import file_sha256
from nepkappa.run_state import RunStateStore


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="Absolute workflow YAML path")
    parser.add_argument("--task-file", required=True, help="JSON list for this array task")
    return parser


def run_task_file(config_path, task_path):
    """Evaluate all force jobs listed in one array-task JSON file."""
    config_path = Path(config_path).resolve()
    task_path = Path(task_path).resolve()
    cfg = parse_workflow_args(config_path, command="fc2fc3")
    workflow = NEPPhononWorkflow(cfg)
    tasks = json.loads(task_path.read_text(encoding="utf-8"))
    if not isinstance(tasks, list) or not tasks:
        raise ValueError(f"Force task file is empty or invalid: {task_path}")
    scopes = {task.get('reuse_scope') for task in tasks if isinstance(task, dict)}
    if len(scopes) != 1:
        raise ValueError("Force tasks must belong to one calculation")
    scope = scopes.pop()
    if scope is not None:
        state = RunStateStore(workflow.output_dir / 'force-slurm' / 'submission.yaml').read()
        if state and scope != state.get('reuse_scope'):
            raise ValueError("Force tasks belong to a different Slurm submission")
        workflow._force_reuse_scope = scope

    jobs = []
    for task in tasks:
        if not isinstance(task, dict):
            raise ValueError(f"Invalid force task entry in {task_path}: {task!r}")
        label = str(task.get("label", "")).lower()
        if label not in {"fc2", "fc3", "hiphive"}:
            raise ValueError(f"Unsupported force task label: {label!r}")
        index = int(task["index"])
        expected = workflow._force_job_dir(label, index) / "POSCAR"
        poscar = Path(task["poscar"]).resolve()
        if poscar != expected.resolve():
            raise ValueError(
                f"Force task POSCAR does not match its cache location: {poscar}"
            )
        atoms = read(str(poscar), format="vasp")
        if "atoms_sha256" in task:
            state = poscar.with_name("atoms.json")
            if file_sha256(state) != task["atoms_sha256"]:
                raise ValueError("Staged ASE state checksum mismatch")
            exact = decode(state.read_text(encoding="utf-8"))
            if not workflow._make_force_store().structures_match(atoms, exact):
                raise ValueError("Staged ASE state differs from its POSCAR")
            atoms = exact
        if cfg.nep_force_backend == "gpu-batch":
            jobs.append((atoms, label, index))
        else:
            workflow._calculate_forces_cached(atoms, label, index)
    if jobs:
        workflow._calculate_nep_batch_cached(jobs)
    print(f"Completed {len(tasks)} force job(s) from {task_path}")


def main(argv=None):
    args = build_parser().parse_args(argv)
    run_task_file(args.config, args.task_file)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
