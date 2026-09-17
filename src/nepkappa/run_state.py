"""Persistent and queryable state for deferred NEP-kappa calculations."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

import yaml


RUN_STATE_SCHEMA_VERSION = 1


def utc_now():
    return datetime.now(timezone.utc).isoformat()


class RunStateStore:
    """Read and atomically update one human-readable run-state manifest."""

    def __init__(self, path):
        self.path = Path(path)

    def read(self):
        if not self.path.is_file():
            return {}
        data = yaml.safe_load(self.path.read_text(encoding="utf-8")) or {}
        if not isinstance(data, dict):
            raise ValueError(f"Run-state manifest must be a mapping: {self.path}")
        return data

    def write(self, state):
        document = dict(state)
        document.setdefault("schema_version", RUN_STATE_SCHEMA_VERSION)
        document["updated_at"] = utc_now()
        self.path.parent.mkdir(parents=True, exist_ok=True)
        temporary = self.path.with_suffix(self.path.suffix + ".tmp")
        temporary.write_text(
            yaml.safe_dump(document, sort_keys=False), encoding="utf-8"
        )
        temporary.replace(self.path)
        return document

    def update(self, **changes):
        state = self.read()
        state.update(changes)
        return self.write(state)

    def mark_submitted(self, job_ids):
        return self.update(
            submitted=True,
            status="submitted",
            submitted_at=utc_now(),
            job_ids={str(name): str(value) for name, value in job_ids.items()},
        )

    def mark_submission_failed(self, error):
        return self.update(status="submission-failed", submission_error=str(error))

    def mark_complete(self):
        return self.update(status="complete", completed_at=utc_now())


def discover_run_states(output_dir, *, scheduler=None):
    """Return all scheduler manifests below one result directory."""
    output_dir = Path(output_dir).resolve()
    records = []
    for path in sorted(output_dir.rglob("submission.yaml")):
        store = RunStateStore(path)
        state = store.read()
        record = {
            "path": str(path),
            "workflow": state.get("workflow", _infer_workflow(path)),
            "backend": state.get("backend", "unknown"),
            "status": _stored_status(state),
            "submitted": bool(state.get("submitted", False)),
            "job_ids": _job_ids(state),
        }
        if scheduler is not None and record["job_ids"]:
            snapshot = scheduler.query(record["job_ids"])
            record["scheduler"] = snapshot.source
            record["jobs"] = snapshot.jobs
            if snapshot.error:
                record["query_error"] = snapshot.error
        records.append(record)
    return records


def format_run_states(records):
    """Format run-state records for the terminal status command."""
    if not records:
        return "No scheduler submissions were found."
    lines = []
    for record in records:
        lines.append(
            f"{record['workflow']}: {record['status']} "
            f"({Path(record['path']).parent})"
        )
        for name, job_id in record["job_ids"].items():
            live = record.get("jobs", {}).get(name)
            suffix = f" [{live}]" if live else ""
            lines.append(f"  - {name}: {job_id}{suffix}")
        if record.get("query_error"):
            lines.append(f"  - scheduler query: {record['query_error']}")
    return "\n".join(lines)


def _job_ids(state):
    job_ids = state.get("job_ids")
    if isinstance(job_ids, dict):
        return {str(name): str(value) for name, value in job_ids.items()}
    if state.get("job_id") is not None:
        return {"main": str(state["job_id"])}
    return {}


def _stored_status(state):
    if state.get("status"):
        return str(state["status"])
    if state.get("submission_error"):
        return "submission-failed"
    if state.get("submitted"):
        return "submitted"
    return "prepared"


def _infer_workflow(path):
    parent = path.parent.name
    return {
        "force-slurm": "force-constants",
        "lbte-slurm": "lbte",
        "fourphonon": "fourphonon",
    }.get(parent, parent)


def main(argv=None):
    """Update a state manifest from a successful scheduler collection job."""
    parser = argparse.ArgumentParser(description="Update NEP-kappa run state")
    parser.add_argument("--complete", required=True, help="submission.yaml path")
    args = parser.parse_args(argv)
    RunStateStore(args.complete).mark_complete()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
