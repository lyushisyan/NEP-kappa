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


def main(argv=None):
    """Update a state manifest from a successful scheduler collection job."""
    parser = argparse.ArgumentParser(description="Update NEP-kappa run state")
    parser.add_argument("--complete", required=True, help="submission.yaml path")
    args = parser.parse_args(argv)
    RunStateStore(args.complete).mark_complete()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
