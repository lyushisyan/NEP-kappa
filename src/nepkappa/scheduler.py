"""Scheduler adapters used by NEP-kappa execution backends."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import subprocess


@dataclass(frozen=True)
class SchedulerSnapshot:
    """One non-mutating view of jobs known to a scheduler."""

    jobs: dict[str, str]
    source: str
    error: str | None = None


class SlurmScheduler:
    """Submit Slurm jobs and query their current queue states."""

    def __init__(self, *, run=None):
        self._run = run or subprocess.run

    def submit(self, path, cwd, dependency=None):
        """Submit one script and return the numeric Slurm job ID."""
        path = Path(path)
        command = ["sbatch", "--parsable"]
        if dependency:
            command.append(f"--dependency=afterok:{dependency}")
        command.append(str(path))
        try:
            result = self._run(
                command,
                cwd=cwd,
                text=True,
                capture_output=True,
                check=False,
            )
        except FileNotFoundError as exc:
            raise FileNotFoundError(
                "Slurm submission requested but 'sbatch' was not found. "
                "Set the relevant parallel.submit option to false to generate "
                "scripts only."
            ) from exc
        if result.returncode != 0:
            detail = result.stderr.strip() or result.stdout.strip()
            raise RuntimeError(f"sbatch failed for {path.name}: {detail}")
        job_id = result.stdout.strip().split(";", 1)[0]
        if not job_id:
            raise RuntimeError(f"sbatch returned no job id for {path.name}")
        return job_id

    def submit_chain(self, scripts, cwd):
        """Submit named scripts sequentially with ``afterok`` dependencies."""
        job_ids = {}
        dependency = None
        for name, script in scripts.items():
            job_id = self.submit(script, cwd, dependency=dependency)
            job_ids[name] = job_id
            dependency = job_id
        return job_ids

    def query(self, jobs):
        """Query ``squeue`` for a mapping of logical job names to IDs."""
        normalized = {str(name): str(job_id) for name, job_id in jobs.items()}
        if not normalized:
            return SchedulerSnapshot({}, source="squeue")
        command = [
            "squeue",
            "--noheader",
            "--format=%i|%T",
            "--jobs",
            ",".join(normalized.values()),
        ]
        try:
            result = self._run(
                command,
                text=True,
                capture_output=True,
                check=False,
            )
        except FileNotFoundError:
            return SchedulerSnapshot(
                {name: "unknown" for name in normalized},
                source="stored",
                error="squeue is not available",
            )
        if result.returncode != 0:
            detail = result.stderr.strip() or result.stdout.strip()
            return SchedulerSnapshot(
                {name: "unknown" for name in normalized},
                source="stored",
                error=f"squeue failed: {detail}",
            )

        by_id = {}
        for line in result.stdout.splitlines():
            job_id, separator, state = line.strip().partition("|")
            if separator and job_id:
                by_id[job_id] = state.lower()
        states = {
            name: by_id.get(job_id, "not-in-queue")
            for name, job_id in normalized.items()
        }
        return SchedulerSnapshot(states, source="squeue")
