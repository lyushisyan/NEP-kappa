"""Slurm orchestration for distributed phono3py LBTE calculations."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import shlex
import subprocess

import yaml


def run_lbte_slurm(
    config,
    output_dir,
    disp_path,
    phono3py_command,
    needs_fc_flags,
    run_command,
):
    """Generate and optionally submit a three-stage Slurm LBTE workflow."""
    output_dir = Path(output_dir).resolve()
    disp_path = Path(disp_path)
    settings = slurm_settings(config)
    mesh_args = ["--mesh", *(str(value) for value in config.mesh)]
    base = [*phono3py_command, disp_path.name]
    if needs_fc_flags:
        base.extend(["--fc2", "--fc3"])

    wgp_command = [*base, "--lbte", *mesh_args, "--wgp"]
    print("  - Discovering irreducible grid points")
    print(f"  - Running command: {shlex.join(wgp_command)}")
    ret = run_command(wgp_command, cwd=output_dir)
    if ret != 0:
        raise RuntimeError(f"Phono3py --wgp failed with return code {ret}")

    grid_file = output_dir / "ir_grid_points.yaml"
    grid_points = load_ir_grid_points(grid_file)
    chunks = split_grid_points(grid_points, settings["jobs"])
    work_dir = output_dir / "lbte-slurm"
    chunk_dir = work_dir / "grid-points"
    log_dir = work_dir / "logs"
    chunk_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    for index, chunk in enumerate(chunks):
        (chunk_dir / f"{index:03d}.txt").write_text(
            ",".join(str(value) for value in chunk) + "\n",
            encoding="utf-8",
        )

    scripts = write_slurm_scripts(
        config,
        settings,
        output_dir,
        work_dir,
        chunks,
        base,
        mesh_args,
    )
    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "backend": "slurm",
        "submitted": False,
        "grid_point_count": len(grid_points),
        "array_jobs": len(chunks),
        "scripts": {key: str(value) for key, value in scripts.items()},
    }
    manifest_path = work_dir / "submission.yaml"
    write_manifest(manifest_path, manifest)

    print(
        f"  - Split {len(grid_points)} irreducible grid points into "
        f"{len(chunks)} Slurm array tasks"
    )
    print(f"  - Slurm files written to {work_dir}")
    if settings["submit"]:
        try:
            job_ids = submit_slurm_pipeline(scripts, output_dir)
        except Exception as exc:
            manifest["submission_error"] = str(exc)
            write_manifest(manifest_path, manifest)
            raise
        manifest["submitted"] = True
        manifest["job_ids"] = job_ids
        print(f"  - Prepare job: {job_ids['prepare']}")
        print(f"  - Array job:   {job_ids['array']}")
        print(f"  - Collect job: {job_ids['collect']}")
        print("  - LBTE calculation submitted; the collect job runs after all array tasks.")
    else:
        print("  - Submission disabled; inspect the generated scripts before enabling submit.")

    write_manifest(manifest_path, manifest)
    return manifest


def write_manifest(path, manifest):
    """Write Slurm workflow state in a human-readable form."""
    path.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")


def slurm_settings(config):
    """Return normalized Slurm settings with conservative defaults."""
    raw = dict(getattr(config, "lbte_parallel", {}) or {})
    settings = {
        "jobs": int(raw.get("jobs", 32)),
        "max_concurrent": raw.get("max_concurrent"),
        "partition": raw.get("partition"),
        "account": raw.get("account"),
        "time": optional_string(raw.get("time")),
        "memory": optional_string(raw.get("memory")),
        "cpus_per_task": optional_int(raw.get("cpus_per_task")),
        "collect_time": optional_string(raw.get("collect_time", raw.get("time"))),
        "collect_memory": optional_string(raw.get("collect_memory", raw.get("memory"))),
        "collect_cpus_per_task": optional_int(
            raw.get("collect_cpus_per_task", raw.get("cpus_per_task"))
        ),
        "job_name": str(raw.get("job_name", "nepkappa-lbte")),
        "preamble": list(raw.get("preamble", [])),
        "extra_sbatch": list(raw.get("extra_sbatch", [])),
        "submit": bool(raw.get("submit", True)),
    }
    if settings["max_concurrent"] is not None:
        settings["max_concurrent"] = int(settings["max_concurrent"])
        if settings["max_concurrent"] < 1:
            raise ValueError("kappa.parallel.max_concurrent must be positive.")
    if settings["collect_cpus_per_task"] is not None and settings["collect_cpus_per_task"] < 1:
        raise ValueError("kappa.parallel.collect_cpus_per_task must be positive.")
    return settings


def optional_string(value):
    """Return a string override or None to preserve the Slurm default."""
    return None if value in (None, "") else str(value)


def optional_int(value):
    """Return an integer override or None to preserve the Slurm default."""
    return None if value is None else int(value)


def load_ir_grid_points(path):
    """Read irreducible grid-point indices written by phono3py --wgp."""
    if not path.exists():
        raise FileNotFoundError(f"Phono3py did not create {path}")
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    entries = data.get("ir_grid_points")
    if not isinstance(entries, list):
        raise ValueError(f"Invalid ir_grid_points data in {path}")
    points = []
    for entry in entries:
        if isinstance(entry, dict) and "grid_point" in entry:
            points.append(int(entry["grid_point"]))
        elif isinstance(entry, int):
            points.append(entry)
        else:
            raise ValueError(f"Invalid grid-point entry in {path}: {entry!r}")
    if not points:
        raise ValueError(f"No irreducible grid points found in {path}")
    return points


def split_grid_points(grid_points, jobs):
    """Split grid points round-robin, matching phono3py's recommended helper."""
    count = min(int(jobs), len(grid_points))
    chunks = [[] for _ in range(count)]
    for index, grid_point in enumerate(grid_points):
        chunks[index % count].append(grid_point)
    return chunks


def write_slurm_scripts(config, settings, output_dir, work_dir, chunks, base, mesh_args):
    """Write prepare, array-worker, and collect Slurm scripts."""
    common = [*base, *mesh_args]
    prepare_command = [*common, "--write-phonon"]
    worker_command = [
        *common,
        "--lbte",
        "--ts",
        str(config.temps[0]),
        "--write-pp",
        "--read-phonon",
    ]
    collect_command = [*common, "--lbte", "--read-pp", "--read-phonon"]
    collect_command.extend(transport_flags(config))
    collect_command.extend(temperature_flags(config.temps))

    array_range = f"0-{len(chunks) - 1}"
    if settings["max_concurrent"] is not None:
        array_range += f"%{settings['max_concurrent']}"

    scripts = {
        "prepare": work_dir / "prepare.sh",
        "array": work_dir / "array.sh",
        "collect": work_dir / "collect.sh",
    }
    scripts["prepare"].write_text(
        render_script(
            settings,
            f"{settings['job_name']}-prepare",
            "prepare-%j.out",
            output_dir,
            [shlex.join(prepare_command)],
        ),
        encoding="utf-8",
    )
    chunk_template = str((work_dir / "grid-points" / "%03d.txt").resolve())
    worker_lines = [
        f'GP_FILE=$(printf {shlex.quote(chunk_template)} "$SLURM_ARRAY_TASK_ID")',
        'GRID_POINTS=$(tr -d "\\n" < "$GP_FILE")',
        shlex.join(worker_command) + ' --gp "$GRID_POINTS"',
    ]
    scripts["array"].write_text(
        render_script(
            settings,
            f"{settings['job_name']}-worker",
            "worker-%A_%a.out",
            output_dir,
            worker_lines,
            array=array_range,
        ),
        encoding="utf-8",
    )
    scripts["collect"].write_text(
        render_script(
            settings,
            f"{settings['job_name']}-collect",
            "collect-%j.out",
            output_dir,
            [shlex.join(collect_command)],
            collect=True,
        ),
        encoding="utf-8",
    )
    for path in scripts.values():
        path.chmod(0o755)
    return scripts


def render_script(settings, job_name, output_name, output_dir, commands, array=None, collect=False):
    """Render one Slurm batch script."""
    cpus = settings["collect_cpus_per_task"] if collect else settings["cpus_per_task"]
    walltime = settings["collect_time"] if collect else settings["time"]
    memory = settings["collect_memory"] if collect else settings["memory"]
    error_name = output_name.removesuffix(".out") + ".err"
    lines = [
        "#!/usr/bin/env bash",
        f"#SBATCH --job-name={job_name}",
        f"#SBATCH --output=lbte-slurm/logs/{output_name}",
        f"#SBATCH --error=lbte-slurm/logs/{error_name}",
    ]
    if cpus is not None:
        lines.append(f"#SBATCH --cpus-per-task={cpus}")
    if walltime is not None:
        lines.append(f"#SBATCH --time={walltime}")
    if memory is not None:
        lines.append(f"#SBATCH --mem={memory}")
    if settings["partition"]:
        lines.append(f"#SBATCH --partition={settings['partition']}")
    if settings["account"]:
        lines.append(f"#SBATCH --account={settings['account']}")
    if array:
        lines.append(f"#SBATCH --array={array}")
    for option in settings["extra_sbatch"]:
        lines.append(f"#SBATCH {option}")
    lines.extend(
        [
            "",
            "set -euo pipefail",
            f"export OMP_NUM_THREADS=${{SLURM_CPUS_PER_TASK:-{cpus or 1}}}",
            *settings["preamble"],
            f"cd {shlex.quote(str(output_dir))}",
            *commands,
            "",
        ]
    )
    return "\n".join(lines)


def transport_flags(config):
    """Return final LBTE transport flags shared with the serial workflow."""
    flags = []
    if config.wigner:
        flags.extend(["--tt", "wte"])
    if config.isotope:
        flags.append("--isotope")
    if config.bfmp is not None:
        flags.extend(["--boundary-mfp", str(config.bfmp)])
    return flags


def temperature_flags(temperatures):
    """Return phono3py temperature command-line flags."""
    if len(temperatures) == 3:
        tmin, tmax, tstep = temperatures
        return ["--tmin", str(tmin), "--tmax", str(tmax), "--tstep", str(tstep)]
    return ["--ts", str(temperatures[0])]


def submit_slurm_pipeline(scripts, cwd):
    """Submit Slurm scripts with afterok dependencies."""
    prepare = submit_script(scripts["prepare"], cwd)
    array = submit_script(scripts["array"], cwd, dependency=prepare)
    collect = submit_script(scripts["collect"], cwd, dependency=array)
    return {"prepare": prepare, "array": array, "collect": collect}


def submit_script(path, cwd, dependency=None):
    """Submit one script with sbatch and return its numeric job id."""
    command = ["sbatch", "--parsable"]
    if dependency:
        command.append(f"--dependency=afterok:{dependency}")
    command.append(str(path))
    try:
        result = subprocess.run(
            command,
            cwd=cwd,
            text=True,
            capture_output=True,
            check=False,
        )
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            "Slurm submission requested but 'sbatch' was not found. "
            "Set kappa.parallel.submit: false to generate scripts only."
        ) from exc
    if result.returncode != 0:
        detail = result.stderr.strip() or result.stdout.strip()
        raise RuntimeError(f"sbatch failed for {path.name}: {detail}")
    job_id = result.stdout.strip().split(";", 1)[0]
    if not job_id:
        raise RuntimeError(f"sbatch returned no job id for {path.name}")
    return job_id
