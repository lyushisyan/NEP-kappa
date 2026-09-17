"""FourPhonon (ShengBTE extension) transport workflow."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

from ase.io import read
import numpy as np
import yaml

from nepkappa.provenance import file_identity
from nepkappa.run_state import RunStateStore
from nepkappa.slurm import submit_script


KAPPA_FILES = {
    "rta": "BTE.KappaTensorVsT_RTA",
    "iterative": "BTE.KappaTensorVsT_CONV",
}


class FourPhononWorkflow:
    """Stage ShengBTE-format IFCs and run 3ph+4ph transport."""

    def __init__(self, config):
        self.cfg = config
        sections = config.sections
        self.settings = sections.fourphonon
        self.structure = sections.structure
        self.output_dir = Path(sections.output.result_dir).resolve()
        self.workdir = self.output_dir / self.settings.workdir
        self.manifest_path = self.workdir / "submission.yaml"

    def run(self):
        self.workdir.mkdir(parents=True, exist_ok=True)
        inputs = self._stage_inputs()
        control = self._write_control()
        command = self._run_command()
        parallel = self.settings.parallel
        scheduler_backend = str(parallel.get("backend", "none")).lower()
        manifest = {
            "created_at": datetime.now(timezone.utc).isoformat(),
            "workflow": "fourphonon",
            "backend": scheduler_backend,
            "status": "prepared",
            "solver": self.settings.solver,
            "submitted": False,
            "command": command,
            "control": file_identity(control),
            "inputs": {name: file_identity(path) for name, path in inputs.items()},
        }
        self._write_manifest(manifest)

        if scheduler_backend == "slurm":
            script = self._write_slurm_script(command)
            manifest["script"] = str(script)
            if bool(parallel.get("submit", True)):
                try:
                    job_id = submit_script(script, self.workdir)
                except Exception as exc:
                    manifest["status"] = "submission-failed"
                    manifest["submission_error"] = str(exc)
                    self._write_manifest(manifest)
                    raise
                manifest["submitted"] = True
                manifest["status"] = "submitted"
                manifest["job_id"] = job_id
                manifest["job_ids"] = {"main": job_id}
                print(f"  - FourPhonon Slurm job: {job_id}")
            else:
                print(f"  - Submission disabled; inspect {script}")
            self._write_manifest(manifest)
            return manifest

        log_path = self.workdir / "fourphonon.log"
        print(f"  - Running: {shlex.join(command)}")
        return_code = run_streaming(command, self.workdir, log_path)
        if return_code != 0:
            raise RuntimeError(
                f"FourPhonon failed with return code {return_code}; see {log_path}."
            )
        summary = collect_outputs(
            self.workdir,
            preferred=self._preferred_solution(),
            allow_no_kappa=self.settings.only_harmonic,
        )
        manifest["status"] = "complete"
        manifest["summary"] = str(self.workdir / "fourphonon-summary.yaml")
        self._write_manifest(manifest)
        return summary

    def _input_path(self, configured, default_name):
        path = Path(configured).expanduser() if configured else self.output_dir / default_name
        if not path.is_absolute():
            path = Path.cwd() / path
        path = path.resolve()
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(
                f"Missing FourPhonon input {path}. Generate ShengBTE-format "
                "FC2/FC3/FC4 first or set its path in the fourphonon section."
            )
        return path

    def _stage_inputs(self):
        harmonic_name = (
            "espresso.ifc2"
            if self.settings.harmonic_format == "espresso"
            else "FORCE_CONSTANTS_2ND"
        )
        sources = {
            harmonic_name: self._input_path(self.settings.fc2, "FORCE_CONSTANTS_2ND"),
            "FORCE_CONSTANTS_3RD": self._input_path(self.settings.fc3, "FORCE_CONSTANTS_3RD"),
            "FORCE_CONSTANTS_4TH": self._input_path(self.settings.fc4, "FORCE_CONSTANTS_4TH"),
        }
        for name, source in sources.items():
            destination = self.workdir / name
            if destination.is_symlink() and destination.resolve() == source:
                continue
            if destination.exists() or destination.is_symlink():
                destination.unlink()
            destination.symlink_to(source)
        return sources

    def _write_control(self):
        destination = self.workdir / "CONTROL"
        if self.settings.control:
            source = Path(self.settings.control).expanduser()
            if not source.is_absolute():
                source = Path.cwd() / source
            if not source.is_file():
                raise FileNotFoundError(f"FourPhonon CONTROL template not found: {source}")
            shutil.copyfile(source, destination)
            print(f"  - Using CONTROL template: {source.resolve()}")
        else:
            atoms = read(self.structure.poscar)
            destination.write_text(render_control(self.settings, atoms), encoding="utf-8")
            print(f"  - Generated CONTROL: {destination}")
        return destination

    def _run_command(self):
        launcher = [
            token.format(nproc=self.settings.mpi_processes)
            for token in self.settings.mpi_launcher
        ]
        return [*launcher, *shlex.split(str(self.settings.command))]

    def _preferred_solution(self):
        return "rta" if self.settings.solver == "rta" else "iterative"

    def _write_slurm_script(self, command):
        settings = self.settings.parallel
        job_name = settings.get("job_name", "nepkappa-kappa4")
        stdout = self.workdir / "fourphonon-%j.out"
        stderr = self.workdir / "fourphonon-%j.err"
        nodes = int(settings.get("nodes", 1))
        ntasks = int(settings.get("ntasks", self.settings.mpi_processes))
        cpus = int(settings.get("cpus_per_task", self.settings.omp_threads))
        lines = [
            "#!/usr/bin/env bash",
            f"#SBATCH --job-name={job_name}",
            f"#SBATCH --output={stdout}",
            f"#SBATCH --error={stderr}",
            f"#SBATCH --nodes={nodes}",
            f"#SBATCH --ntasks={ntasks}",
            f"#SBATCH --cpus-per-task={cpus}",
        ]
        for key, option in (("time", "time"), ("memory", "mem"), ("partition", "partition"), ("account", "account")):
            if settings.get(key) not in (None, ""):
                lines.append(f"#SBATCH --{option}={settings[key]}")
        lines.extend(f"#SBATCH {item}" for item in settings.get("extra_sbatch", []))
        collect = [
            str(Path(sys.executable).resolve()),
            "-m",
            "nepkappa.fourphonon",
            "--collect",
            str(self.workdir),
            "--preferred",
            self._preferred_solution(),
        ]
        if self.settings.only_harmonic:
            collect.append("--allow-no-kappa")
        lines.extend(
            [
                "",
                "set -euo pipefail",
                f"export OMP_NUM_THREADS=${{SLURM_CPUS_PER_TASK:-{cpus}}}",
                f"export OMP_STACKSIZE={shlex.quote(str(self.settings.omp_stacksize))}",
                *settings.get("preamble", []),
                f"cd {shlex.quote(str(self.workdir))}",
                shlex.join(command),
                shlex.join(collect),
                "",
            ]
        )
        script = self.workdir / "run.sh"
        script.write_text("\n".join(lines), encoding="utf-8")
        script.chmod(0o755)
        print(f"  - FourPhonon Slurm script: {script}")
        return script

    def _write_manifest(self, manifest):
        return RunStateStore(self.manifest_path).write(manifest)


def render_control(config, atoms):
    """Render a non-polar FourPhonon CONTROL from an ASE unit cell."""
    symbols = atoms.get_chemical_symbols()
    elements = list(dict.fromkeys(symbols))
    types = [elements.index(symbol) + 1 for symbol in symbols]
    mesh = config.mesh
    scell = config.supercell
    positions = atoms.get_scaled_positions(wrap=True)
    convergence = config.solver != "rta"
    four_iteration = config.solver == "full-iterative"
    lines = [
        "&allocations",
        f"  nelements={len(elements)},",
        f"  natoms={len(atoms)},",
        f"  ngrid(:)={_numbers(mesh)}",
        "&end",
        "&crystal",
        "  lfactor=0.1,",
    ]
    for index, vector in enumerate(np.asarray(atoms.cell), 1):
        lines.append(f"  lattvec(:,{index})={_numbers(vector)},")
    lines.append("  elements=" + " ".join(f"'{item}'" for item in elements) + ",")
    lines.append(f"  types={_numbers(types)},")
    for index, position in enumerate(positions, 1):
        lines.append(f"  positions(:,{index})={_numbers(position)},")
    lines.extend(
        [
            f"  scell(:)={_numbers(scell)}",
            "&end",
            "&parameters",
        ]
    )
    if len(config.temperatures) == 1:
        lines.append(f"  T={config.temperatures[0]:g},")
    else:
        tmin, tmax, tstep = config.temperatures
        lines.extend(
            [f"  T_min={tmin:g},", f"  T_max={tmax:g},", f"  T_step={tstep:g},"]
        )
    lines.extend(
        [
            f"  scalebroad={config.scalebroad:g},",
            f"  num_sample_process_3ph={config.sample_3ph},",
            f"  num_sample_process_3ph_phase_space={config.sample_3ph_phase_space},",
            f"  num_sample_process_4ph={config.sample_4ph},",
            f"  num_sample_process_4ph_phase_space={config.sample_4ph_phase_space}",
            "&end",
            "&flags",
            f"  nonanalytic={_logical(config.nonanalytic)},",
            f"  convergence={_logical(convergence)},",
            f"  isotopes={_logical(config.isotopes)},",
            "  autoisotopes=.TRUE.,",
            f"  onlyharmonic={_logical(config.only_harmonic)},",
            "  espresso=.FALSE.,",
            "  tdep=.FALSE.,",
            "  four_phonon=.TRUE.,",
            f"  four_phonon_iteration={_logical(four_iteration)}",
            "&end",
            "",
        ]
    )
    return "\n".join(lines)


def _numbers(values):
    return " ".join(f"{value:g}" if isinstance(value, (float, np.floating)) else str(value) for value in values)


def _logical(value):
    return ".TRUE." if value else ".FALSE."


def run_streaming(command, cwd, log_path):
    """Run FourPhonon while mirroring output to a persistent log."""
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = str(environment.get("OMP_NUM_THREADS", 1))
    process = subprocess.Popen(
        command,
        cwd=cwd,
        env=environment,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    with Path(log_path).open("w", encoding="utf-8") as handle:
        for line in process.stdout or []:
            print(line, end="")
            handle.write(line)
    return process.wait()


def read_kappa_tensor(path):
    """Read temperature plus nine tensor components from a FourPhonon output."""
    data = np.loadtxt(path, comments="#", ndmin=2)
    if data.shape[1] < 10:
        raise ValueError(f"Expected at least 10 columns in {path}, got {data.shape[1]}.")
    if not np.all(np.isfinite(data[:, :10])):
        raise ValueError(f"Non-finite FourPhonon conductivity data in {path}.")
    return data[:, :10]


def collect_outputs(workdir, preferred=None, allow_no_kappa=False):
    """Normalize FourPhonon tensor outputs and write a YAML summary."""
    workdir = Path(workdir).resolve()
    solutions = {}
    for name, filename in KAPPA_FILES.items():
        path = workdir / filename
        if not path.is_file() or path.stat().st_size == 0:
            continue
        data = read_kappa_tensor(path)
        normalized = workdir / f"kappa4-{name}.dat"
        np.savetxt(
            normalized,
            data[:, [0, 1, 5, 9]],
            header="T(K) kxx(W/mK) kyy(W/mK) kzz(W/mK)",
            fmt="%.10g",
        )
        solutions[name] = {
            "source": str(path),
            "normalized": str(normalized),
            "rows": int(data.shape[0]),
            "temperatures": data[:, 0].tolist(),
            "diagonal": data[:, [1, 5, 9]].tolist(),
        }
    if not solutions and not allow_no_kappa:
        expected = ", ".join(KAPPA_FILES.values())
        raise FileNotFoundError(
            f"FourPhonon finished without a conductivity tensor ({expected}) in {workdir}."
        )
    primary = preferred if preferred in solutions else next(iter(solutions), None)
    summary = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "primary": primary,
        "solutions": solutions,
        "phase_space": {
            name: str(workdir / name)
            for name in ("BTE.P4", "BTE.P4_total", "BTE.Numprocess_4ph")
            if (workdir / name).is_file()
        },
    }
    (workdir / "fourphonon-summary.yaml").write_text(
        yaml.safe_dump(summary, sort_keys=False), encoding="utf-8"
    )
    manifest = workdir / "submission.yaml"
    if manifest.is_file():
        RunStateStore(manifest).mark_complete()
    print(f"  - FourPhonon summary: {workdir / 'fourphonon-summary.yaml'}")
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description="Collect FourPhonon outputs")
    parser.add_argument("--collect", required=True, help="FourPhonon work directory")
    parser.add_argument("--preferred", choices=["rta", "iterative"], default=None)
    parser.add_argument("--allow-no-kappa", action="store_true")
    args = parser.parse_args(argv)
    collect_outputs(args.collect, args.preferred, args.allow_no_kappa)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
