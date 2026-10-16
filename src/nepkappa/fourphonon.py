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
import time

from ase.io import read
import numpy as np
import yaml

from nepkappa.provenance import file_identity
from nepkappa.coordinates import canonicalize_atoms, COORDINATE_TOLERANCE
from nepkappa.run_state import RunStateStore
from nepkappa.output import activity_status, print_detail
from nepkappa.slurm import submit_script


KAPPA_FILES = {
    "rta": "BTE.KappaTensorVsT_RTA",
    "iterative": "BTE.KappaTensorVsT_CONV",
}
WIGNER_RTA_FILES = {
    "coherence": "BTE.KappaCohTensorVsT_RTA",
    "total": "BTE.KappaTotalTensorVsT_RTA",
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
        if self.settings.compare_3ph:
            if self.settings.control or self.settings.wigner or self.settings.only_harmonic:
                raise ValueError("compare-3ph requires generated CONTROL and ordinary transport.")
            # Never mix an old tensor with a partially completed comparison.
            for directory in (self.workdir, self.workdir / "threephonon"):
                if any((directory / name).exists() for name in KAPPA_FILES.values()):
                    raise FileExistsError(
                        f"Comparison results already exist in {directory}; use a new result-dir."
                    )
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
            "scattering": "3ph+4ph",
            "wigner": self.settings.wigner,
            "submitted": False,
            "command": command,
            "control": file_identity(control),
            "inputs": {name: file_identity(path) for name, path in inputs.items()},
        }
        if self.settings.compare_3ph:
            reference = self.workdir / "threephonon"
            reference.mkdir(parents=True, exist_ok=True)
            reference_inputs = {name: path for name, path in inputs.items()
                                if name != "FORCE_CONSTANTS_4TH"}
            for name, path in reference_inputs.items():
                target = reference / name
                if target.exists() or target.is_symlink():
                    target.unlink()
                target.symlink_to(path)
            reference_control = reference / "CONTROL"
            reference_control.write_text(
                render_control(self.settings, read(self.structure.poscar), four_phonon=False),
                encoding="utf-8",
            )
            reference_manifest = {
                **manifest,
                "scattering": "3ph",
                "solver": "rta" if self.settings.solver == "rta" else "3ph-iterative",
                "control": file_identity(reference_control),
                "inputs": {name: file_identity(path) for name, path in reference_inputs.items()},
            }
            RunStateStore(reference / "submission.yaml").write(reference_manifest)
            manifest["three_phonon_reference"] = "threephonon"
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
        if self.settings.compare_3ph:
            print("  - Reference: three-phonon transport with the same FourPhonon backend")
            self._run_checked(command, self.workdir / "threephonon" / "fourphonon.log",
                              cwd=self.workdir / "threephonon")
            collect_outputs(self.workdir / "threephonon", preferred=self._preferred_solution())
            print("  - Comparison: three-phonon plus four-phonon transport")
        self._run_checked(command, log_path)
        summary = collect_outputs(
            self.workdir,
            preferred=self._preferred_solution(),
            allow_no_kappa=self.settings.only_harmonic,
            wigner=self.settings.wigner,
        )
        manifest["status"] = "complete"
        manifest["summary"] = str(self.workdir / "fourphonon-summary.yaml")
        self._write_manifest(manifest)
        return summary

    def _run_checked(self, command, log_path, *, cwd=None):
        print("  - Running FourPhonon transport calculation")
        print_detail(f"  - Running: {shlex.join(command)}")
        compact = getattr(sys.stdout, "detail_stream", sys.stdout) is not sys.stdout
        started = time.monotonic()
        with activity_status(
            "FourPhonon", enabled=compact and getattr(self.cfg, "progress", True)
        ):
            return_code = run_streaming(command, cwd or self.workdir, log_path)
        if return_code != 0:
            raise RuntimeError(
                f"FourPhonon failed with return code {return_code}; see {log_path}."
            )
        print(f"[Done] FourPhonon finished in {time.monotonic() - started:.2f}s")

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
        if not self.settings.control and self.settings.harmonic_format == "shengbte":
            self._check_coordinate_basis(sources[harmonic_name])
        for name, source in sources.items():
            destination = self.workdir / name
            if destination.is_symlink() and destination.resolve() == source:
                continue
            if destination.exists() or destination.is_symlink():
                destination.unlink()
            destination.symlink_to(source)
        return sources

    def _check_coordinate_basis(self, fc2_path):
        """Reject legacy FC2 indexed in a different unit-cell image.

        Never silently wrap CONTROL while leaving existing IFCs unchanged.
        Explicit external CONTROL files retain their supplied convention.
        """
        poscar = Path(self.structure.poscar)
        if not poscar.is_file():
            return  # The normal structure reader reports a missing POSCAR.
        atoms = read(poscar)
        canonical = canonicalize_atoms(atoms).get_scaled_positions(wrap=False)
        positions = atoms.get_scaled_positions()
        metadata = None
        for name in ("phono3py_disp.yaml", "phono3py.yaml", "phonopy.yaml"):
            candidate = fc2_path.parent / name
            if candidate.is_file():
                data = yaml.safe_load(candidate.read_text()) or {}
                cell = data.get("unit_cell", {})
                points = cell.get("points", [])
                if points:
                    if ([p["symbol"] for p in points] != atoms.get_chemical_symbols()
                            or not np.allclose(cell["lattice"], atoms.cell, atol=1e-8, rtol=0)):
                        raise ValueError(f"FC2 metadata {candidate} does not match structure.poscar.")
                    positions = np.array([p["coordinates"] for p in points])
                    metadata = candidate
                    break
        if not np.allclose(positions, canonical, atol=COORDINATE_TOLERANCE, rtol=0):
            origin = str(metadata) if metadata else str(poscar)
            raise ValueError(
                f"Inconsistent force-constant coordinate basis in {origin}: "
                "existing FC2 may use different lattice images from FC3/FC4. "
                "Regenerate matched FC2/FC3/FC4 with canonical coordinates in a new "
                "result-dir, or explicitly remap all IFCs with their original metadata. "
                "Do not fix this by editing CONTROL alone."
            )

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
        if self.settings.wigner:
            collect.append("--wigner")
        following_steps = []
        plan = list(self.cfg.sections.workflow.steps or [])
        if "kappa4" in plan:
            following_steps = plan[plan.index("kappa4") + 1 :]
        reference_commands = []
        if self.settings.compare_3ph:
            reference = self.workdir / "threephonon"
            reference_collect = list(collect)
            reference_collect[reference_collect.index("--collect") + 1] = str(reference)
            reference_commands = [
                f"cd {shlex.quote(str(reference))}",
                f"{shlex.join(command)} 2>&1 | tee fourphonon.log",
                shlex.join(reference_collect),
            ]
        lines.extend(
            [
                "",
                "set -euo pipefail",
                f"export OMP_NUM_THREADS=${{SLURM_CPUS_PER_TASK:-{cpus}}}",
                f"export OMP_STACKSIZE={shlex.quote(str(self.settings.omp_stacksize))}",
                *settings.get("preamble", []),
                *reference_commands,
                f"cd {shlex.quote(str(self.workdir))}",
                f"{shlex.join(command)} 2>&1 | tee fourphonon.log",
                shlex.join(collect),
                "",
            ]
        )
        if following_steps:
            config_path = getattr(self.cfg, "config_path", None)
            if not config_path:
                raise ValueError(
                    "A deferred FourPhonon stage with following stages requires "
                    "the original YAML config path."
                )
            resume = [
                str(Path(sys.executable).resolve()),
                "-m",
                "nepkappa",
                "run",
                str(Path(config_path).resolve()),
            ]
            lines.extend(
                [f"NEPKAPPA_RESUME_AFTER=kappa4 {shlex.join(resume)}", ""]
            )
        script = self.workdir / "run.sh"
        script.write_text("\n".join(lines), encoding="utf-8")
        script.chmod(0o755)
        print(f"  - FourPhonon Slurm script: {script}")
        return script

    def _write_manifest(self, manifest):
        return RunStateStore(self.manifest_path).write(manifest)


def render_control(config, atoms, *, four_phonon=True):
    """Render a non-polar FourPhonon CONTROL from an ASE unit cell."""
    symbols = atoms.get_chemical_symbols()
    elements = list(dict.fromkeys(symbols))
    types = [elements.index(symbol) + 1 for symbol in symbols]
    mesh = config.mesh
    scell = config.supercell
    positions = canonicalize_atoms(atoms).get_scaled_positions(wrap=False)
    convergence = config.solver != "rta"
    four_iteration = four_phonon and config.solver == "full-iterative"
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
    lines.append(f"  scalebroad={config.scalebroad:g},")
    if not config.wigner:
        lines.extend(
            [
                f"  num_sample_process_3ph={config.sample_3ph},",
                f"  num_sample_process_3ph_phase_space={config.sample_3ph_phase_space},",
                f"  num_sample_process_4ph={config.sample_4ph},",
                f"  num_sample_process_4ph_phase_space={config.sample_4ph_phase_space}",
            ]
        )
    lines.extend(
        [
            "&end",
            "&flags",
            f"  nonanalytic={_logical(config.nonanalytic)},",
            f"  convergence={_logical(convergence)},",
            f"  isotopes={_logical(config.isotopes)},",
            "  autoisotopes=.TRUE.,",
            f"  onlyharmonic={_logical(config.only_harmonic)},",
            "  espresso=.FALSE.,",
        ]
    )
    if not config.wigner:
        lines.append("  tdep=.FALSE.,")
    lines.extend(
        [
            f"  four_phonon={_logical(four_phonon)},",
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
            print_detail(line, end="", flush=True)
            handle.write(line)
            handle.flush()
    return process.wait()


def read_kappa_tensor(path):
    """Read temperature plus nine tensor components from a FourPhonon output."""
    data = np.loadtxt(path, comments="#", ndmin=2)
    if data.shape[1] < 10:
        raise ValueError(f"Expected at least 10 columns in {path}, got {data.shape[1]}.")
    if not np.all(np.isfinite(data[:, :10])):
        raise ValueError(f"Non-finite FourPhonon conductivity data in {path}.")
    return data[:, :10]


def collect_outputs(workdir, preferred=None, allow_no_kappa=False, wigner=False):
    """Normalize FourPhonon tensors and optional Wigner_Park components."""
    workdir = Path(workdir).resolve()
    fresh_after = None
    manifest = workdir / "submission.yaml"
    run = RunStateStore(manifest).read() if manifest.is_file() else {}
    reference_name = run.get("three_phonon_reference")
    scattering = run.get("scattering", "3ph+4ph")
    if wigner or reference_name or scattering == "3ph":
        created_at = run.get("created_at")
        if created_at:
            fresh_after = datetime.fromisoformat(created_at).timestamp() - 2
    reference_summary = None
    if reference_name:
        if reference_name != "threephonon":
            raise ValueError("Invalid FourPhonon reference directory.")
        reference_dir = workdir / reference_name
        reference_run = RunStateStore(reference_dir / "submission.yaml").read()
        if reference_run.get("scattering") != "3ph" or reference_run.get("command") != run.get("command"):
            raise ValueError("Three-phonon reference must use the same FourPhonon command.")
        # Validate the actual staged inputs as well as their recorded identities.
        for case_dir, case in ((workdir, run), (reference_dir, reference_run)):
            for name, identity in {"CONTROL": case["control"], **case["inputs"]}.items():
                if file_identity(case_dir / name).get("sha256") != identity.get("sha256"):
                    raise ValueError(f"Changed comparison input: {case_dir / name}")
        for name, identity in reference_run["inputs"].items():
            if identity.get("sha256") != run["inputs"].get(name, {}).get("sha256"):
                raise ValueError(f"Three-phonon reference uses different {name}.")
        expected = (workdir / "CONTROL").read_text().replace(
            "four_phonon=.TRUE.", "four_phonon=.FALSE."
        ).replace("four_phonon_iteration=.TRUE.", "four_phonon_iteration=.FALSE.")
        if (reference_dir / "CONTROL").read_text() != expected:
            raise ValueError("Comparison CONTROL settings differ beyond four-phonon flags.")
        reference_summary = collect_outputs(reference_dir, preferred=preferred)
    wigner_data = {}
    if wigner:
        # Validate all raw tensors before writing any normalized Wigner results.
        for component, filename in {
            "population": KAPPA_FILES["rta"],
            **WIGNER_RTA_FILES,
        }.items():
            path = workdir / filename
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(
                    f"Missing Wigner_Park output {path}. Set fourphonon.command "
                    "to an executable built from FourPhonon's Wigner_Park branch."
                )
            if fresh_after is not None and path.stat().st_mtime < fresh_after:
                raise ValueError(f"Stale Wigner_Park output: {path}.")
            wigner_data[component] = read_kappa_tensor(path)
        population = wigner_data["population"]
        for component in ("coherence", "total"):
            data = wigner_data[component]
            if data.shape[0] != population.shape[0] or not np.allclose(
                data[:, 0], population[:, 0], rtol=0, atol=1e-6
            ):
                raise ValueError(
                    f"{WIGNER_RTA_FILES[component]} temperatures do not match "
                    f"{KAPPA_FILES['rta']}."
                )
        if not np.allclose(
            wigner_data["total"][:, 1:10],
            population[:, 1:10] + wigner_data["coherence"][:, 1:10],
            rtol=5e-4,
            atol=5e-4,
        ):
            raise ValueError(
                "Wigner_Park total conductivity does not equal its "
                "population plus coherence tensors."
            )
    solutions = {}
    for name, filename in KAPPA_FILES.items():
        if wigner and name != "rta":
            continue
        path = workdir / filename
        if not path.is_file() or path.stat().st_size == 0:
            continue
        if fresh_after is not None and path.stat().st_mtime < fresh_after:
            raise ValueError(f"Stale FourPhonon output: {path}.")
        data = wigner_data["population"] if wigner else read_kappa_tensor(path)
        prefix = "kappa3" if scattering == "3ph" else "kappa4"
        normalized = workdir / f"{prefix}-{name}.dat"
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
        if wigner and name == "rta":
            for component, filename in WIGNER_RTA_FILES.items():
                component_path = workdir / filename
                component_data = wigner_data[component]
                component_normalized = workdir / f"kappa4-{name}-{component}.dat"
                np.savetxt(
                    component_normalized,
                    component_data[:, [0, 1, 5, 9]],
                    header="T(K) kxx(W/mK) kyy(W/mK) kzz(W/mK)",
                    fmt="%.10g",
                )
                solutions[name][component] = {
                    "source": str(component_path),
                    "normalized": str(component_normalized),
                    "diagonal": component_data[:, [1, 5, 9]].tolist(),
                }
    if not solutions and not allow_no_kappa:
        expected = ", ".join(KAPPA_FILES.values())
        raise FileNotFoundError(
            f"FourPhonon finished without a conductivity tensor ({expected}) in {workdir}."
        )
    if (reference_name or scattering == "3ph") and preferred and preferred not in solutions:
        raise FileNotFoundError(f"Missing requested {preferred} conductivity in {workdir}.")
    if reference_summary:
        if set(solutions) != set(reference_summary["solutions"]):
            raise ValueError("Comparison solver outputs do not match.")
        for name, solution in solutions.items():
            if solution["temperatures"] != reference_summary["solutions"][name]["temperatures"]:
                raise ValueError("Comparison temperatures do not match.")
    primary = preferred if preferred in solutions else next(iter(solutions), None)
    summary = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "primary": primary,
        "scattering": scattering,
        "solver": run.get("solver"),
        "solutions": solutions,
        "wigner": {
            "requested": bool(wigner),
            "components_verified": bool(wigner),
            "implementation": "FourPhonon Wigner_Park" if wigner else None,
        },
        "phase_space": {
            name: str(workdir / name)
            for name in ("BTE.P4", "BTE.P4_total", "BTE.Numprocess_4ph")
            if (workdir / name).is_file()
        },
    }
    if reference_summary:
        summary["three_phonon_reference"] = {
            "workdir": "threephonon", "backend": "fourphonon",
            "solutions": reference_summary["solutions"],
            "primary": reference_summary["primary"],
            "matched_inputs_verified": True,
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
    parser.add_argument("--wigner", action="store_true")
    args = parser.parse_args(argv)
    collect_outputs(args.collect, args.preferred, args.allow_no_kappa, args.wigner)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
