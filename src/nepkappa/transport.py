"""Phono3py three-phonon thermal-transport adapter."""

from __future__ import annotations

import contextlib
from pathlib import Path
import shlex
import shutil
import sys
import time
from types import SimpleNamespace

import phono3py as phono3py_module

from nepkappa.output import activity_status, print_detail
from nepkappa.provenance import file_sha256


FORCE_CONSTANT_FILES = ("phono3py_disp.yaml", "fc2.hdf5", "fc3.hdf5")


class Phono3pyTransportWorkflow:
    """Validate inputs and execute serial or Slurm phono3py transport."""

    def __init__(
        self,
        config,
        output_dir,
        *,
        run_command,
        command_resolver=None,
        needs_fc_flags=None,
    ):
        self.config = config
        self.settings = _kappa_settings(config)
        self.output_dir = Path(output_dir).resolve()
        source = self.settings.force_constants
        self.source_dir = Path(source).expanduser().resolve() if source else self.output_dir
        self.disp_path = self.output_dir / "phono3py_disp.yaml"
        self.fc2_path = self.output_dir / "fc2.hdf5"
        self.fc3_path = self.output_dir / "fc3.hdf5"
        self.run_command = run_command
        self.command_resolver = command_resolver or (
            lambda: resolve_phono3py_command(config)
        )
        self.needs_fc_flags = needs_fc_flags or phono3py_needs_fc_flags

    def run(self):
        """Compute conductivity from existing FC2/FC3 and displacement metadata."""
        print("\n[Step 3] Compute Kappa with Phono3py CLI")
        self._validate_inputs()
        self._prepare_inputs()

        if self.settings.command:
            return self._run_custom_command()

        parallel = self.settings.parallel
        if (
            self.settings.method == "lbte"
            and str(parallel.get("backend", "none")).lower() == "slurm"
        ):
            print("  - Method: distributed LBTE with Slurm")
            from nepkappa.slurm import run_lbte_slurm

            return run_lbte_slurm(
                self.settings,
                self.output_dir,
                self.disp_path,
                self.command_resolver(),
                self.needs_fc_flags(),
                self.run_command,
            )

        return self._run_serial()

    def _validate_inputs(self):
        required = [self.source_dir / name for name in FORCE_CONSTANT_FILES]
        missing = [
            path for path in required if not path.is_file()
        ]
        if missing:
            paths = ", ".join(str(path) for path in missing)
            raise FileNotFoundError(
                f"Missing required phono3py file(s): {paths}. "
                "Run `nepkappa fc2fc3` first or place fc2.hdf5, fc3.hdf5, "
                f"and {self.disp_path.name} in {self.source_dir}, or set "
                "kappa.force-constants to their directory."
            )

    def _prepare_inputs(self):
        """Keep a local copy for the CLI and later plotting; never overwrite inputs."""
        if self.source_dir == self.output_dir:
            return
        names = list(FORCE_CONSTANT_FILES)
        if (self.source_dir / "BORN").is_file():
            names.append("BORN")
        elif (self.output_dir / "BORN").exists() or (self.output_dir / "BORN").is_symlink():
            raise FileExistsError(
                "The result directory contains BORN but the force-constant source does not. "
                "Choose a new output.result-dir to avoid mixing NAC data."
            )
        # Check all destinations before copying any files, including dangling links.
        for name in names:
            source, destination = self.source_dir / name, self.output_dir / name
            if destination.exists() or destination.is_symlink():
                if (
                    not destination.is_file()
                    or destination.is_symlink()
                    or destination.samefile(source)
                    or file_sha256(source) != file_sha256(destination)
                ):
                    raise FileExistsError(
                        f"Conflicting input {destination}; choose a new output.result-dir."
                    )
        self.output_dir.mkdir(parents=True, exist_ok=True)
        for name in names:
            destination = self.output_dir / name
            if not destination.exists():
                shutil.copy2(self.source_dir / name, destination)
        print_detail(f"  - Force-constant source: {self.source_dir}")

    def _run_custom_command(self):
        command = shlex.split(self.settings.command)
        print("  - Method: custom phono3py command")
        self._run_checked(command, label="Custom phono3py command")

    def _run_serial(self):
        if self.settings.method == "lbte":
            method_flags = ["--lbte"]
            print("  - Method: LBTE (Linearized Boltzmann Transport Equation)")
        else:
            method_flags = ["--br", "--nu"]
            print("  - Method: RTA (Relaxation Time Approximation)")

        command = [*self.command_resolver(), self.disp_path.name]
        if self.needs_fc_flags():
            command.extend(["--fc2", "--fc3"])
        command.extend([*method_flags, "--mesh", *(str(v) for v in self.settings.mesh)])

        if self.settings.wigner:
            print("  - Wigner transport: enabled via phono3py SMM19 (--tt smm19)")
            # Phono3py v4 provides this in-tree implementation of the unified
            # Wigner transport formulation. The external phono3py-wte 0.1
            # plugin currently imports APIs removed from released v4 builds.
            command.extend(["--tt", "smm19"])
        if self.settings.isotope:
            print("  - Isotope scattering: enabled")
            command.append("--isotope")
        if self.settings.boundary_mfp is not None:
            print(
                "  - Boundary mean free path: "
                f"{self.settings.boundary_mfp:g} micrometer"
            )
            command.extend(["--boundary-mfp", str(self.settings.boundary_mfp)])
        command.extend(temperature_flags(self.settings.temperatures))

        self._run_checked(command, label="Phono3py")
        print(f"  - Result: {self.output_dir / self.expected_kappa_name()}")

    def _run_checked(self, command, *, label):
        print_detail(f"  - Output directory: {self.output_dir}")
        print_detail(f"  - Running command: {shlex.join(command)}")
        detail_stream = getattr(sys.stdout, "detail_stream", sys.stdout)
        compact = detail_stream is not sys.stdout
        started = time.monotonic()
        with activity_status(
            label, enabled=compact and getattr(self.config, "progress", True)
        ), contextlib.redirect_stdout(detail_stream):
            return_code = self.run_command(command, cwd=self.output_dir)
        if return_code != 0:
            log_stream = getattr(sys.stdout, "log_stream", None)
            log_path = getattr(log_stream, "name", self.output_dir / "run.log")
            raise RuntimeError(
                f"{label} failed with return code {return_code}. See {log_path} for details."
            )
        print(f"[Done] {label} finished in {time.monotonic() - started:.2f}s")

    def expected_kappa_name(self):
        return expected_kappa_name(self.settings.mesh)


def _kappa_settings(config):
    """Return a typed section, with a narrow legacy adapter for callers/tests."""
    if hasattr(config, "sections"):
        return config.sections.kappa
    return SimpleNamespace(
        mesh=tuple(config.mesh),
        temperatures=tuple(config.temps),
        method=config.method,
        command=getattr(config, "kappa_command", None),
        isotope=config.isotope,
        boundary_mfp=config.bfmp,
        wigner=config.wigner,
        parallel=dict(getattr(config, "lbte_parallel", {}) or {}),
        force_constants=getattr(config, "kappa_force_constants", None),
    )


def resolve_phono3py_command(config):
    """Return a robust command for the configured phono3py installation."""
    configured = getattr(config, "phono3py_command", None)
    if configured:
        return shlex.split(str(configured))
    executable = shutil.which("phono3py")
    if executable:
        return [executable]
    sibling = Path(sys.executable).with_name("phono3py")
    return [str(sibling)] if sibling.exists() else ["phono3py"]


def phono3py_needs_fc_flags():
    """Return whether the installed phono3py needs explicit FC file flags."""
    try:
        major = int(str(phono3py_module.__version__).split(".", 1)[0])
    except (AttributeError, TypeError, ValueError):
        return True
    return major < 4


def expected_kappa_name(mesh):
    """Return phono3py's conventional conductivity filename for a mesh."""
    return "kappa-m" + "".join(str(value) for value in mesh) + ".hdf5"


def temperature_flags(temperatures):
    """Return phono3py temperature flags for one value or a range."""
    if len(temperatures) == 3:
        tmin, tmax, tstep = temperatures
        return ["--tmin", str(tmin), "--tmax", str(tmax), "--tstep", str(tstep)]
    return ["--ts", str(temperatures[0])]
