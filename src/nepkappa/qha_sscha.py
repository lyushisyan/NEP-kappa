"""Run SSCHA at the temperature-dependent equilibrium volumes from QHA."""

from __future__ import annotations

import copy
from pathlib import Path

import numpy as np
from ase.io import read, write
from phonopy import Phonopy
from phonopy.file_IO import write_FORCE_CONSTANTS
import yaml

from nepkappa.approximations import sscha_approximation
from nepkappa.config_validation import expand_transport_temperatures, qha_sscha_transport_flags, validate_scph
from nepkappa.qha import (
    load_qha_equilibrium_volumes,
    scale_atoms_to_primitive_volume,
    validate_temperature_range,
)
from nepkappa.runtime import ase_to_phonopy
from nepkappa.temperatures import temperature_directory_name, temperature_points
from nepkappa.sscha import (
    load_full_force_constants,
)


class QHASSCHAWorkflow:
    """Couple isotropic QHA expansion to independent Phonopy SSCHA runs."""

    def __init__(self, config, *, execution, runner_factory=None):
        self.cfg = config
        self.execution = execution
        self.runner_factory = runner_factory
        self.output_dir = Path(config.result_dir).expanduser().resolve()
        self.workdir = self.output_dir / "qha-sscha"

    def run(self):
        """Generate a volume-matched structure, FCs, and SSCHA result at each T."""
        error = validate_scph(self.cfg)
        if error:
            raise ValueError(error)
        if not self.cfg.qha_enabled or not self.cfg.scph_enabled:
            raise ValueError(
                "qha-sscha requires both `qha:` and `scph:` YAML sections."
            )
        parallel = self.cfg.sections.force_constants.parallel
        if str(parallel.get("backend", "none")).lower() == "slurm":
            raise ValueError(
                "qha-sscha currently runs its per-temperature force constants "
                "inside one process. Submit `nepkappa run` itself as one Slurm "
                "job, and set force-constant.parallel.backend: none."
            )
        self.three_phonon, self.four_phonon = qha_sscha_transport_flags(self.cfg)
        self.transport_temperatures = (
            expand_transport_temperatures(self.cfg.temps)
            if self.three_phonon or self.four_phonon else []
        )

        qha_temperatures, qha_volumes, multiplicity = load_qha_equilibrium_volumes(
            self.output_dir
        )
        requested = np.asarray(
            temperature_points(self.cfg.scph_temps), dtype=float
        )
        validate_temperature_range(requested, qha_temperatures)
        volumes = np.interp(requested, qha_temperatures, qha_volumes)

        self.workdir.mkdir(parents=True, exist_ok=True)
        base = read(self.cfg.poscar)
        cases = []
        print("\n[SSCHA at QHA volumes] Sequential thermal-expansion approximation")
        print("  - No SSCHA volume optimization; bubble is optional postprocessing")
        for temperature, primitive_volume in zip(requested, volumes):
            cases.append(
                self._run_temperature(
                    base,
                    float(temperature),
                    float(primitive_volume),
                    multiplicity,
                )
            )

        summary = {
            "status": "complete",
            "method": "Phonopy SSCHA at QHA equilibrium volumes",
            "approximation": sscha_approximation(qha_volume=True),
            "qha_volume_source": str(self.output_dir / "volume-temperature.dat"),
            "normalization_multiplicity": multiplicity,
            "bubble_postprocessing": getattr(self.cfg, "scph_bubble", False),
            "transport": {
                "three_phonon": self.three_phonon,
                "four_phonon": self.four_phonon,
                "temperatures_K": self.transport_temperatures,
                "force_constant_treatment": {
                    "fc2": "SSCHA auxiliary harmonic FC2 at T and QHA equilibrium volume",
                    "fc3": "finite-displacement at QHA equilibrium volume",
                    "fc4": (
                        "finite-displacement at QHA equilibrium volume"
                        if self.four_phonon
                        else "not calculated"
                    ),
                },
            },
            "cases": cases,
        }
        path = self.workdir / "qha-sscha-summary.yaml"
        path.write_text(yaml.safe_dump(summary, sort_keys=False), encoding="utf-8")
        print(f"\n[Done] SSCHA-at-QHA-volumes summary written to {path}")
        return summary

    def _run_temperature(self, base, temperature, primitive_volume, multiplicity):
        from nepkappa.qha import QHAWorkflow

        if self.runner_factory is None:
            raise ValueError("A runner_factory is required to execute temperature cases.")

        label = temperature_directory_name(temperature)
        case_dir = self.workdir / label
        case_dir.mkdir(parents=True, exist_ok=True)
        atoms, scale = scale_atoms_to_primitive_volume(
            base, primitive_volume, multiplicity
        )

        case_config = copy.copy(self.cfg)
        case_config.result_dir = str(case_dir)
        case_config.poscar = str((case_dir / "POSCAR_qha_volume").resolve())
        case_config.do_relax = False
        transport_temperatures = getattr(self, "transport_temperatures", None)
        if transport_temperatures is None:
            transport_temperatures = (
                [temperature] if self.three_phonon or self.four_phonon else []
            )
        transport_here = any(
            np.isclose(temperature, point, rtol=1.0e-10, atol=1.0e-8)
            for point in transport_temperatures
        )
        three_phonon_here = self.three_phonon and transport_here
        four_phonon_here = self.four_phonon and transport_here
        needs_fc3 = (
            three_phonon_here or four_phonon_here
            or getattr(self.cfg, "scph_bubble", False)
        )
        case_config.workflow_steps = ["fc2fc3" if needs_fc3 else "fc2"]
        if four_phonon_here:
            case_config.workflow_steps.append("fc4")
            case_config.fc_format = "both"
        case_config.workflow_steps.append("scph")
        case_config.scph_temps = [temperature, temperature, 1.0]
        case_config.scph_initial_fc2 = None
        case_config.scph_transport_fc3 = None
        case_config.scph_transport_metadata = None
        case_config.scph_random_seed = int(self.cfg.scph_random_seed) + int(
            round(temperature)
        )
        case_config.scph_run_transport = three_phonon_here
        case_config.temps = [temperature]

        runner = self.runner_factory(case_config)
        if self.cfg.qha_relax_internal:
            atoms = QHAWorkflow(case_config, workflow=runner.core)._relax_internal(
                atoms, case_dir
            )
        write(
            case_config.poscar,
            atoms,
            format="vasp",
            direct=True,
            vasp5=True,
        )
        print(
            f"\n  - {temperature:g} K: V_QHA={primitive_volume:.8f} "
            f"A^3/primitive cell, isotropic scale={scale:.8f}"
        )
        force_step = "fc2fc3" if needs_fc3 else "fc2"
        outcome = runner.run(force_step)
        if outcome is not None:
            raise RuntimeError(
                "qha-sscha cannot continue from deferred force calculations. "
                "Use serial force generation inside an outer Slurm allocation."
            )
        if four_phonon_here:
            runner.run("fc4")
        runner.run("scph")
        four_phonon_result = None
        if four_phonon_here:
            four_phonon_result = self._run_four_phonon_transport(
                runner, case_config, atoms, temperature, case_dir
            )
        return {
            "temperature_K": temperature,
            "qha_primitive_volume_A3": primitive_volume,
            "unit_cell_volume_A3": float(atoms.get_volume()),
            "isotropic_scale": scale,
            "directory": str(case_dir),
            "three_phonon_transport": three_phonon_here,
            "four_phonon_transport": four_phonon_here,
            "four_phonon_result": four_phonon_result,
        }

    def _run_four_phonon_transport(
        self, runner, case_config, atoms, temperature, case_dir
    ):
        """Run FourPhonon with FC2(T), FC3(V(T)), and FC4(V(T))."""
        temperature_dir = (
            case_dir
            / case_config.scph_workdir
            / temperature_directory_name(temperature)
        )
        fc2_hdf5 = temperature_dir / "fc2.hdf5"
        fc2_text = temperature_dir / "FORCE_CONSTANTS_2ND_SSCHA"
        export_sscha_fc2(
            fc2_hdf5,
            fc2_text,
            atoms,
            case_config.dim_fc2,
        )

        case_config.fp_fc2 = str(fc2_text.resolve())
        case_config.fp_fc3 = str((case_dir / "FORCE_CONSTANTS_3RD").resolve())
        case_config.fp_fc4 = str((case_dir / "FORCE_CONSTANTS_4TH").resolve())
        case_config.fp_temps = [temperature]
        case_config.fp_parallel = {
            **dict(case_config.fp_parallel or {}),
            "backend": "none",
        }
        print(
            f"  - Running FourPhonon with FC2_SSCHA({temperature:g} K), "
            "FC3(V_QHA), and FC4(V_QHA)"
        )
        result = runner.run("kappa4")
        if isinstance(result, dict) and result.get("submitted"):
            raise RuntimeError(
                "qha-sscha FourPhonon transport must run inside the outer job; "
                "nested submission is not supported."
            )
        return result


def export_sscha_fc2(source, destination, atoms, dim_fc2):
    """Export compact/full SSCHA FC2 as a full ShengBTE force-constant file."""
    source = Path(source)
    destination = Path(destination)
    if not source.is_file():
        raise FileNotFoundError(f"SSCHA FC2 not found: {source}")
    phonon = Phonopy(
        ase_to_phonopy(atoms),
        supercell_matrix=np.diag(dim_fc2),
        primitive_matrix="auto",
    )
    full_fc2 = load_full_force_constants(source, phonon)
    write_FORCE_CONSTANTS(full_fc2, filename=str(destination))
    print(f"  - Exported SSCHA FC2 for FourPhonon: {destination}")
    return destination
