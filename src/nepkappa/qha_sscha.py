"""Run SSCHA at the temperature-dependent equilibrium volumes from QHA."""

from __future__ import annotations

import copy
from pathlib import Path

import numpy as np
from ase.io import read, write
from phonopy import Phonopy
from phonopy.file_IO import write_FORCE_CONSTANTS
import yaml

from nepkappa.config import qha_sscha_transport_flags
from nepkappa.runtime import ase_to_phonopy
from nepkappa.sscha import (
    load_full_force_constants,
    temperature_directory_name,
    temperature_points,
)


class QHASSCHAWorkflow:
    """Couple isotropic QHA expansion to independent Phonopy SSCHA runs."""

    def __init__(self, config, *, execution):
        self.cfg = config
        self.execution = execution
        self.output_dir = Path(config.result_dir).expanduser().resolve()
        self.workdir = self.output_dir / "qha-sscha"

    def run(self):
        """Generate a volume-matched structure, FCs, and SSCHA result at each T."""
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
        print("\n[QHA+SSCHA] SSCHA at QHA equilibrium volumes")
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
            "method": "QHA equilibrium volume + Phonopy SSCHA",
            "qha_volume_source": str(self.output_dir / "volume-temperature.dat"),
            "normalization_multiplicity": multiplicity,
            "transport": {
                "three_phonon": self.three_phonon,
                "four_phonon": self.four_phonon,
                "force_constant_treatment": {
                    "fc2": "SSCHA-renormalized at T and QHA equilibrium volume",
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
        print(f"\n[Done] QHA+SSCHA summary written to {path}")
        return summary

    def _run_temperature(self, base, temperature, primitive_volume, multiplicity):
        from nepkappa.application import WorkflowStageRunner
        from nepkappa.qha import QHAWorkflow

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
        case_config.workflow_preset = "custom"
        needs_fc3 = self.three_phonon or self.four_phonon
        case_config.workflow_steps = ["fc2fc3" if needs_fc3 else "fc2"]
        if self.four_phonon:
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
        case_config.scph_run_transport = self.three_phonon
        case_config.temps = [temperature]

        runner = WorkflowStageRunner(case_config, self.execution)
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
                "QHA+SSCHA cannot continue from deferred force calculations. "
                "Use serial force generation inside an outer Slurm allocation."
            )
        if self.four_phonon:
            runner.run("fc4")
        runner.run("scph")
        four_phonon_result = None
        if self.four_phonon:
            four_phonon_result = self._run_four_phonon_transport(
                runner, case_config, atoms, temperature, case_dir
            )
        return {
            "temperature_K": temperature,
            "qha_primitive_volume_A3": primitive_volume,
            "unit_cell_volume_A3": float(atoms.get_volume()),
            "isotropic_scale": scale,
            "directory": str(case_dir),
            "three_phonon_transport": self.three_phonon,
            "four_phonon_transport": self.four_phonon,
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
                "QHA+SSCHA FourPhonon transport must run inside the outer job; "
                "nested submission is not supported."
            )
        return result


def load_qha_equilibrium_volumes(result_dir):
    """Load QHA temperatures, primitive volumes, and normalization metadata."""
    result_dir = Path(result_dir)
    volume_path = result_dir / "volume-temperature.dat"
    summary_path = result_dir / "qha-summary.yaml"
    missing = [path for path in (volume_path, summary_path) if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "QHA+SSCHA needs completed QHA outputs. Missing: "
            + ", ".join(str(path) for path in missing)
        )
    values = np.loadtxt(volume_path, comments="#", ndmin=2)
    if values.shape[1] < 2 or values.shape[0] < 2:
        raise ValueError(f"Invalid QHA volume table: {volume_path}")
    temperatures = np.asarray(values[:, 0], dtype=float)
    volumes = np.asarray(values[:, 1], dtype=float)
    if not np.all(np.isfinite(temperatures)) or not np.all(np.isfinite(volumes)):
        raise ValueError(f"Non-finite values in QHA volume table: {volume_path}")
    if np.any(np.diff(temperatures) <= 0) or np.any(volumes <= 0):
        raise ValueError(
            "QHA volume-temperature data require increasing temperatures and "
            "positive volumes."
        )
    summary = yaml.safe_load(summary_path.read_text(encoding="utf-8")) or {}
    multiplicity = int(summary.get("normalization_multiplicity", 0))
    if multiplicity <= 0:
        raise ValueError(
            f"Missing positive normalization_multiplicity in {summary_path}"
        )
    return temperatures, volumes, multiplicity


def validate_temperature_range(requested, available):
    """Ensure all SSCHA temperatures can be interpolated from the QHA result."""
    if requested.size == 0:
        raise ValueError("No SSCHA temperatures were requested.")
    if requested.min() < available.min() or requested.max() > available.max():
        raise ValueError(
            "scph.temps must lie inside the completed QHA range "
            f"[{available.min():g}, {available.max():g}] K."
        )


def scale_atoms_to_primitive_volume(atoms, primitive_volume, multiplicity):
    """Isotropically scale a cell to a target primitive-cell volume."""
    if primitive_volume <= 0 or multiplicity <= 0:
        raise ValueError("Target volume and normalization multiplicity must be positive.")
    target_volume = float(primitive_volume) * int(multiplicity)
    scale = (target_volume / float(atoms.get_volume())) ** (1.0 / 3.0)
    scaled = atoms.copy()
    scaled.set_cell(atoms.cell * scale, scale_atoms=True)
    return scaled, float(scale)


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
