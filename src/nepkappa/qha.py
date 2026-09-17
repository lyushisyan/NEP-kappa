"""Isotropic quasi-harmonic approximation workflow."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import warnings

import numpy as np
import yaml
from ase.optimize import FIRE
from ase.io import read, write
from phonopy import Phonopy
from phonopy.api_qha import PhonopyQHA
from phonopy.file_IO import write_force_constants_to_hdf5

from nepkappa.runtime import ase_to_phonopy, phonopy_to_ase, progress_iter


@dataclass
class QHAPoint:
    """Calculated data at one isotropically scaled cell volume."""

    ratio: float
    volume: float
    electronic_energy: float
    temperatures: np.ndarray
    free_energy: np.ndarray
    heat_capacity: np.ndarray
    entropy: np.ndarray
    normalization_multiplicity: int
    directory: Path


class QHAWorkflow:
    """Calculate harmonic thermodynamics over volume and fit isotropic QHA."""

    def __init__(self, config, *, workflow):
        self.cfg = config
        self.workflow = workflow
        self.output_dir = Path(config.result_dir).resolve()
        self.points_dir = self.output_dir / "qha-points"
        self.points_dir.mkdir(parents=True, exist_ok=True)

    def run(self):
        """Run all volume points and the final equation-of-state fit."""
        if not self.cfg.qha_enabled:
            raise ValueError("The `qha` command requires a `qha:` YAML section.")
        print("\n[QHA] Isotropic quasi-harmonic approximation")
        base = read(self.cfg.poscar)
        points = []
        for index, ratio in enumerate(
            progress_iter(
                self.cfg.qha_volume_ratios,
                enabled=self.workflow.show_progress,
                desc="QHA volume points",
                unit="volume",
            ),
            1,
        ):
            points.append(self._run_volume_point(base, index, float(ratio)))
        self._fit_qha(points)
        print(f"\n[Done] QHA results written to {self.output_dir}")

    def _run_volume_point(self, base, index, ratio):
        label = f"{index:03d}-ratio-{ratio:.6f}"
        point_dir = self.points_dir / label
        point_dir.mkdir(parents=True, exist_ok=True)
        atoms = scale_atoms_to_volume_ratio(base, ratio)
        write(point_dir / "POSCAR_scaled", atoms, format="vasp", direct=True, vasp5=True)
        print(f"\n  - QHA volume point {index}: ratio={ratio:.6f}")

        if self.cfg.qha_relax_internal:
            atoms = self._relax_internal(atoms, point_dir)
        write(point_dir / "POSCAR", atoms, format="vasp", direct=True, vasp5=True)

        electronic_energy = self._static_energy(atoms, point_dir)
        phonon, thermal = self._calculate_fc2_thermal(atoms, point_dir, label)
        multiplicity = len(atoms) / len(phonon.primitive)
        if multiplicity <= 0 or not np.isclose(multiplicity, round(multiplicity)):
            raise RuntimeError(
                "Could not normalize QHA energy and volume to the primitive cell."
            )
        multiplicity = int(round(multiplicity))
        volume = float(atoms.get_volume() / multiplicity)
        energy = float(electronic_energy / multiplicity)
        metadata = {
            "status": "complete",
            "volume_ratio": ratio,
            "unit_cell_atoms": len(atoms),
            "primitive_atoms": len(phonon.primitive),
            "normalization_multiplicity": multiplicity,
            "unit_cell_volume_A3": float(atoms.get_volume()),
            "primitive_volume_A3": volume,
            "unit_cell_energy_eV": float(electronic_energy),
            "primitive_energy_eV": energy,
            "minimum_mesh_frequency_THz": float(
                np.min(phonon.get_mesh_dict()["frequencies"])
            ),
        }
        (point_dir / "point.yaml").write_text(
            yaml.safe_dump(metadata, sort_keys=False), encoding="utf-8"
        )
        return QHAPoint(
            ratio=ratio,
            volume=volume,
            electronic_energy=energy,
            temperatures=np.asarray(thermal["temperatures"], dtype=float),
            free_energy=np.asarray(thermal["free_energy"], dtype=float),
            heat_capacity=np.asarray(thermal["heat_capacity"], dtype=float),
            entropy=np.asarray(thermal["entropy"], dtype=float),
            normalization_multiplicity=multiplicity,
            directory=point_dir,
        )

    def _calculator(self, atoms):
        name = self.workflow._calculator_name()
        if name == "nep":
            atoms.calc = self.workflow._make_nep_calculator()
        elif name == "vasp":
            raise RuntimeError("VASP is executed through its configured command.")
        else:
            atoms.calc = self.workflow._make_external_backend().calculator()
        return atoms.calc

    def _relax_internal(self, atoms, point_dir):
        name = self.workflow._calculator_name()
        relax_dir = point_dir / "relax"
        if name == "vasp":
            overrides = {
                "isif": 2,
                "ibrion": 2,
                "nsw": int(self.cfg.qha_relax_steps),
                "ediffg": -float(self.cfg.qha_relax_fmax),
                **self.cfg.qha_vasp_relax_kwargs,
            }
            self.workflow._write_vasp_run_inputs(
                relax_dir,
                atoms,
                incar_overrides=overrides,
                system="NEP-kappa QHA fixed-cell relaxation",
            )
            ret = self.workflow._run_command(
                self.workflow._vasp_command(), cwd=relax_dir
            )
            if ret != 0:
                raise RuntimeError(
                    f"VASP QHA relaxation failed in {relax_dir} with return code {ret}"
                )
            relaxed = self.workflow._read_vasp_relaxed_structure(relax_dir)
            # Guard against a calculator changing the QHA volume by mistake.
            relaxed.set_cell(atoms.cell, scale_atoms=False)
            return relaxed

        self._calculator(atoms)
        optimizer = FIRE(atoms, logfile=str(relax_dir.with_suffix(".log")))
        optimizer.run(
            fmax=float(self.cfg.qha_relax_fmax),
            steps=int(self.cfg.qha_relax_steps),
        )
        return atoms

    def _static_energy(self, atoms, point_dir):
        name = self.workflow._calculator_name()
        if name == "vasp":
            static_dir = point_dir / "static"
            overrides = {
                "ibrion": -1,
                "nsw": 0,
                **self.cfg.qha_vasp_static_kwargs,
            }
            self.workflow._write_vasp_run_inputs(
                static_dir,
                atoms,
                incar_overrides=overrides,
                system="NEP-kappa QHA static energy",
            )
            ret = self.workflow._run_command(
                self.workflow._vasp_command(), cwd=static_dir
            )
            if ret != 0:
                raise RuntimeError(
                    f"VASP QHA static calculation failed in {static_dir} "
                    f"with return code {ret}"
                )
            energy = self._read_vasp_energy(static_dir)
        else:
            self._calculator(atoms)
            energy = atoms.get_potential_energy()
        if not np.isfinite(energy):
            raise RuntimeError(f"Non-finite QHA static energy in {point_dir}")
        return float(energy)

    @staticmethod
    def _read_vasp_energy(run_dir):
        """Read the final VASP energy from an output carrying calculator results."""
        last_error = None
        for filename in ("vasprun.xml", "OUTCAR"):
            path = run_dir / filename
            if not path.is_file() or path.stat().st_size == 0:
                continue
            try:
                return float(read(str(path), index=-1).get_potential_energy())
            except Exception as exc:
                last_error = exc
        raise RuntimeError(f"Could not read a final VASP energy from {run_dir}") from last_error

    def _calculate_fc2_thermal(self, atoms, point_dir, label):
        dim = self.cfg.qha_dim_fc2 or self.cfg.dim_fc2
        phonon = Phonopy(
            ase_to_phonopy(atoms),
            supercell_matrix=np.diag(dim),
            primitive_matrix="auto",
        )
        phonon.generate_displacements(
            distance=float(self.cfg.qha_displacement_distance)
        )
        supercells = phonon.supercells_with_displacements
        print(f"    FC2 displaced supercells: {len(supercells)}")
        forces = []
        for number, supercell in enumerate(supercells, 1):
            displaced = phonopy_to_ase(supercell)
            if self.workflow._calculator_name() == "vasp":
                force = self.workflow._run_vasp_forces_in_dir(
                    point_dir / "fc2" / f"{number:05d}",
                    displaced,
                    f"qha-{label}",
                    number,
                )
            else:
                force = self.workflow._calculate_forces(
                    displaced, f"qha-{label}", number
                )
            forces.append(force)
        phonon.forces = np.asarray(forces)
        phonon.produce_force_constants()
        phonon.symmetrize_force_constants(show_drift=False)
        write_force_constants_to_hdf5(
            phonon.force_constants,
            filename=str(point_dir / "force_constants.hdf5"),
        )
        phonon.save(
            point_dir / "phonopy_params.yaml",
            settings={"force_constants": True},
        )
        phonon.run_mesh(self.cfg.qha_mesh)
        minimum = float(np.min(phonon.get_mesh_dict()["frequencies"]))
        if minimum < float(self.cfg.qha_imaginary_frequency_tolerance):
            raise RuntimeError(
                f"QHA volume ratio has an unstable mesh frequency: {minimum:.6f} "
                f"THz < {self.cfg.qha_imaginary_frequency_tolerance:.6f} THz "
                f"in {point_dir}"
            )
        tmin, tmax, tstep = self.cfg.qha_temps
        # Phonopy QHA needs two points above the requested output maximum for
        # temperature derivatives such as thermal expansion and Cp.
        phonon.run_thermal_properties(
            t_min=tmin,
            t_max=tmax + 2 * tstep,
            t_step=tstep,
            cutoff_frequency=float(self.cfg.qha_cutoff_frequency),
        )
        phonon.write_yaml_thermal_properties(point_dir / "thermal_properties.yaml")
        return phonon, phonon.get_thermal_properties_dict()

    def _fit_qha(self, points):
        temperatures = points[0].temperatures
        multiplicity = points[0].normalization_multiplicity
        for point in points[1:]:
            if not np.array_equal(point.temperatures, temperatures):
                raise RuntimeError("QHA volume points have inconsistent temperatures.")
            if point.normalization_multiplicity != multiplicity:
                raise RuntimeError(
                    "QHA primitive-cell detection changed between volume points; "
                    "use a symmetry-consistent input structure and relaxation."
                )
        volumes = np.asarray([point.volume for point in points])
        energies = np.asarray([point.electronic_energy for point in points])
        # PhonopyQHA expects temperature-major arrays: (n_temperatures, n_volumes).
        free_energy = np.asarray([point.free_energy for point in points]).T
        cv = np.asarray([point.heat_capacity for point in points]).T
        entropy = np.asarray([point.entropy for point in points]).T
        write_energy_volume(self.output_dir / "e-v.dat", volumes, energies)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            qha = PhonopyQHA(
                volumes=volumes,
                electronic_energies=energies,
                temperatures=temperatures,
                free_energy=free_energy,
                cv=cv,
                entropy=entropy,
                pressure=float(self.cfg.qha_pressure),
                eos=self.cfg.qha_eos,
                t_max=float(self.cfg.qha_temps[1]),
                verbose=self.workflow.show_progress,
            )
        outputs = {
            "volume-temperature.dat": qha.write_volume_temperature,
            "thermal_expansion.dat": qha.write_thermal_expansion,
            "bulk_modulus-temperature.dat": qha.write_bulk_modulus_temperature,
            "gibbs-temperature.dat": qha.write_gibbs_temperature,
            "gruneisen-temperature.dat": qha.write_gruneisen_temperature,
            "helmholtz-volume.dat": qha.write_helmholtz_volume,
        }
        for filename, writer in outputs.items():
            writer(filename=self.output_dir / filename)
        qha.write_heat_capacity_P_polyfit(
            filename=self.output_dir / "Cp-temperature_polyfit.dat",
            filename_ev=self.output_dir / "entropy-volume.dat",
            filename_cvv=self.output_dir / "Cv-volume.dat",
            filename_dsdvt=self.output_dir / "dsdv-temperature.dat",
        )
        qha.plot_pdf_volume_temperature(
            filename=self.output_dir / "volume-temperature.pdf"
        )
        qha.plot_pdf_thermal_expansion(
            filename=self.output_dir / "thermal_expansion.pdf"
        )
        qha.plot_pdf_bulk_modulus_temperature(
            filename=self.output_dir / "bulk_modulus-temperature.pdf"
        )
        qha.plot_pdf_heat_capacity_P_polyfit(
            filename=self.output_dir / "Cp-temperature_polyfit.pdf"
        )
        summary = {
            "status": "complete",
            "eos": self.cfg.qha_eos,
            "pressure_GPa": float(self.cfg.qha_pressure),
            "normalization_multiplicity": multiplicity,
            "points": [
                {
                    "volume_ratio": float(point.ratio),
                    "primitive_volume_A3": float(point.volume),
                    "primitive_energy_eV": float(point.electronic_energy),
                    "directory": str(point.directory),
                }
                for point in points
            ],
        }
        (self.output_dir / "qha-summary.yaml").write_text(
            yaml.safe_dump(summary, sort_keys=False), encoding="utf-8"
        )


def scale_atoms_to_volume_ratio(atoms, ratio):
    """Return a copy isotropically scaled to a target volume ratio."""
    if ratio <= 0:
        raise ValueError("Volume ratio must be positive.")
    scaled = atoms.copy()
    scaled.set_cell(atoms.cell * float(ratio) ** (1.0 / 3.0), scale_atoms=True)
    return scaled


def write_energy_volume(path, volumes, energies):
    """Write the two-column volume-energy input understood by phonopy QHA."""
    values = np.column_stack([volumes, energies])
    np.savetxt(
        path,
        values,
        fmt="%.12f",
        header="primitive_cell_volume_A3 primitive_cell_energy_eV",
    )
