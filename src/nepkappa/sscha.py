"""Stochastic self-consistent harmonic approximation with ASE calculators.

This module follows Phonopy's MLP-SSCHA iteration while evaluating the sampled
supercells through NEP-kappa's reusable ASE calculator backends.  It deliberately
keeps the transport step separate: SSCHA produces temperature-dependent FC2,
and phono3py consumes each FC2 together with the original FC3.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass
from pathlib import Path
import shutil

import h5py
import numpy as np
from phonopy import Phonopy, __version__ as phonopy_version
from phonopy.file_IO import (
    parse_BORN,
    read_force_constants_hdf5,
    write_force_constants_to_hdf5,
)
from phonopy.harmonic.force_constants import (
    compact_fc_to_full_fc,
    full_fc_to_compact_fc,
)
from phonopy.physical_units import get_physical_units
import yaml

from nepkappa.provenance import canonical_data, data_sha256, file_identity
from nepkappa.runtime import ase_to_phonopy, phonopy_to_ase, progress_iter


TRACE_TYPE = "SSCHATrace"


@dataclass
class SSCHATraceData:
    """In-memory data compatible with Phonopy 4.5's SSCHATrace HDF5 schema."""

    temperature: float
    free_energies: list[float]
    errors: list[float]
    potential_energies: list[float]
    harmonic_potential_energies: list[float]
    reference_energy: float
    lattice_lengths: np.ndarray
    force_constants: np.ndarray
    force_constants_history: list[np.ndarray]
    p2s_map: np.ndarray

    @property
    def iterations(self):
        return len(self.free_energies)


class PhonopySSCHAWorkflow:
    """Compute temperature-dependent FC2 using Phonopy sampling and ASE forces."""

    def __init__(self, config, workflow):
        self.cfg = config
        self.settings = config.sections.scph
        self.workflow = workflow
        self.output_dir = Path(config.result_dir).resolve()
        self.workdir = (self.output_dir / self.settings.workdir).resolve()
        if self.output_dir != self.workdir and self.output_dir not in self.workdir.parents:
            raise ValueError("scph.workdir must stay inside output.result-dir.")
        self._calculator_instance = None

    def run(self):
        """Run one independent SSCHA calculation at every configured temperature."""
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.workflow.load_force_constant_structure()
        initial_fc2 = self._initial_fc2_path()
        if not initial_fc2.is_file():
            raise FileNotFoundError(
                f"Initial FC2 not found: {initial_fc2}. Run `nepkappa fc2` first, "
                "or set scph.initial-fc2."
            )
        print("\n[Phonopy SSCHA] Stochastic self-consistent harmonic approximation")
        print(f"  - Phonopy version: {phonopy_version}")
        print(f"  - Initial FC2: {initial_fc2}")
        print(f"  - Supercell: {list(self.cfg.dim_fc2)}")
        print(f"  - Snapshots per iteration: {self.settings.snapshots}")
        print(f"  - Iterations: {self.settings.iterations}")

        results = []
        for temperature in temperature_points(self.settings.temperatures):
            results.append(self._run_temperature(float(temperature), initial_fc2))

        index = {
            "backend": "phonopy-sscha",
            "calculator": self.workflow._calculator_name(),
            "initial_fc2": str(initial_fc2),
            "temperatures": [result["temperature"] for result in results],
            "directories": [result["directory"] for result in results],
            "transport": self.settings.run_transport,
        }
        (self.workdir / "sscha-summary.yaml").write_text(
            yaml.safe_dump(index, sort_keys=False), encoding="utf-8"
        )
        print(f"\n[Done] Phonopy SSCHA results written below {self.workdir}")
        return index

    def _run_temperature(self, temperature, initial_fc2):
        label = temperature_directory_name(temperature)
        temperature_dir = self.workdir / label
        temperature_dir.mkdir(parents=True, exist_ok=True)
        phonon = self._new_phonopy()
        initial_full = load_full_force_constants(initial_fc2, phonon)
        phonon.force_constants = initial_full
        reference_energy = self._evaluate_reference_energy(phonon)
        fingerprint = self._fingerprint(temperature, initial_fc2, phonon)
        trace_path = temperature_dir / "mlpsscha.hdf5"
        state_path = temperature_dir / "sscha-state.yaml"
        trace = self._resume_trace(
            trace_path,
            state_path,
            fingerprint,
            temperature,
            reference_energy,
            phonon,
        )
        if trace.iterations:
            phonon.force_constants = compact_fc_to_full_fc(
                phonon.primitive, trace.force_constants
            )
            print(
                f"\n[SSCHA {temperature:g} K] Resuming after "
                f"{trace.iterations}/{self.settings.iterations} iterations"
            )
        else:
            print(f"\n[SSCHA {temperature:g} K] Starting")

        for iteration in range(trace.iterations + 1, self.settings.iterations + 1):
            print(f"  - Iteration {iteration}/{self.settings.iterations}")
            sampled_fc = np.asarray(phonon.force_constants, dtype=float).copy()
            seed = iteration_seed(self.settings.random_seed, iteration)
            try:
                phonon.generate_displacements(
                    number_of_snapshots=self.settings.snapshots,
                    temperature=temperature,
                    random_seed=seed,
                    cutoff_frequency=self.settings.cutoff_frequency,
                )
            except Exception as exc:
                raise RuntimeError(
                    f"Could not sample the {temperature:g} K SSCHA ensemble. "
                    "The initial/current FC2 may contain unstable modes; inspect the "
                    "harmonic dispersion or adjust scph.cutoff-frequency."
                ) from exc
            displacements = np.asarray(phonon.displacements, dtype=float)
            supercells = phonon.supercells_with_displacements
            forces, energies = self._evaluate_supercells(
                supercells, temperature, iteration
            )
            phonon.force_constants = sampled_fc
            harmonic_fe, potential_e, harmonic_potential_e, error = free_energy_terms(
                phonon,
                sampled_fc,
                displacements,
                energies,
                reference_energy,
                temperature,
                self.settings.sscha_mesh,
            )
            anharmonic = potential_e - harmonic_potential_e
            trace.free_energies.append(harmonic_fe + anharmonic)
            trace.errors.append(error)
            trace.potential_energies.append(potential_e)
            trace.harmonic_potential_energies.append(harmonic_potential_e)
            trace.force_constants_history.append(
                full_fc_to_compact_fc(phonon.primitive, sampled_fc)
            )

            phonon.produce_force_constants(
                forces=forces,
                calculate_full_force_constants=True,
                fc_calculator=self.settings.fc_calculator,
                fc_calculator_options=self.settings.fc_calculator_options,
                show_drift=False,
                fc_calculator_log_level=0,
            )
            trace.force_constants = full_fc_to_compact_fc(
                phonon.primitive, np.asarray(phonon.force_constants, dtype=float)
            )
            write_trace_hdf5(trace, trace_path)
            if self.settings.save_datasets:
                np.savez_compressed(
                    temperature_dir / f"dataset-iter-{iteration:03d}.npz",
                    displacements=displacements,
                    forces=forces,
                    energies=energies,
                )
            print(
                "    F = "
                f"{trace.free_energies[-1] * 1000:.6f} +/- "
                f"{trace.errors[-1] * 1000:.6f} meV/primitive cell"
            )

        averaged_fc = np.mean(
            np.asarray(trace.force_constants_history)[self.settings.transient :], axis=0
        )
        translational_drift = float(np.max(np.abs(np.sum(averaged_fc, axis=1))))
        write_force_constants_to_hdf5(
            averaged_fc,
            filename=str(temperature_dir / "force_constants.hdf5"),
            p2s_map=trace.p2s_map,
            physical_unit="eV/angstrom^2",
        )
        write_force_constants_to_hdf5(
            averaged_fc,
            filename=str(temperature_dir / "fc2.hdf5"),
            p2s_map=trace.p2s_map,
            physical_unit="eV/angstrom^2",
        )
        phonon.force_constants = compact_fc_to_full_fc(phonon.primitive, averaged_fc)
        phonon.save(
            filename=temperature_dir / "phonopy_sscha.yaml",
            settings={"force_constants": True},
        )
        average = averaged_trace_values(trace, self.settings.transient)
        departure = maximum_kept_departure(trace, self.settings.transient)
        summary = {
            "status": "complete",
            "temperature": temperature,
            "snapshots_per_iteration": self.settings.snapshots,
            "iterations": self.settings.iterations,
            "transient": self.settings.transient,
            "kept_iterations": self.settings.iterations - self.settings.transient,
            "free_energy_ev_per_primitive_cell": average[0],
            "free_energy_error_ev_per_primitive_cell": average[1],
            "maximum_kept_free_energy_departure_sigma": departure,
            "reference_energy_ev_per_primitive_cell": trace.reference_energy,
            "maximum_translational_drift_ev_per_angstrom2": translational_drift,
            "trace": trace_path.name,
            "force_constants": "force_constants.hdf5",
            "phono3py_fc2": "fc2.hdf5",
        }
        (temperature_dir / "summary.yaml").write_text(
            yaml.safe_dump(summary, sort_keys=False), encoding="utf-8"
        )
        print(
            f"  - Averaged FC2 translational drift: {translational_drift:.3e} "
            "eV/angstrom^2"
        )
        if self.settings.run_transport:
            self._run_transport(temperature, temperature_dir)
        return {"temperature": temperature, "directory": label}

    def _new_phonopy(self):
        phonon = Phonopy(
            ase_to_phonopy(self.workflow.prim),
            supercell_matrix=np.diag(self.cfg.dim_fc2),
            primitive_matrix="auto",
        )
        if self.settings.born:
            born = Path(self.settings.born).expanduser().resolve()
            if not born.is_file():
                raise FileNotFoundError(f"BORN file not found: {born}")
            phonon.nac_params = parse_BORN(phonon.primitive, filename=born)
        return phonon

    def _calculator(self):
        if self._calculator_instance is None:
            name = self.workflow._calculator_name()
            if name == "nep":
                self._calculator_instance = self.workflow._make_nep_calculator()
            elif name == "vasp":
                raise ValueError(
                    "Phonopy SSCHA requires an ASE calculator and does not use the "
                    "command-driven VASP backend."
                )
            else:
                self._calculator_instance = self.workflow._make_external_backend().calculator()
        return self._calculator_instance

    def _evaluate_reference_energy(self, phonon):
        atoms = phonopy_to_ase(phonon.supercell)
        atoms.calc = self._calculator()
        try:
            energy = float(atoms.get_potential_energy())
        except Exception as exc:
            raise RuntimeError(
                "The selected SSCHA calculator must provide total energies as well as forces."
            ) from exc
        if not np.isfinite(energy):
            raise RuntimeError("The SSCHA calculator returned a non-finite energy.")
        return energy

    def _evaluate_supercells(self, supercells, temperature, iteration):
        forces = []
        energies = []
        valid_supercells = [supercell for supercell in supercells if supercell is not None]
        for index, supercell in enumerate(
            progress_iter(
                valid_supercells,
                enabled=self.workflow.show_progress,
                total=len(valid_supercells),
                desc=f"SSCHA {temperature:g} K iter {iteration}",
                unit="structure",
            ),
            1,
        ):
            atoms = phonopy_to_ase(supercell)
            atoms.calc = self._calculator()
            try:
                energy = float(atoms.get_potential_energy())
                force = np.asarray(atoms.get_forces(), dtype=float)
            except Exception as exc:
                raise RuntimeError(
                    f"SSCHA calculator failed at {temperature:g} K, iteration "
                    f"{iteration}, snapshot {index}."
                ) from exc
            expected = (len(atoms), 3)
            if force.shape != expected or not np.all(np.isfinite(force)):
                raise RuntimeError(
                    f"Invalid SSCHA forces at snapshot {index}: expected {expected}, "
                    f"got {force.shape}."
                )
            if not np.isfinite(energy):
                raise RuntimeError(f"Non-finite SSCHA energy at snapshot {index}.")
            forces.append(force)
            energies.append(energy)
        return np.asarray(forces), np.asarray(energies)

    def _initial_fc2_path(self):
        if self.settings.initial_fc2:
            return Path(self.settings.initial_fc2).expanduser().resolve()
        return self.output_dir / "fc2.hdf5"

    def _fingerprint(self, temperature, initial_fc2, phonon):
        payload = {
            "schema": 1,
            "temperature": temperature,
            "snapshots": self.settings.snapshots,
            "mesh": list(self.settings.sscha_mesh),
            "random_seed": self.settings.random_seed,
            "cutoff_frequency": self.settings.cutoff_frequency,
            "fc_calculator": self.settings.fc_calculator,
            "fc_calculator_options": self.settings.fc_calculator_options,
            "initial_fc2": file_identity(initial_fc2),
            "cell": np.asarray(phonon.unitcell.cell),
            "positions": np.asarray(phonon.unitcell.scaled_positions),
            "symbols": list(phonon.unitcell.symbols),
        }
        if self.workflow._calculator_name() == "nep":
            payload["potential"] = file_identity(self.cfg.nep_model)
        else:
            payload["potential"] = self.workflow._make_external_backend().cache_inputs()
        return data_sha256(canonical_data(payload))

    def _resume_trace(
        self,
        trace_path,
        state_path,
        fingerprint,
        temperature,
        reference_energy,
        phonon,
    ):
        if trace_path.is_file() or state_path.is_file():
            if not trace_path.is_file() or not state_path.is_file():
                raise RuntimeError(
                    f"Incomplete SSCHA checkpoint in {trace_path.parent}; move it aside "
                    "or restore both mlpsscha.hdf5 and sscha-state.yaml."
                )
            state = yaml.safe_load(state_path.read_text(encoding="utf-8")) or {}
            if state.get("fingerprint") != fingerprint:
                raise RuntimeError(
                    f"SSCHA checkpoint settings changed in {trace_path.parent}. Use a "
                    "new scph.workdir or restore the original configuration."
                )
            trace = read_trace_hdf5(trace_path)
            if trace.iterations > self.settings.iterations:
                raise RuntimeError("SSCHA checkpoint has more iterations than requested.")
            return trace

        state_path.write_text(
            yaml.safe_dump({"fingerprint": fingerprint}, sort_keys=False),
            encoding="utf-8",
        )
        n_cell = len(phonon.supercell) // len(phonon.primitive)
        compact = full_fc_to_compact_fc(phonon.primitive, phonon.force_constants)
        trace = SSCHATraceData(
            temperature=temperature,
            free_energies=[],
            errors=[],
            potential_energies=[],
            harmonic_potential_energies=[],
            reference_energy=reference_energy / n_cell,
            lattice_lengths=np.linalg.norm(phonon.unitcell.cell, axis=1),
            force_constants=compact,
            force_constants_history=[],
            p2s_map=np.asarray(phonon.primitive.p2s_map, dtype=int),
        )
        write_trace_hdf5(trace, trace_path)
        return trace

    def _run_transport(self, temperature, temperature_dir):
        metadata = (
            Path(self.settings.transport_metadata).expanduser().resolve()
            if self.settings.transport_metadata
            else self.output_dir / "phono3py_disp.yaml"
        )
        fc3 = (
            Path(self.settings.transport_fc3).expanduser().resolve()
            if self.settings.transport_fc3
            else self.output_dir / "fc3.hdf5"
        )
        missing = [path for path in (metadata, fc3) if not path.is_file()]
        if missing:
            raise FileNotFoundError(
                "SSCHA transport needs the original phono3py_disp.yaml and fc3.hdf5. "
                f"Missing: {', '.join(str(path) for path in missing)}"
            )
        shutil.copy2(metadata, temperature_dir / metadata.name)
        shutil.copy2(fc3, temperature_dir / "fc3.hdf5")
        born_source = (
            Path(self.settings.born).expanduser().resolve()
            if self.settings.born
            else self.output_dir / "BORN"
        )
        if born_source.is_file():
            shutil.copy2(born_source, temperature_dir / "BORN")

        from nepkappa.transport import Phono3pyTransportWorkflow

        transport_config = copy.copy(self.cfg)
        transport_config.temps = [temperature]
        print(f"  - Running phono3py transport with FC2({temperature:g} K)")
        Phono3pyTransportWorkflow(
            transport_config,
            temperature_dir,
            run_command=self.workflow._run_command,
        ).run()


def load_full_force_constants(path, phonon):
    """Load compact or full Phonopy FC2 and normalize it to full form."""
    fc = np.asarray(
        read_force_constants_hdf5(path, p2s_map=phonon.primitive.p2s_map),
        dtype=float,
    )
    n_supercell = len(phonon.supercell)
    if fc.ndim != 4 or fc.shape[-2:] != (3, 3) or fc.shape[1] != n_supercell:
        raise ValueError(
            f"FC2 shape {fc.shape} is incompatible with a {n_supercell}-atom "
            "SSCHA supercell. Check force-constant.dim-fc2."
        )
    if fc.shape[0] == n_supercell:
        return fc
    if fc.shape[0] != len(phonon.primitive):
        raise ValueError(
            f"Compact FC2 first dimension is {fc.shape[0]}, but the primitive cell "
            f"contains {len(phonon.primitive)} atoms."
        )
    return compact_fc_to_full_fc(phonon.primitive, fc)


def free_energy_terms(
    phonon,
    force_constants,
    displacements,
    energies,
    reference_energy,
    temperature,
    mesh,
):
    """Return harmonic FE, sampled potential terms, and statistical error."""
    phonon.force_constants = force_constants
    phonon.run_mesh(mesh=mesh)
    phonon.run_thermal_properties(temperatures=[temperature])
    thermal = phonon.thermal_properties
    if hasattr(thermal, "free_energy"):
        free_energy_kjmol = thermal.free_energy[0]
    else:
        # Phonopy <= 4.4 exposes (T, F, S, Cv) through the nested property;
        # Phonopy 4.5 exposes the individual arrays directly.
        free_energy_kjmol = thermal.thermal_properties[1][0]
    harmonic_fe = float(free_energy_kjmol) / float(get_physical_units().EvTokJmol)
    n_cell = len(phonon.supercell) // len(phonon.primitive)
    n_cart = 3 * force_constants.shape[0]
    phi = force_constants.transpose(0, 2, 1, 3).reshape(n_cart, n_cart)
    u = np.asarray(displacements, dtype=float).reshape(-1, n_cart)
    harmonic_samples = (np.dot(u, phi) * u).sum(axis=1) / 2
    potential_samples = np.asarray(energies, dtype=float) - reference_energy
    anharmonic_samples = (potential_samples - harmonic_samples) / n_cell
    potential = float(np.mean(potential_samples) / n_cell)
    harmonic_potential = float(np.mean(harmonic_samples) / n_cell)
    error = (
        float(np.std(anharmonic_samples, ddof=1) / np.sqrt(len(anharmonic_samples)))
        if len(anharmonic_samples) > 1
        else float("nan")
    )
    return harmonic_fe, potential, harmonic_potential, error


def write_trace_hdf5(trace, filename):
    """Atomically write a Phonopy-compatible SSCHA trace checkpoint."""
    filename = Path(filename)
    temporary = filename.with_suffix(filename.suffix + ".tmp")
    with h5py.File(temporary, "w") as handle:
        handle.attrs["creator"] = "nepkappa"
        handle.attrs["phonopy_version"] = phonopy_version
        handle.attrs["type"] = TRACE_TYPE
        handle.attrs["unit"] = "eV/primitive_cell"
        handle.create_dataset("temperature", data=trace.temperature)
        handle.create_dataset("free_energies", data=np.asarray(trace.free_energies))
        handle.create_dataset("errors", data=np.asarray(trace.errors))
        handle.create_dataset(
            "potential_energies", data=np.asarray(trace.potential_energies)
        )
        handle.create_dataset(
            "harmonic_potential_energies",
            data=np.asarray(trace.harmonic_potential_energies),
        )
        handle.create_dataset("reference_energy", data=trace.reference_energy)
        handle.create_dataset("lattice_lengths", data=trace.lattice_lengths)
        handle.create_dataset("force_constants", data=trace.force_constants)
        handle.create_dataset(
            "force_constants_history", data=np.asarray(trace.force_constants_history)
        )
        handle.create_dataset("p2s_map", data=trace.p2s_map)
    temporary.replace(filename)


def read_trace_hdf5(filename):
    """Read a trace checkpoint written by this module or Phonopy 4.5."""
    with h5py.File(filename, "r") as handle:
        trace_type = handle.attrs.get("type", "")
        if isinstance(trace_type, bytes):
            trace_type = trace_type.decode()
        if trace_type != TRACE_TYPE:
            raise ValueError(f"{filename} is not an {TRACE_TYPE} file.")
        return SSCHATraceData(
            temperature=float(handle["temperature"][()]),
            free_energies=list(np.asarray(handle["free_energies"])),
            errors=list(np.asarray(handle["errors"])),
            potential_energies=list(np.asarray(handle["potential_energies"])),
            harmonic_potential_energies=list(
                np.asarray(handle["harmonic_potential_energies"])
            ),
            reference_energy=float(handle["reference_energy"][()]),
            lattice_lengths=np.asarray(handle["lattice_lengths"]),
            force_constants=np.asarray(handle["force_constants"]),
            force_constants_history=list(
                np.asarray(handle["force_constants_history"])
            ),
            p2s_map=np.asarray(handle["p2s_map"], dtype=int),
        )


def averaged_trace_values(trace, transient):
    """Return free-energy mean and quadrature-combined sampling error."""
    energies = np.asarray(trace.free_energies[transient:], dtype=float)
    errors = np.asarray(trace.errors[transient:], dtype=float)
    return float(np.mean(energies)), float(np.sqrt(np.square(errors).sum()) / len(errors))


def maximum_kept_departure(trace, transient):
    """Return the largest kept free-energy departure in its own standard errors."""
    energies = np.asarray(trace.free_energies, dtype=float)
    errors = np.asarray(trace.errors, dtype=float)
    kept = slice(transient, None)
    mean = float(np.mean(energies[kept]))
    with np.errstate(divide="ignore", invalid="ignore"):
        departure = np.abs((energies[kept] - mean) / errors[kept])
    finite = departure[np.isfinite(departure)]
    return float(np.max(finite)) if finite.size else None


def temperature_points(values):
    """Expand inclusive [minimum, maximum, step] temperatures robustly."""
    minimum, maximum, step = (float(value) for value in values)
    count = int(np.floor((maximum - minimum) / step + 1.0e-10))
    return [minimum + index * step for index in range(count + 1)]


def temperature_directory_name(temperature):
    """Return a sortable, filesystem-safe temperature directory name."""
    if float(temperature).is_integer():
        return f"T{int(temperature):04d}K"
    value = f"{temperature:.6f}".rstrip("0").rstrip(".").replace(".", "p")
    return f"T{value}K"


def iteration_seed(random_seed, iteration):
    """Derive an independent deterministic seed for one iteration."""
    sequence = np.random.SeedSequence([int(random_seed), int(iteration)])
    return int(sequence.generate_state(1, dtype=np.uint32)[0])
