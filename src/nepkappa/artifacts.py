"""Persistent scientific artifacts and audited force-job cache handling."""

from __future__ import annotations

import json
from pathlib import Path

from ase.io import read, write
import numpy as np


class ForceArtifactStore:
    """Persist displacement structures, forces, and their input identities."""

    def __init__(self, *, vasp_root, generic_root):
        self.vasp_root = Path(vasp_root)
        self.generic_root = Path(generic_root)

    def job_dir(self, label, index, *, calculator):
        """Return the stable directory for one displaced-structure job."""
        root = self.vasp_root if calculator == "vasp" else self.generic_root
        return root / label.lower() / f"{index:05d}"

    @staticmethod
    def validate_forces(forces, atom_count, source):
        """Normalize forces and reject incomplete or corrupt results."""
        array = np.asarray(forces, dtype=float)
        expected_shape = (atom_count, 3)
        if array.shape != expected_shape:
            raise RuntimeError(
                f"Invalid force array from {source}: expected {expected_shape}, "
                f"got {array.shape}."
            )
        if not np.all(np.isfinite(array)):
            raise RuntimeError(f"Non-finite forces found in {source}.")
        return array

    @staticmethod
    def structures_match(first, second, atol=1.0e-8):
        """Return whether two cached force-job structures are equivalent."""
        if len(first) != len(second):
            return False
        if not np.array_equal(
            first.get_atomic_numbers(), second.get_atomic_numbers()
        ):
            return False
        if not np.allclose(
            first.cell.array, second.cell.array, rtol=0.0, atol=atol
        ):
            return False
        return np.allclose(first.positions, second.positions, rtol=0.0, atol=atol)

    @staticmethod
    def matches_inputs(job_dir, input_digest):
        """Return whether a job marker matches the current scientific inputs."""
        marker = Path(job_dir) / "job-input.sha256"
        if not marker.is_file():
            return False
        return marker.read_text(encoding="utf-8").strip() == input_digest

    @staticmethod
    def record_inputs(job_dir, input_digest, payload):
        """Record a fast marker and the auditable payload behind its digest."""
        job_dir = Path(job_dir)
        marker = job_dir / "job-input.sha256"
        marker.parent.mkdir(parents=True, exist_ok=True)
        marker.write_text(input_digest + "\n", encoding="utf-8")
        (job_dir / "job-input.json").write_text(
            json.dumps(
                {"sha256": input_digest, "inputs": payload},
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )

    def stage_structure(self, atoms, label, index, *, calculator):
        """Persist a displaced structure and return its scheduler task record."""
        job_dir = self.job_dir(label, index, calculator=calculator)
        job_dir.mkdir(parents=True, exist_ok=True)
        poscar = job_dir / "POSCAR"
        write(str(poscar), atoms, format="vasp", direct=True, vasp5=True)
        return {"label": label.lower(), "index": index, "poscar": str(poscar)}

    def calculate_or_reuse(
        self,
        atoms,
        label,
        index,
        *,
        calculator,
        fingerprint,
        calculate,
    ):
        """Return audited cached forces or calculate and atomically store them."""
        job_dir = self.job_dir(label, index, calculator=calculator)
        job_dir.mkdir(parents=True, exist_ok=True)
        poscar = job_dir / "POSCAR"
        force_path = job_dir / "forces.npy"
        backend = f"finite-{label.lower()}"

        cached = self._load_matching_structure_cache(
            atoms,
            poscar,
            force_path,
            job_dir,
            backend,
            fingerprint,
            label,
            index,
        )
        if cached is not None:
            return cached

        write(str(poscar), atoms, format="vasp", direct=True, vasp5=True)
        input_digest, input_payload = fingerprint(poscar, backend)
        cached = self._load_matching_digest_cache(
            atoms, force_path, job_dir, input_digest, label, index
        )
        if cached is not None:
            return cached

        forces = self.validate_forces(
            calculate(atoms, label, index),
            len(atoms),
            f"{calculator} {label.upper()} #{index}",
        )
        temporary = force_path.with_suffix(".npy.tmp")
        with temporary.open("wb") as handle:
            np.save(handle, forces, allow_pickle=False)
        temporary.replace(force_path)
        self.record_inputs(job_dir, input_digest, input_payload)
        return forces

    def _load_matching_structure_cache(
        self,
        atoms,
        poscar,
        force_path,
        job_dir,
        backend,
        fingerprint,
        label,
        index,
    ):
        if not (force_path.is_file() and poscar.is_file()):
            return None
        try:
            cached_atoms = read(str(poscar), format="vasp")
            cached_digest, _ = fingerprint(poscar, backend)
            if self.structures_match(atoms, cached_atoms) and self.matches_inputs(
                job_dir, cached_digest
            ):
                forces = self.validate_forces(
                    np.load(force_path, allow_pickle=False), len(atoms), force_path
                )
                print(f"    Reusing {label.upper()} #{index}: {force_path}")
                return forces
        except (OSError, ValueError, RuntimeError) as exc:
            print(f"    Ignoring invalid force cache {force_path}: {exc}")
        return None

    def _load_matching_digest_cache(
        self, atoms, force_path, job_dir, input_digest, label, index
    ):
        if not (
            force_path.is_file()
            and self.matches_inputs(job_dir, input_digest)
        ):
            return None
        try:
            forces = self.validate_forces(
                np.load(force_path, allow_pickle=False), len(atoms), force_path
            )
        except (OSError, ValueError, RuntimeError) as exc:
            print(f"    Ignoring invalid force cache {force_path}: {exc}")
            return None
        print(f"    Reusing {label.upper()} #{index}: {force_path}")
        return forces
