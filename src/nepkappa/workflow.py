# -*- coding: utf-8 -*-
# ==============================================================================
# Author: Shixian Liu, Fei Yin
# Date: 2025
# Description: Core workflow logic for Thermal Conductivity calculations 
#              using NEP (Neuroevolution Potential), HiPhive (optional), 
#              and Phono3py.
# ==============================================================================

import subprocess
import hashlib
import shlex
import shutil
import time
import re
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np

# ASE & Calorine
from ase.io import read, write
from ase import Atoms
from calorine.calculators import CPUNEP
from calorine.tools import relax_structure

# HiPhive
from hiphive import ForceConstants
from hiphive.structure_generation import generate_mc_rattled_structures

# Phonopy / Phono3py
from phonopy.file_IO import write_FORCE_CONSTANTS
from phonopy.harmonic.force_constants import compact_fc_to_full_fc
from phono3py import Phono3py
from phono3py.file_IO import write_fc2_to_hdf5, write_fc3_to_hdf5
from phono3py.phonon3.fc3 import compact_fc3_to_full_fc3

from nepkappa import __version__
from nepkappa.calculators import calculator_backend_from_config
from nepkappa.provenance import (
    canonical_data,
    data_sha256,
    file_identity,
    installed_versions,
)
from nepkappa.runtime import (
    ase_to_phonopy,
    format_duration,
    phonopy_to_ase,
)

# --- Main Logic Class ---
class NEPPhononWorkflow:
    """
    Handler for NEP + (HiPhive/FiniteDisp) + Phono3py Workflow.
    """
    def __init__(self, config):
        """
        Initialize the workflow.
        """
        self.cfg = config
        self.show_progress = getattr(config, "progress", True)
        self.stage_timings = []
        self.output_dir = Path(getattr(config, "result_dir", "result")).resolve()
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.fc2_path = self.output_dir / "fc2.hdf5"
        self.fc3_path = self.output_dir / "fc3.hdf5"
        self.shengbte_fc2_path = self.output_dir / "FORCE_CONSTANTS_2ND"
        self.shengbte_fc3_path = self.output_dir / "FORCE_CONSTANTS_3RD"
        self.fc4_path = self.output_dir / "FORCE_CONSTANTS_4TH"
        self.disp_path = self.output_dir / "phono3py_disp.yaml"
        self.relaxed_poscar_path = self.output_dir / "POSCAR_relaxed"
        self.hiphive_model_path = self.output_dir / "hiphive_model.fcp"
        self.vasp_root = self.output_dir / getattr(config, "vasp_workdir", "vasp-runs")
        self.fc3_root = self.output_dir / getattr(
            config, "fc3_workdir", "fc3-thirdorder-runs"
        )
        self.fc4_root = self.output_dir / getattr(config, "fc4_workdir", "fc4-runs")
        self.vasp_relax_root = self.output_dir / getattr(
            config, "vasp_relax_workdir", "vasp-relax"
        )
        self._nep_calc = None
        self._external_backend = None
        self._vasp_backend = None
        self._force_store = None
        self._cache_software_versions = None
        self.prim = None 

    def _run_timed_stage(self, label, func):
        """Run a workflow stage and record its elapsed time."""
        start = time.time()
        try:
            return func()
        finally:
            elapsed = time.time() - start
            self.stage_timings.append((label, elapsed))
            print(f"  - {label} elapsed time: {format_duration(elapsed)}")

    def _print_timing_summary(self):
        if not self.stage_timings:
            return
        print("\n[Timing Summary]")
        for label, elapsed in self.stage_timings:
            print(f"  - {label:<24}: {format_duration(elapsed)}")

    def _run_command(self, cmd, cwd=None):
        """Run a subprocess while forwarding output to terminal and run.log."""
        process = subprocess.Popen(
            cmd,
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        if process.stdout is not None:
            for line in process.stdout:
                print(line, end="")
        return process.wait()

    def _run_command_with_input(self, cmd, input_text, cwd=None):
        """Run a subprocess with stdin while forwarding output to terminal/log."""
        process = subprocess.Popen(
            cmd,
            cwd=cwd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        output, _ = process.communicate(input_text)
        if output:
            print(output, end="")
        return process.returncode

    def _phono3py_command(self):
        """Return a robust phono3py executable command."""
        from nepkappa.transport import resolve_phono3py_command

        return resolve_phono3py_command(self.cfg)

    def _fourthorder_command(self):
        """Return a robust Fourthorder executable command."""
        configured = getattr(self.cfg, "fourthorder_command", None)
        if configured:
            return shlex.split(str(configured))

        for executable_name in ("Fourthorder_vasp.py", "fourthorder_vasp.py"):
            executable = shutil.which(executable_name)
            if executable:
                return [executable]

        return ["Fourthorder_vasp.py"]

    def _thirdorder_command(self):
        """Return a robust thirdorder executable command."""
        configured = getattr(self.cfg, "thirdorder_command", None)
        if configured:
            return shlex.split(str(configured))

        executable = shutil.which("thirdorder_vasp.py")
        if executable:
            return [executable]

        return ["thirdorder_vasp.py"]

    def _format_fourthorder_value(self, value):
        """Format Fourthorder numeric arguments without unnecessary .0 suffixes."""
        value = float(value)
        if value.is_integer():
            return str(int(value))
        return f"{value:g}"

    def _fourthorder_args(self, mode):
        """Build Fourthorder_vasp.py sow/reap command arguments."""
        nx, ny, nz = self.cfg.dim_fc4
        return [
            *self._fourthorder_command(),
            mode,
            str(nx),
            str(ny),
            str(nz),
            self._format_fourthorder_value(self.cfg.cutoff_fc4),
        ]

    def _thirdorder_args(self, mode):
        """Build thirdorder_vasp.py sow/reap command arguments."""
        nx, ny, nz = self.cfg.dim_fc3
        return [
            *self._thirdorder_command(),
            mode,
            str(nx),
            str(ny),
            str(nz),
            self._format_fourthorder_value(self.cfg.cutoff_fc3),
        ]

    def _fc3_backend(self):
        """Return the selected FC3 generation backend."""
        return getattr(self.cfg, "fc3_backend", "phono3py")

    def _fc_format(self):
        """Return requested FC2/FC3 export format."""
        return getattr(self.cfg, "fc_format", "phono3py")

    def _write_phono3py_fc2(self, fc2, p2s_map=None):
        """Write FC2 in phono3py/phonopy HDF5 format."""
        write_fc2_to_hdf5(fc2, filename=str(self.fc2_path), p2s_map=p2s_map)
        print(f"  - Generated: {self.fc2_path}")

    def _write_phono3py_fc3(self, fc3, p2s_map=None, fc3_nonzero_indices=None):
        """Write FC3 in phono3py HDF5 format."""
        write_fc3_to_hdf5(
            fc3,
            fc3_nonzero_indices=fc3_nonzero_indices,
            filename=str(self.fc3_path),
            p2s_map=p2s_map,
        )
        print(f"  - Generated: {self.fc3_path}")

    def _write_shengbte_fc2(self, fc2_full):
        """Write FC2 in ShengBTE/phonopy text format."""
        write_FORCE_CONSTANTS(fc2_full, filename=str(self.shengbte_fc2_path))
        print(f"  - Generated: {self.shengbte_fc2_path}")

    def _write_shengbte_fc3(self, fc3_full, supercell, primitive):
        """Write FC3 in ShengBTE text format."""
        fcs = ForceConstants.from_arrays(supercell, fc3_array=fc3_full)
        fcs.write_to_shengBTE(str(self.shengbte_fc3_path), primitive)
        print(f"  - Generated: {self.shengbte_fc3_path}")

    def _write_fc2_exports(
        self,
        fc2,
        *,
        p2s_map=None,
        fc2_full=None,
        phonopy_primitive=None,
    ):
        """Write FC2 in the requested output format(s)."""
        fc_format = self._fc_format()
        if fc_format in {"phono3py", "both"}:
            self._write_phono3py_fc2(fc2, p2s_map=p2s_map)
        if fc_format in {"shengbte", "both"}:
            if fc2_full is None:
                if p2s_map is not None and phonopy_primitive is not None:
                    fc2_full = compact_fc_to_full_fc(phonopy_primitive, fc2)
                else:
                    fc2_full = fc2
            self._write_shengbte_fc2(fc2_full)

    def _write_fc3_exports(
        self,
        fc3,
        *,
        p2s_map=None,
        fc3_nonzero_indices=None,
        fc3_full=None,
        supercell_atoms=None,
        primitive_atoms=None,
        phono3py_primitive=None,
    ):
        """Write FC3 in the requested output format(s)."""
        fc_format = self._fc_format()
        if fc_format in {"phono3py", "both"}:
            self._write_phono3py_fc3(
                fc3,
                p2s_map=p2s_map,
                fc3_nonzero_indices=fc3_nonzero_indices,
            )
        if fc_format in {"shengbte", "both"}:
            if supercell_atoms is None or primitive_atoms is None:
                raise ValueError("ShengBTE FC3 export requires supercell and primitive")
            if fc3_full is None:
                if p2s_map is not None:
                    if phono3py_primitive is None:
                        raise ValueError(
                            "Compact ShengBTE FC3 export requires phono3py primitive"
                        )
                    fc3_full = compact_fc3_to_full_fc3(phono3py_primitive, fc3)
                else:
                    fc3_full = fc3
            self._write_shengbte_fc3(fc3_full, supercell_atoms, primitive_atoms)

    def _phono3py_needs_fc_flags(self):
        """Return True when the phono3py CLI needs explicit --fc2/--fc3 flags."""
        from nepkappa.transport import phono3py_needs_fc_flags

        return phono3py_needs_fc_flags()

    def _make_phono3py(self):
        """Create a Phono3py object with the workflow's structure settings."""
        return Phono3py(
            ase_to_phonopy(self.prim),
            supercell_matrix=self.cfg.dim_fc3,
            phonon_supercell_matrix=self.cfg.dim_fc2,
            primitive_matrix="auto",
        )

    def _save_phono3py_metadata(self):
        """Write the phono3py yaml file needed by phono3py v4 calculation."""
        ph3 = self._make_phono3py()
        ph3.save(str(self.disp_path))
        print(f"  - Phono3py metadata saved to {self.disp_path}")

    def _calculator_name(self):
        return getattr(self.cfg, "calculator", "nep").lower()

    def _force_desc(self, label):
        return f"{self._calculator_name().upper()} forces ({label})"

    def _make_nep_calculator(self):
        if self._nep_calc is None:
            self._nep_calc = CPUNEP(self.cfg.nep_model)
        return self._nep_calc

    def _make_external_backend(self):
        """Return the configured generic ASE/plugin calculator backend."""
        if self._external_backend is None:
            self._external_backend = calculator_backend_from_config(self.cfg)
        return self._external_backend

    def _make_vasp_backend(self):
        """Create and cache the extracted VASP adapter."""
        if self._vasp_backend is None:
            from nepkappa.adapters.vasp import VaspBackend

            self._vasp_backend = VaspBackend(
                self.cfg,
                output_dir=self.output_dir,
                force_root=self.vasp_root,
                relax_root=self.vasp_relax_root,
                relaxed_structure=self.relaxed_poscar_path,
                run_command=lambda command, cwd=None: self._run_command(
                    command, cwd=cwd
                ),
            )
        return self._vasp_backend

    def _make_force_store(self):
        """Create and cache the force-job artifact store."""
        if self._force_store is None:
            from nepkappa.artifacts import ForceArtifactStore

            self._force_store = ForceArtifactStore(
                vasp_root=self.vasp_root,
                generic_root=self.output_dir / "force-jobs",
            )
        return self._force_store

    def _vasp_command(self):
        """Compatibility proxy for the VASP adapter command resolver."""
        return self._make_vasp_backend().command()

    def _detect_vasp_path(self):
        """Compatibility proxy for VASP executable discovery."""
        return self._make_vasp_backend().detect_path()

    def _resolve_potcar_path(self, atoms):
        """Compatibility proxy for POTCAR resolution."""
        return self._make_vasp_backend().resolve_potcar(atoms)

    def _detect_potcar_path(self, symbols):
        """Compatibility proxy for POTCAR discovery."""
        return self._make_vasp_backend().detect_potcar(symbols)

    def _assemble_potcar_from_library(self, library, symbols):
        """Compatibility proxy for multi-species POTCAR assembly."""
        return self._make_vasp_backend().assemble_potcar(library, symbols)

    def _find_potcar_for_symbol(self, library, symbol):
        return self._make_vasp_backend().find_potcar(library, symbol)

    def _write_vasp_run_inputs(self, run_dir, atoms, incar_overrides=None, system=None):
        """Compatibility proxy for VASP input generation."""
        return self._make_vasp_backend().write_run_inputs(
            run_dir,
            atoms,
            incar_overrides=incar_overrides,
            system=system,
        )

    def _run_vasp_forces(self, atoms, label, index):
        """Compatibility proxy for one VASP force calculation."""
        return self._make_vasp_backend().run_forces(atoms, label, index)

    def _run_vasp_forces_in_dir(self, run_dir, atoms, label, index):
        """Compatibility proxy for a VASP force calculation directory."""
        return self._make_vasp_backend().run_forces_in_dir(
            run_dir, atoms, label, index
        )

    def _run_vasp_fc4_job(self, atoms, job_dir, index):
        """Compatibility proxy for a Fourthorder VASP job."""
        return self._make_vasp_backend().run_order_job(
            atoms, job_dir, "fc4", consumer="Fourthorder reap"
        )

    def _run_vasp_thirdorder_job(self, atoms, job_dir, index):
        """Compatibility proxy for a thirdorder VASP job."""
        return self._make_vasp_backend().run_order_job(
            atoms, job_dir, "fc3-thirdorder", consumer="thirdorder reap"
        )

    def _run_vasp_force_xml_job(self, atoms, job_dir, label):
        """Compatibility proxy for an externally consumed VASP XML job."""
        consumer = "Fourthorder reap" if label == "fc4" else "thirdorder reap"
        return self._make_vasp_backend().run_order_job(
            atoms, job_dir, label, consumer=consumer
        )

    def _write_vasp_force_xml(self, filename, forces):
        """Write the minimal VASP XML force block parsed by *order scripts."""
        force_array = np.asarray(forces, dtype="double")
        if force_array.ndim != 2 or force_array.shape[1] != 3:
            raise ValueError(
                "forces must have shape (number_of_atoms, 3)"
            )
        if not np.all(np.isfinite(force_array)):
            raise ValueError("forces contain non-finite values")

        root = ET.Element("modeling")
        calculation = ET.SubElement(root, "calculation")
        varray = ET.SubElement(calculation, "varray", {"name": "forces"})
        for force in force_array:
            vector = ET.SubElement(varray, "v")
            vector.text = " ".join(f"{component:.16e}" for component in force)

        filename = Path(filename)
        filename.parent.mkdir(parents=True, exist_ok=True)
        ET.ElementTree(root).write(filename, encoding="utf-8", xml_declaration=True)
        return filename

    def _run_nep_thirdorder_job(self, atoms, job_dir, index):
        """Evaluate NEP forces and expose them through thirdorder's VASP parser."""
        atoms.calc = self._make_nep_calculator()
        forces = atoms.get_forces()
        vasprun_path = job_dir / "vasprun.xml"
        return self._write_vasp_force_xml(vasprun_path, forces)

    def _run_nep_fc4_job(self, atoms, job_dir, index):
        """Evaluate NEP forces and expose them through Fourthorder's parser."""
        atoms.calc = self._make_nep_calculator()
        forces = atoms.get_forces()
        vasprun_path = job_dir / "vasprun.xml"
        return self._write_vasp_force_xml(vasprun_path, forces)

    def _run_external_order_job(self, atoms, job_dir):
        """Expose external ASE calculator forces through the order-script XML API."""
        forces = self._make_external_backend().calculate_forces(atoms)
        return self._write_vasp_force_xml(job_dir / "vasprun.xml", forces)

    def _combined_vasp_kwargs(self, overrides=None, system=None):
        """Compatibility proxy for normalized VASP parameters."""
        return self._make_vasp_backend().combined_kwargs(overrides, system)

    def _write_vasp_incar(self, path, kwargs=None):
        return self._make_vasp_backend().write_incar(path, kwargs)

    def _write_vasp_kpoints(self, path, kwargs=None):
        return self._make_vasp_backend().write_kpoints(path, kwargs)

    def _format_vasp_value(self, value):
        return self._make_vasp_backend().format_value(value)

    def _read_vasp_forces(self, run_dir):
        return self._make_vasp_backend().read_forces(run_dir)

    def _default_vasp_relax_stages(self):
        return self._make_vasp_backend().default_relax_stages()

    def _vasp_relax_stages(self):
        return self._make_vasp_backend().relax_stages()

    def _read_vasp_relaxed_structure(self, run_dir):
        return self._make_vasp_backend().read_relaxed_structure(run_dir)

    def _run_vasp_relax_stage(self, atoms, stage_name, stage_index, overrides):
        return self._make_vasp_backend().run_relax_stage(
            atoms, stage_name, stage_index, overrides
        )

    def _run_vasp_relaxation(self, atoms):
        return self._make_vasp_backend().relax(atoms)

    def _calculate_forces(self, atoms, label, index):
        """Calculate forces with the configured backend."""
        calculator = self._calculator_name()
        if calculator == "nep":
            atoms.calc = self._make_nep_calculator()
            return atoms.get_forces()
        elif calculator == "vasp":
            return self._run_vasp_forces(atoms, label, index)
        else:
            return self._make_external_backend().calculate_forces(atoms)

    def _force_job_dir(self, label, index):
        """Return the stable work directory for one displaced structure."""
        return self._make_force_store().job_dir(
            label,
            index,
            calculator=self._calculator_name(),
        )

    def _validated_force_array(self, forces, atom_count, source):
        """Normalize forces and reject incomplete or corrupt cache entries."""
        return self._make_force_store().validate_forces(
            forces, atom_count, source
        )

    def _force_structures_match(self, first, second, atol=1.0e-8):
        """Return whether two force-job structures are numerically equivalent."""
        return self._make_force_store().structures_match(first, second, atol)

    def _calculate_forces_cached(self, atoms, label, index):
        """Calculate or reuse forces for one structure using audited inputs."""
        return self._make_force_store().calculate_or_reuse(
            atoms,
            label,
            index,
            calculator=self._calculator_name(),
            fingerprint=self._job_input_fingerprint,
            calculate=self._calculate_forces,
        )

    def _stage_force_job(self, atoms, label, index):
        """Persist one displaced structure and return its Slurm task record."""
        return self._make_force_store().stage_structure(
            atoms,
            label,
            index,
            calculator=self._calculator_name(),
        )

    def _prepare_force_slurm_jobs(self, include_fc3=True):
        """Generate all structures needed by a native FC or HiPhive route."""
        cfg = self.cfg
        # The collection job reconstructs displacements from POSCAR_relaxed.
        # Round-trip through that persisted structure before staging the array
        # jobs so both processes hash byte-identical displaced POSCAR files.
        if self.relaxed_poscar_path.is_file():
            self.prim = read(str(self.relaxed_poscar_path))
        jobs = []
        if cfg.use_hiphive:
            np.random.seed(42)
            self._save_phono3py_metadata()
            dimensions = cfg.dim_fc3 if include_fc3 else cfg.dim_fc2
            ideal = self.prim.repeat(tuple(dimensions))
            structures = generate_mc_rattled_structures(
                ideal, cfg.n_structures, cfg.rattle_std, cfg.min_dist
            )
            for index, atoms in enumerate(structures, 1):
                jobs.append(self._stage_force_job(atoms, "hiphive", index))
            return jobs

        ph3 = self._make_phono3py()
        ph3.generate_fc2_displacements()
        ph3.save(str(self.disp_path))
        for index, supercell in enumerate(ph3.phonon_supercells_with_displacements, 1):
            jobs.append(
                self._stage_force_job(phonopy_to_ase(supercell), "fc2", index)
            )
        if include_fc3:
            pair_cutoff = getattr(self.cfg, "pair_cutoff_fc3", None)
            if pair_cutoff is None:
                ph3.generate_displacements()
            else:
                ph3.generate_displacements(cutoff_pair_distance=pair_cutoff)
            ph3.save(str(self.disp_path))
            for index, supercell in enumerate(ph3.supercells_with_displacements, 1):
                jobs.append(
                    self._stage_force_job(phonopy_to_ase(supercell), "fc3", index)
                )
        return jobs

    def _submit_force_slurm(self, include_fc3=True):
        """Prepare displaced structures and dispatch their force calculations."""
        from nepkappa.slurm import run_force_slurm

        jobs = self._prepare_force_slurm_jobs(include_fc3=include_fc3)
        return run_force_slurm(self.cfg, self.output_dir, jobs)

    def relax_structure_stage(self):
        """Load and relax through the extracted structure stage."""
        from nepkappa.stages.structure import StructureRelaxationStage

        return StructureRelaxationStage(self, ase_relax=relax_structure).run()

    def load_force_constant_structure(self):
        """Load the structure used by the force-constant stage."""
        if self.prim is not None:
            print(f"  - Using in-memory structure from the relax stage")
            return

        if self.cfg.do_relax:
            if not self.relaxed_poscar_path.exists():
                raise FileNotFoundError(
                    f"Relaxed structure not found: {self.relaxed_poscar_path}. "
                    "Run `nepkappa relax` first, or use `nepkappa run` for the full workflow."
                )
            print(f"  - Reading relaxed structure from {self.relaxed_poscar_path}")
            self.prim = read(str(self.relaxed_poscar_path))
        else:
            print(f"  - Reading input structure from {self.cfg.poscar}")
            self.prim = read(self.cfg.poscar)

    def run_hiphive_fitting(self, include_fc3=True):
        """Fit force constants through the extracted HiPhive stage."""
        from nepkappa.stages.force_constants import HiPhiveStage

        return HiPhiveStage(self, include_fc3=include_fc3).run()

    def run_finite_disp_fitting(self, include_fc3=True):
        """Generate FC2/FC3 through the extracted finite-displacement stage."""
        from nepkappa.stages.force_constants import FiniteDisplacementStage

        return FiniteDisplacementStage(self, include_fc3=include_fc3).run()

    def _displacement_pattern_sort_key(self, path):
        """Sort thirdorder/Fourthorder pattern files by embedded numbers."""
        return [
            int(part) if part.isdigit() else part
            for part in re.split(r"(\d+)", path.name)
        ]

    def _clear_displacement_patterns(self, root, pattern):
        """Remove generated sow patterns so changed settings cannot leave stale files."""
        for path in root.glob(pattern):
            if path.is_file():
                path.unlink()

    def _pattern_digest(self, pattern):
        """Return a stable fingerprint for one displaced structure."""
        return hashlib.sha256(Path(pattern).read_bytes()).hexdigest()

    def _command_cache_identity(self, command):
        """Describe a configured external command and its executable, if resolvable."""
        tokens = list(command)
        executable = None
        for token in reversed(tokens):
            candidate = Path(token).expanduser()
            if candidate.is_file():
                executable = candidate
                break
            resolved = shutil.which(token)
            if resolved:
                executable = Path(resolved)
                break
        return {
            "tokens": tokens,
            "executable": file_identity(executable) if executable else None,
        }

    def _potcar_cache_inputs(self, atoms):
        """Return identities of the exact POTCAR sources selected for a VASP job."""
        configured = getattr(self.cfg, "potcar_path", None)
        symbols = list(dict.fromkeys(atoms.get_chemical_symbols()))
        sources = []
        if configured:
            path = Path(configured).expanduser()
            if path.is_file():
                sources = [path]
            elif path.is_dir():
                sources = [self._find_potcar_for_symbol(path, symbol) for symbol in symbols]
        else:
            detected = self._detect_potcar_path(symbols)
            if detected is not None:
                _, detected_sources = detected
                sources = list(detected_sources)
        return [
            file_identity(source) if source is not None else {"exists": False}
            for source in sources
        ] or [{"path": str(configured) if configured else None, "exists": False}]

    def _job_input_payload(self, pattern, backend):
        """Build all scientifically relevant inputs used to calculate one force job."""
        calculator = self._calculator_name()
        payload = {
            "schema_version": 2,
            "code": {
                "nepkappa_version": __version__,
                "workflow_sha256": file_identity(__file__)["sha256"],
            },
            "backend": backend,
            "displacement_sha256": self._pattern_digest(pattern),
            "calculator": calculator,
        }
        if self._cache_software_versions is None:
            self._cache_software_versions = installed_versions()
        payload["software"] = self._cache_software_versions
        if backend == "thirdorder":
            payload["order_settings"] = {
                "command": self._command_cache_identity(self._thirdorder_command()),
                "dimension": list(self.cfg.dim_fc3),
                "cutoff": float(self.cfg.cutoff_fc3),
            }
        elif backend == "fourthorder":
            payload["order_settings"] = {
                "command": self._command_cache_identity(self._fourthorder_command()),
                "dimension": list(self.cfg.dim_fc4),
                "cutoff": float(self.cfg.cutoff_fc4),
            }

        if calculator == "nep":
            model = getattr(self.cfg, "nep_model", None)
            payload["calculator_inputs"] = {
                "model": file_identity(model) if model else None,
            }
        elif calculator == "vasp":
            command = getattr(self.cfg, "vasp_command", None)
            if command:
                command = shlex.split(str(command))
            elif getattr(self.cfg, "vasp_path", None):
                command = [str(self.cfg.vasp_path)]
            else:
                detected = self._detect_vasp_path()
                command = [str(detected)] if detected else ["vasp_std"]
            atoms = read(str(pattern), format="vasp")
            payload["calculator_inputs"] = {
                "command": self._command_cache_identity(command),
                "incar": canonical_data(self._combined_vasp_kwargs()),
                "potcar_sources": self._potcar_cache_inputs(atoms),
            }
        else:
            payload["calculator_inputs"] = (
                self._make_external_backend().cache_inputs()
            )
        return canonical_data(payload)

    def _job_input_fingerprint(self, pattern, backend):
        """Return the full cache digest and its auditable input payload."""
        payload = self._job_input_payload(pattern, backend)
        return data_sha256(payload), payload

    def _job_matches_inputs(self, job_dir, input_digest):
        return self._make_force_store().matches_inputs(job_dir, input_digest)

    def _record_job_inputs(self, job_dir, input_digest, payload):
        """Record both a fast cache marker and the inputs used to derive it."""
        return self._make_force_store().record_inputs(
            job_dir, input_digest, payload
        )

    def _thirdorder_pattern_paths(self):
        return sorted(
            self.fc3_root.glob("3RD.POSCAR.*"),
            key=self._displacement_pattern_sort_key,
        )

    def _fc4_pattern_paths(self):
        return sorted(
            self.fc4_root.glob("4TH.POSCAR.*"),
            key=self._displacement_pattern_sort_key,
        )

    def _export_thirdorder_fc3(self, generated):
        """Export a thirdorder FORCE_CONSTANTS_3RD in requested format(s)."""
        fc_format = self._fc_format()
        if fc_format in {"shengbte", "both"}:
            if generated.resolve() != self.shengbte_fc3_path.resolve():
                shutil.copyfile(generated, self.shengbte_fc3_path)
            print(f"  - Generated: {self.shengbte_fc3_path}")

        if fc_format in {"phono3py", "both"}:
            ph3 = self._make_phono3py()
            supercell = phonopy_to_ase(ph3.supercell)
            try:
                fcs = ForceConstants.read_shengBTE(
                    supercell,
                    str(generated),
                    self.prim,
                )
            except Exception as exc:
                raise RuntimeError(
                    "Could not convert thirdorder FORCE_CONSTANTS_3RD to "
                    "phono3py fc3.hdf5. Check that dim-fc3 is large enough "
                    "for cutoff-fc3 and matches the thirdorder sow/reap run."
                ) from exc
            self._write_phono3py_fc3(fcs.get_fc_array(3))

    def run_thirdorder_fc3_fitting(self):
        """Generate FC3 through the extracted thirdorder stage."""
        from nepkappa.stages.force_constants import ThirdOrderStage

        return ThirdOrderStage(self).run()

    def run_fc4_fitting(self):
        """Generate FC4 through the extracted Fourthorder stage."""
        from nepkappa.stages.force_constants import FourthOrderStage

        return FourthOrderStage(self).run()


    def compute_kappa(self):
        """Run the extracted phono3py thermal-transport stage."""
        from nepkappa.transport import Phono3pyTransportWorkflow

        transport = Phono3pyTransportWorkflow(
            self.cfg,
            self.output_dir,
            run_command=self._run_command,
            command_resolver=self._phono3py_command,
            needs_fc_flags=self._phono3py_needs_fc_flags,
        )
        return transport.run()

    def generate_force_constants(self, include_fc3=True):
        """Generate FC2/FC3 through the extracted strategy layer."""
        from nepkappa.force_constants import ForceConstantsWorkflow

        workflow = ForceConstantsWorkflow(
            self.cfg,
            prepare_structure=self.load_force_constant_structure,
            timed_stage=self._run_timed_stage,
            submit_slurm=self._submit_force_slurm,
            run_finite_displacement=self.run_finite_disp_fitting,
            run_hiphive=self.run_hiphive_fitting,
            run_thirdorder=self.run_thirdorder_fc3_fitting,
        )
        return workflow.run(include_fc3=include_fc3)

    def generate_fc4_force_constants(self):
        """Generate fourth-order force constants."""
        print("\n[Step 2] Generate Fourth-Order Force Constants")
        self.load_force_constant_structure()
        self._run_timed_stage("FourPhonon FC4", self.run_fc4_fitting)

    def calculate_kappa(self):
        """Compute thermal conductivity from existing force constants."""
        self._run_timed_stage("Phono3py kappa", self.compute_kappa)

    def run_relax(self):
        """Execute only the structure relaxation stage."""
        self._run_timed_stage("Relax structure", self.relax_structure_stage)
        self._print_timing_summary()

    def run_force_constants(self, include_fc3=True):
        """Execute only the force-constant generation stage."""
        result = self.generate_force_constants(include_fc3=include_fc3)
        self._print_timing_summary()
        return result

    def run_fc4(self):
        """Execute only fourth-order force-constant generation."""
        self.generate_fc4_force_constants()
        self._print_timing_summary()

    def run_kappa(self):
        """Execute only the thermal conductivity stage."""
        self.calculate_kappa()
        self._print_timing_summary()

    def run(self):
        """Execute relaxation, force-constant generation, and kappa."""
        self._run_timed_stage("Relax structure", self.relax_structure_stage)
        force_result = self.generate_force_constants()
        if force_result is not None:
            print(
                "  - Kappa is deferred until the Slurm force array and FC fitting "
                "job complete."
            )
            self._print_timing_summary()
            return force_result
        self.calculate_kappa()
        self._print_timing_summary()
