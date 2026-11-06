# -*- coding: utf-8 -*-
# ==============================================================================
# Author: Shixian Liu, Fei Yin
# Date: 2025
# Description: Shared scientific context, force caches and export helpers.
#              Command orchestration is owned by application.py.
# ==============================================================================

import subprocess
import hashlib
import shlex
import shutil
import re
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np

# ASE & Calorine
from ase.io import read
from calorine.calculators import CPUNEP

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
from nepkappa.coordinates import canonicalize_atoms
from nepkappa.output import print_detail
from nepkappa.calculators import calculator_backend_from_config
from nepkappa.provenance import (
    canonical_data,
    data_sha256,
    file_identity,
    installed_versions,
)
from nepkappa.runtime import (
    ase_to_phonopy,
    phonopy_to_ase,
)

# --- Main Logic Class ---
class NEPPhononWorkflow:
    """Shared structures, calculators and force artifacts for scientific stages.

    Execution order and stage timing belong to application.py. The context is
    also used directly by the Slurm force worker; it does not execute a plan.
    """
    def __init__(self, config):
        """
        Initialize the workflow.
        """
        self.cfg = config
        self.show_progress = getattr(config, "progress", True)
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
        self._nep_batch_identity = None
        self._external_backend = None
        self._vasp_backend = None
        self._force_store = None
        self._cache_software_versions = None
        self.prim = None 


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
                print(line, end="", flush=True)
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
        # HiPhive requires fractional coordinates away from the upper boundary.
        # Wrap copies only: atom order and FC indices must remain unchanged.
        supercell = supercell.copy()
        primitive = primitive.copy()
        supercell.wrap(eps=1e-6)
        primitive.wrap(eps=1e-6)
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
        self.prim = canonicalize_atoms(self.prim)
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


    def _calculate_forces(self, atoms, label, index):
        """Calculate forces with the configured backend."""
        calculator = self._calculator_name()
        if calculator == "nep":
            if getattr(self.cfg, "nep_force_backend", "cpu") == "gpu-batch":
                raise RuntimeError("NEP GPU forces must be evaluated through the batch force path.")
            atoms.calc = self._make_nep_calculator()
            return atoms.get_forces()
        elif calculator == "vasp":
            return self._make_vasp_backend().run_forces(atoms, label, index)
        else:
            return self._make_external_backend().calculate_forces(atoms)

    def _force_job_dir(self, label, index):
        """Return the stable work directory for one displaced structure."""
        return self._make_force_store().job_dir(
            label,
            index,
            calculator=self._calculator_name(),
        )

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

    def _calculate_nep_batch_cached(self, jobs):
        """Evaluate ordered (atoms, label, index) jobs, batching only cache misses."""
        from nepkappa.nep_batch import predict_forces, resolve_nep_executable

        command = resolve_nep_executable(self.cfg.nep_command)
        print(f"  - NEP GPU batches: at most {self.cfg.nep_batch_size} structures per call")
        return self._make_force_store().calculate_batch_or_reuse(
            jobs, calculator="nep", fingerprint=self._job_input_fingerprint,
            batch_size=self.cfg.nep_batch_size,
            calculate=lambda structures: predict_forces(
                structures, model=self.cfg.nep_model, command=command
            ),
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
            self.prim = canonicalize_atoms(read(str(self.relaxed_poscar_path)))
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
            pair_cutoff = getattr(self.cfg, "cutoff_fc3", None)
            if pair_cutoff is None:
                ph3.generate_displacements()
            else:
                ph3.generate_displacements(cutoff_pair_distance=pair_cutoff)
            ph3.save(str(self.disp_path))
            for index, supercell in enumerate(ph3.supercells_with_displacements, 1):
                if supercell is None:  # Keep original indices for cutoff-excluded pairs.
                    continue
                jobs.append(
                    self._stage_force_job(phonopy_to_ase(supercell), "fc3", index)
                )
        return jobs

    def _submit_force_slurm(self, include_fc3=True):
        """Prepare displaced structures and dispatch their force calculations."""
        from nepkappa.slurm import run_force_slurm

        jobs = self._prepare_force_slurm_jobs(include_fc3=include_fc3)
        return run_force_slurm(self.cfg, self.output_dir, jobs)


    def load_force_constant_structure(self):
        """Load the structure used by the force-constant stage."""
        if self.prim is not None:
            self.prim = canonicalize_atoms(self.prim)
            print_detail(f"  - Using in-memory structure from the relax stage")
            return

        if self.cfg.do_relax:
            if not self.relaxed_poscar_path.exists():
                raise FileNotFoundError(
                    f"Relaxed structure not found: {self.relaxed_poscar_path}. "
                    "Run `nepkappa relax` first, or use `nepkappa run` for the full workflow."
                )
            print_detail(f"  - Reading relaxed structure from {self.relaxed_poscar_path}")
            self.prim = canonicalize_atoms(read(str(self.relaxed_poscar_path)))
        else:
            print_detail(f"  - Reading input structure from {self.cfg.poscar}")
            self.prim = canonicalize_atoms(read(self.cfg.poscar))


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
                sources = [self._make_vasp_backend().find_potcar(path, symbol) for symbol in symbols]
        else:
            detected = self._make_vasp_backend().detect_potcar(symbols)
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
            if getattr(self.cfg, "nep_force_backend", "cpu") == "gpu-batch":
                from nepkappa import nep_batch

                if self._nep_batch_identity is None:
                    executable = nep_batch.resolve_nep_executable(self.cfg.nep_command)
                    self._nep_batch_identity = {
                        "executable": file_identity(executable),
                        "adapter": file_identity(nep_batch.__file__),
                    }
                payload["calculator_inputs"].update({
                    "force_backend": "gpu-batch",
                    "batch_size": self.cfg.nep_batch_size,
                    **self._nep_batch_identity,
                })
        elif calculator == "vasp":
            command = getattr(self.cfg, "vasp_command", None)
            if command:
                command = shlex.split(str(command))
            elif getattr(self.cfg, "vasp_path", None):
                command = [str(self.cfg.vasp_path)]
            else:
                detected = self._make_vasp_backend().detect_path()
                command = [str(detected)] if detected else ["vasp_std"]
            atoms = read(str(pattern), format="vasp")
            payload["calculator_inputs"] = {
                "command": self._command_cache_identity(command),
                "incar": canonical_data(self._make_vasp_backend().combined_kwargs()),
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
