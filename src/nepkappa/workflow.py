# -*- coding: utf-8 -*-
# ==============================================================================
# Author: Shixian Liu, Fei Yin
# Date: 2025
# Description: Core workflow logic for Thermal Conductivity calculations 
#              using NEP (Neuroevolution Potential), HiPhive (optional), 
#              and Phono3py.
# ==============================================================================

import subprocess
import shlex
import shutil
import sys
import threading
import time
from pathlib import Path
import numpy as np

# ASE & Calorine
from ase.io import read, write
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from calorine.calculators import CPUNEP
from calorine.tools import relax_structure

# HiPhive
from hiphive import ClusterSpace, StructureContainer, ForceConstantPotential, ForceConstants
from hiphive.structure_generation import generate_mc_rattled_structures
from hiphive.utilities import prepare_structures
from hiphive import enforce_rotational_sum_rules
from trainstation import Optimizer

# Phonopy / Phono3py
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
import phono3py as phono3py_module
from phono3py import Phono3py
from phono3py.file_IO import write_fc2_to_hdf5, write_fc3_to_hdf5

try:
    from tqdm.auto import tqdm
except ImportError:
    tqdm = None

# --- Helper Functions ---
def format_duration(seconds):
    """Format elapsed seconds as a compact human-readable duration."""
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    secs = seconds % 60
    if hours:
        return f"{hours}h {minutes}m {secs:.2f}s"
    if minutes:
        return f"{minutes}m {secs:.2f}s"
    return f"{secs:.2f}s"

def format_activity_elapsed(seconds):
    """Format elapsed seconds like tqdm's compact clock."""
    total = int(seconds)
    hours, remainder = divmod(total, 3600)
    minutes, secs = divmod(remainder, 60)
    if hours:
        return f"{hours}:{minutes:02d}:{secs:02d}"
    return f"{minutes:02d}:{secs:02d}"

def progress_iter(iterable, enabled=True, **kwargs):
    """Return a tqdm iterator when available and enabled."""
    if enabled and tqdm is not None:
        return tqdm(iterable, **kwargs)
    return iterable

def run_activity_task(label, func, enabled=True, update_interval=1.0, log_interval=30.0):
    """Run a blocking task while showing elapsed activity feedback."""
    if not enabled:
        return func()

    result = {}

    def worker():
        try:
            result["value"] = func()
        except BaseException as exc:
            result["error"] = exc

    thread = threading.Thread(target=worker)
    thread.start()
    start = time.time()
    progress_stream = getattr(sys.stderr, "terminal_stream", sys.stderr)

    if tqdm is not None:
        with tqdm(
            total=None,
            desc=label,
            unit="s",
            bar_format="{desc}: {elapsed} elapsed",
            file=progress_stream,
            leave=False,
        ) as bar:
            while thread.is_alive():
                thread.join(update_interval)
                bar.update(update_interval)
    else:
        last_log = start
        while thread.is_alive():
            thread.join(update_interval)
            now = time.time()
            if now - last_log >= log_interval:
                print(
                    f"\r{label}: {format_activity_elapsed(now - start)} elapsed",
                    end="",
                    file=progress_stream,
                    flush=True,
                )
                last_log = now

    thread.join()
    if "error" in result:
        raise result["error"]
    print(f"{label}: {format_activity_elapsed(time.time() - start)} elapsed")
    return result.get("value")

def ase_to_phonopy(ase_atoms):
    """Convert ASE Atoms object to PhonopyAtoms object."""
    return PhonopyAtoms(
        symbols=ase_atoms.get_chemical_symbols(),
        scaled_positions=ase_atoms.get_scaled_positions(),
        cell=ase_atoms.cell,
    )

def phonopy_to_ase(ph_atoms):
    """Convert PhonopyAtoms object to ASE Atoms object."""
    return Atoms(
        symbols=ph_atoms.symbols,
        positions=ph_atoms.positions,
        cell=ph_atoms.cell,
        pbc=True
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
        self.disp_path = self.output_dir / "phono3py_disp.yaml"
        self.relaxed_poscar_path = self.output_dir / "POSCAR_relaxed"
        self.hiphive_model_path = self.output_dir / "hiphive_model.fcp"
        self.vasp_root = self.output_dir / getattr(config, "vasp_workdir", "vasp-runs")
        self.vasp_relax_root = self.output_dir / getattr(
            config, "vasp_relax_workdir", "vasp-relax"
        )
        self._nep_calc = None
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

    def _phono3py_command(self):
        """Return a robust phono3py executable command."""
        configured = getattr(self.cfg, "phono3py_command", None)
        if configured:
            return shlex.split(str(configured))

        executable = shutil.which("phono3py")
        if executable:
            return [executable]

        sibling_script = Path(sys.executable).with_name("phono3py")
        if sibling_script.exists():
            return [str(sibling_script)]

        return ["phono3py"]

    def _phono3py_needs_fc_flags(self):
        """Return True when the phono3py CLI needs explicit --fc2/--fc3 flags."""
        try:
            major = int(str(phono3py_module.__version__).split(".", 1)[0])
        except (AttributeError, TypeError, ValueError):
            return True
        return major < 4

    def _expected_kappa_name(self):
        mx, my, mz = self.cfg.mesh
        return f"kappa-m{mx}{my}{mz}.hdf5"

    def _make_phono3py(self):
        """Create a Phono3py object with the workflow's structure settings."""
        return Phono3py(
            ase_to_phonopy(self.prim),
            supercell_matrix=self.cfg.dim,
            phonon_supercell_matrix=self.cfg.dim,
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

    def _vasp_command(self):
        """Resolve the VASP command from YAML or common server locations."""
        if getattr(self.cfg, "vasp_command", None):
            return shlex.split(self.cfg.vasp_command)
        if getattr(self.cfg, "vasp_path", None):
            return [self.cfg.vasp_path]
        detected = self._detect_vasp_path()
        if detected is not None:
            print(f"  - Detected VASP executable: {detected}")
            return [str(detected)]
        raise FileNotFoundError(
            "VASP executable not found. Set calculator.vasp_path or "
            "calculator.vasp_command in the YAML file."
        )

    def _detect_vasp_path(self):
        """Find a likely VASP executable on common local/server paths."""
        candidates = [
            Path("/root/software/vasp.6.4.3/bin/vasp_std"),
            Path("/root/software/vasp.6.4.3/build/std/vasp"),
        ]
        for candidate in candidates:
            if candidate.is_file():
                return candidate
        for candidate in sorted(Path("/root/software").glob("vasp*/bin/vasp_std")):
            if candidate.is_file():
                return candidate
        resolved = shutil.which("vasp_std")
        if resolved:
            return Path(resolved)
        return None

    def _resolve_potcar_path(self, atoms):
        """Resolve a POTCAR file or assemble one from a potential directory."""
        configured = getattr(self.cfg, "potcar_path", None)
        symbols = []
        for symbol in atoms.get_chemical_symbols():
            if symbol not in symbols:
                symbols.append(symbol)

        if configured:
            path = Path(configured)
            if path.is_file():
                return path, (path,)
            if path.is_dir():
                return self._assemble_potcar_from_library(path, symbols)
            raise FileNotFoundError(f"POTCAR path not found: {path}")

        detected = self._detect_potcar_path(symbols)
        if detected is not None:
            print(f"  - Detected POTCAR source: {detected}")
            return detected
        raise FileNotFoundError(
            "POTCAR not found. Set calculator.potcar_path in the YAML file."
        )

    def _detect_potcar_path(self, symbols):
        """Find a likely POTCAR source under /root/software."""
        if len(symbols) == 1:
            symbol = symbols[0]
            candidates = [
                Path(f"/root/software/potpaw_PBE.64/{symbol}/POTCAR"),
                Path(f"/root/software/vasp.6.4.3/testsuite/POTCARS/POTCAR.{symbol}"),
                Path(f"/root/software/vasp.6.4.3/potpaw_PBE/{symbol}/POTCAR"),
                Path(f"/root/software/potpaw_PBE/{symbol}/POTCAR"),
            ]
            for candidate in candidates:
                if candidate.is_file():
                    return candidate, (candidate,)
        for library in (
            Path("/root/software/potpaw_PBE.64"),
            Path("/root/software/vasp.6.4.3/potpaw_PBE"),
            Path("/root/software/potpaw_PBE"),
            Path("/root/software/vasp.6.4.3/testsuite/POTCARS"),
        ):
            if library.is_dir():
                try:
                    return self._assemble_potcar_from_library(library, symbols)
                except FileNotFoundError:
                    pass
        return None

    def _assemble_potcar_from_library(self, library, symbols):
        """Return a single POTCAR file or create a combined POTCAR in result_dir."""
        chunks = []
        chosen = []
        for symbol in symbols:
            potcar = self._find_potcar_for_symbol(library, symbol)
            if potcar is None:
                raise FileNotFoundError(
                    f"missing POTCAR for {symbol} under {library}"
                )
            chunks.append(potcar.read_bytes())
            chosen.append(str(potcar))

        if len(chunks) == 1:
            path = Path(chosen[0])
            return path, (path,)

        combined = self.output_dir / "POTCAR.combined"
        combined.write_bytes(b"\n".join(chunks))
        spec = self.output_dir / "POTCAR.spec"
        spec.write_text("\n".join(chosen) + "\n", encoding="utf-8")
        return combined, tuple(Path(path) for path in chosen)

    def _find_potcar_for_symbol(self, library, symbol):
        candidates = [
            library / symbol / "POTCAR",
            library / f"{symbol}" / "POTCAR",
            library / f"{symbol}_sv" / "POTCAR",
            library / f"{symbol}_pv" / "POTCAR",
            library / f"POTCAR.{symbol}",
            library / "POTCAR",
        ]
        for candidate in candidates:
            if candidate.is_file():
                return candidate
        return None

    def _write_vasp_run_inputs(self, run_dir, atoms, incar_overrides=None, system=None):
        """Write POSCAR, INCAR, KPOINTS, POTCAR, and POTCAR.spec for one VASP run."""
        run_dir.mkdir(parents=True, exist_ok=True)
        write(str(run_dir / "POSCAR"), atoms, format="vasp", direct=True, vasp5=True)
        vasp_kwargs = self._combined_vasp_kwargs(incar_overrides, system=system)
        self._write_vasp_incar(run_dir / "INCAR", vasp_kwargs)
        self._write_vasp_kpoints(run_dir / "KPOINTS", vasp_kwargs)
        potcar_path, potcar_sources = self._resolve_potcar_path(atoms)
        shutil.copyfile(potcar_path, run_dir / "POTCAR")
        (run_dir / "POTCAR.spec").write_text(
            "\n".join(str(source) for source in potcar_sources) + "\n",
            encoding="utf-8",
        )

    def _run_vasp_forces(self, atoms, label, index):
        """Run one VASP single-point force calculation and return forces."""
        run_dir = self.vasp_root / label.lower() / f"{index:05d}"
        self._write_vasp_run_inputs(run_dir, atoms, system="NEP-kappa VASP single point")

        cmd = self._vasp_command()
        print(f"    VASP {label.upper()} #{index}: {run_dir}")
        ret = self._run_command(cmd, cwd=run_dir)
        if ret != 0:
            raise RuntimeError(f"VASP failed in {run_dir} with return code {ret}")
        return self._read_vasp_forces(run_dir)

    def _combined_vasp_kwargs(self, overrides=None, system=None):
        """Merge global VASP kwargs, workflow defaults, and stage-specific overrides."""
        kwargs = dict(getattr(self.cfg, "vasp_kwargs", {}) or {})
        defaults = {
            "system": system or "PESMaker single point",
            "gga": "PE",
            "lreal": "Auto",
            "ibrion": -1,
            "nsw": 0,
            "algo": "Normal",
            "ediff": 1.0e-6,
            "sigma": 0.02,
            "ismear": 0,
            "prec": "Accurate",
            "nelm": 150,
            "lwave": False,
            "lcharg": False,
        }
        for key, value in defaults.items():
            kwargs.setdefault(key, value)
        if overrides:
            kwargs.update(overrides)
        return kwargs

    def _write_vasp_incar(self, path, kwargs=None):
        if kwargs is None:
            kwargs = self._combined_vasp_kwargs()
        skip = {"kpts", "gamma", "txt", "directory", "command", "xc"}
        with path.open("w", encoding="utf-8") as handle:
            for key, value in kwargs.items():
                if key.lower() in skip:
                    continue
                handle.write(f"{key.upper()} = {self._format_vasp_value(value)}\n")

    def _write_vasp_kpoints(self, path, kwargs=None):
        if kwargs is None:
            kwargs = self._combined_vasp_kwargs()
        if "kspacing" in {str(key).lower() for key in kwargs}:
            if path.exists():
                path.unlink()
            return
        kpts = kwargs.get("kpts", [1, 1, 1])
        if isinstance(kpts, int):
            kpts = [kpts, kpts, kpts]
        if len(kpts) != 3:
            raise ValueError("vasp_kwargs.kpts must be an integer or a length-3 list")
        mode = "Gamma" if kwargs.get("gamma", True) else "Monkhorst-Pack"
        path.write_text(
            "Automatic mesh\n"
            "0\n"
            f"{mode}\n"
            f"{int(kpts[0])} {int(kpts[1])} {int(kpts[2])}\n"
            "0 0 0\n",
            encoding="utf-8",
        )

    def _format_vasp_value(self, value):
        if isinstance(value, bool):
            return ".TRUE." if value else ".FALSE."
        if isinstance(value, (list, tuple)):
            return " ".join(self._format_vasp_value(item) for item in value)
        return str(value)

    def _read_vasp_forces(self, run_dir):
        for filename in ("vasprun.xml", "OUTCAR"):
            path = run_dir / filename
            if not path.exists():
                continue
            try:
                atoms = read(str(path), index=-1)
                return atoms.get_forces()
            except Exception as exc:
                last_error = exc
        raise RuntimeError(
            f"Could not read VASP forces from {run_dir / 'vasprun.xml'} "
            f"or {run_dir / 'OUTCAR'}"
        ) from locals().get("last_error")

    def _default_vasp_relax_stages(self):
        return {
            "coarse": {
                "system": "NEP-kappa VASP coarse relaxation",
                "nsw": 80,
                "ibrion": 2,
                "isif": 3,
                "ediff": 1.0e-5,
                "ediffg": -0.05,
                "prec": "Normal",
            },
            "fine": {
                "system": "NEP-kappa VASP fine relaxation",
                "nsw": 150,
                "ibrion": 2,
                "isif": 3,
                "ediff": 1.0e-6,
                "ediffg": -0.01,
                "prec": "Accurate",
            },
        }

    def _vasp_relax_stages(self):
        stages = getattr(self.cfg, "vasp_relax_stages", None)
        if not stages:
            stages = self._default_vasp_relax_stages()
        if not isinstance(stages, dict):
            raise ValueError("vasp_relax.stages must be a mapping")
        ordered = []
        for name in ("coarse", "fine"):
            if name in stages:
                value = stages[name] or {}
                if not isinstance(value, dict):
                    raise ValueError(f"vasp_relax.stages.{name} must be a mapping")
                ordered.append((name, value))
        for name, value in stages.items():
            if name in {"coarse", "fine"}:
                continue
            if not isinstance(value, dict):
                raise ValueError(f"vasp_relax.stages.{name} must be a mapping")
            ordered.append((name, value))
        return ordered

    def _read_vasp_relaxed_structure(self, run_dir):
        for filename in ("CONTCAR", "vasprun.xml", "OUTCAR"):
            path = run_dir / filename
            if not path.exists() or path.stat().st_size == 0:
                continue
            try:
                return read(str(path), index=-1)
            except Exception as exc:
                last_error = exc
        raise RuntimeError(
            f"Could not read relaxed VASP structure from {run_dir}"
        ) from locals().get("last_error")

    def _run_vasp_relax_stage(self, atoms, stage_name, stage_index, overrides):
        run_dir = self.vasp_relax_root / f"{stage_index:02d}-{stage_name}"
        system = overrides.get(
            "system", f"NEP-kappa VASP {stage_name} relaxation"
        )
        self._write_vasp_run_inputs(
            run_dir,
            atoms,
            incar_overrides=overrides,
            system=system,
        )
        cmd = self._vasp_command()
        print(f"    VASP relax {stage_name}: {run_dir}")
        ret = self._run_command(cmd, cwd=run_dir)
        if ret != 0:
            raise RuntimeError(f"VASP relaxation failed in {run_dir} with return code {ret}")
        relaxed = self._read_vasp_relaxed_structure(run_dir)
        write(str(run_dir / "POSCAR_next"), relaxed, format="vasp", direct=True, vasp5=True)
        return relaxed

    def _run_vasp_relaxation(self, atoms):
        print("  - Relaxing structure using VASP")
        current = atoms.copy()
        for index, (stage_name, overrides) in enumerate(self._vasp_relax_stages(), 1):
            print(f"  - VASP relaxation stage {index}: {stage_name}")
            current = self._run_vasp_relax_stage(current, stage_name, index, overrides)
        write(str(self.relaxed_poscar_path), current, format="vasp", direct=True, vasp5=True)
        print(f"  - Relaxed structure saved to {self.relaxed_poscar_path}")
        return current

    def _calculate_forces(self, atoms, label, index):
        """Calculate forces with the configured backend."""
        calculator = self._calculator_name()
        if calculator == "nep":
            atoms.calc = self._make_nep_calculator()
            return atoms.get_forces()
        elif calculator == "vasp":
            return self._run_vasp_forces(atoms, label, index)
        else:
            raise ValueError(f"Unsupported calculator: {calculator}")

    def relax_structure_stage(self):
        """Step 1: Load and optionally relax the primitive cell."""
        print("\n[Step 1] Relax Structure")
        print(f"  - Reading from {self.cfg.poscar}")
        self.prim = read(self.cfg.poscar)

        if self.cfg.do_relax:
            if self._calculator_name() == "vasp":
                self.prim = self._run_vasp_relaxation(self.prim)
            else:
                print(f"  - Relaxing structure using NEP: {self.cfg.nep_model}")
                self.prim.calc = self._make_nep_calculator()
                relax_structure(self.prim, fmax=1e-3)
                relax_structure(self.prim, fmax=1e-5)
                write(str(self.relaxed_poscar_path), self.prim)
                print(f"  - Relaxed structure saved to {self.relaxed_poscar_path}")
        else:
            write(str(self.relaxed_poscar_path), self.prim, format="vasp", direct=True, vasp5=True)
            print("  - Relaxation disabled; copied input structure for downstream steps")
            print(f"  - Structure saved to {self.relaxed_poscar_path}")

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

    def run_hiphive_fitting(self):
        """Step 1 (Path A): Fit Force Constants using HiPhive (Compressive Sensing)."""
        print("\n[Step 2 - HiPhive] Generating Training Data & Fitting")
        cfg = self.cfg
        np.random.seed(42)
        self._save_phono3py_metadata()

        # 1. Create supercell
        nx, ny, nz = cfg.dim
        atoms_ideal = self.prim.repeat((nx, ny, nz))
        
        # 2. Generate rattled structures
        print(f"  - Generating {cfg.n_structures} rattled structures (std={cfg.rattle_std}, min_dist={cfg.min_dist})")   
        structures = generate_mc_rattled_structures(
            atoms_ideal, 
            cfg.n_structures, 
            cfg.rattle_std, 
            cfg.min_dist
        )

        # 3. Compute forces
        print(
            f"  - Computing forces with {self._calculator_name().upper()} "
            f"for {len(structures)} structures..."
        )
        for i, at in enumerate(progress_iter(
            structures,
            enabled=self.show_progress,
            total=len(structures),
            desc=self._force_desc("HiPhive"),
            unit="structure"
        )):
            f = self._calculate_forces(at, "hiphive", i + 1)
            at.calc = SinglePointCalculator(at, forces=f)
            at.get_forces()
            if (i+1) % 50 == 0:
                print(f"    Processed {i+1}/{len(structures)}")

        # 4. Train HiPhive Potential
        print(f"  - Fitting Force Constants (Cutoffs: {cfg.cutoffs})")
        cs = ClusterSpace(self.prim, cfg.cutoffs)
        sc = StructureContainer(cs)
        for s in prepare_structures(structures, atoms_ideal):
            sc.add_structure(s)
            
        opt = Optimizer(sc.get_fit_data())
        opt.train()
        print(f"    RMSE: {opt.rmse_train:.6f}")

        # 5. Enforce sum rules and save
        params = enforce_rotational_sum_rules(cs, opt.parameters, ['Huang', 'Born-Huang'])
        fcp = ForceConstantPotential(cs, params)
        fcp.write(str(self.hiphive_model_path))
        print(f"  - HiPhive model saved to {self.hiphive_model_path}")

        # 6. Export to Phono3py FC2 and FC3
        print("  - Exporting FC2 and FC3 from HiPhive model")
        phonopy_obj = Phonopy(ase_to_phonopy(self.prim), supercell_matrix=np.diag(cfg.dim))
        supercell_ase = phonopy_to_ase(phonopy_obj.supercell)
        fcs = fcp.get_force_constants(supercell_ase)
        
        fcs.write_to_phonopy(str(self.fc2_path))
        fcs.write_to_phono3py(str(self.fc3_path))
        print(f"  - Generated: {self.fc2_path}, {self.fc3_path}")

    def run_finite_disp_fitting(self):
        """Step 1 (Path B): Standard Finite Displacement Method using Phono3py."""
        print("\n[Step 2 - FiniteDisp] Phono3py Finite Displacement Method")
        cfg = self.cfg
        
        # 1. Initialize Phono3py
        ph3 = self._make_phono3py()
        
        # 2. Generate displacements
        ph3.generate_displacements()
        ph3.save(str(self.disp_path))
        
        fc3_scs = ph3.supercells_with_displacements
        fc2_scs = ph3.phonon_supercells_with_displacements
        print(f"  - Generated {len(fc3_scs)} FC3 supercells and {len(fc2_scs)} FC2 supercells")

        # 3. Compute forces for FC2
        print("  - Computing forces for FC2...")
        forces_fc2 = []
        for i, sc in enumerate(progress_iter(
            fc2_scs,
            enabled=self.show_progress,
            total=len(fc2_scs),
            desc=self._force_desc("FC2"),
            unit="structure"
        )):
            atoms = phonopy_to_ase(sc)
            forces_fc2.append(self._calculate_forces(atoms, "fc2", i + 1))
        ph3.phonon_forces = np.array(forces_fc2)

        # 4. Compute forces for FC3
        print("  - Computing forces for FC3...")
        forces_fc3 = []
        for i, sc in enumerate(progress_iter(
            fc3_scs,
            enabled=self.show_progress,
            total=len(fc3_scs),
            desc=self._force_desc("FC3"),
            unit="structure"
        )):
            atoms = phonopy_to_ase(sc)
            forces_fc3.append(self._calculate_forces(atoms, "fc3", i + 1))
        ph3.forces = np.array(forces_fc3)

        # 5. Produce and save FCs
        fc_calculator = getattr(cfg, "fc_calculator", "traditional")
        solver_name = None if fc_calculator == "traditional" else fc_calculator
        solver_options = getattr(cfg, "fc_calculator_options", None)
        symmetrize_traditional = fc_calculator == "traditional"
        produce_kwargs = {
            "fc_calculator": solver_name,
            "fc_calculator_options": solver_options,
        }
        solver_label = fc_calculator if fc_calculator else "traditional"
        print(f"  - Force constants solver: {solver_label}")
        print("  - Producing FC2 and FC3...")
        run_activity_task(
            "Producing FC2",
            lambda: ph3.produce_fc2(**produce_kwargs),
            enabled=self.show_progress,
        )
        if symmetrize_traditional:
            print("  - Symmetrizing FC2...")
            ph3.symmetrize_fc2()
        write_fc2_to_hdf5(ph3.fc2, filename=str(self.fc2_path))
        run_activity_task(
            "Producing FC3",
            lambda: ph3.produce_fc3(**produce_kwargs),
            enabled=self.show_progress,
        )
        if symmetrize_traditional:
            print("  - Symmetrizing FC3...")
            ph3.symmetrize_fc3()
        write_fc3_to_hdf5(ph3.fc3, filename=str(self.fc3_path))
        print(f"  - Generated: {self.fc2_path}, {self.fc3_path}")
        

    def compute_kappa(self):
        """Step 2: Run phono3py CLI to compute thermal conductivity."""
        print("\n[Step 3] Compute Kappa with Phono3py CLI")
        cfg = self.cfg
        
        mx, my, mz = cfg.mesh

        missing_fc = [
            path for path in (self.fc2_path, self.fc3_path, self.disp_path)
            if not path.exists()
        ]
        if missing_fc:
            missing = ", ".join(str(path) for path in missing_fc)
            raise FileNotFoundError(
                f"Missing required phono3py file(s): {missing}. "
                f"Run `nepkappa fc` first or place fc2.hdf5, fc3.hdf5, "
                f"and {self.disp_path.name} in {self.output_dir}."
            )
        
        if cfg.method == 'lbte':
            method_flags = ["--lbte"]
            print("  - Method: LBTE (Linearized Boltzmann Transport Equation)")
        else:
            method_flags = ["--br", "--nu"]
            print("  - Method: RTA (Relaxation Time Approximation)")

        cmd = [
            *self._phono3py_command(),
            self.disp_path.name,
        ]
        if self._phono3py_needs_fc_flags():
            cmd.extend(["--fc2", "--fc3"])
        cmd.extend([*method_flags, "--mesh", str(mx), str(my), str(mz)])

        if cfg.wigner:
            print("  - Wigner transport: enabled via phono3py-wte (--tt wte)")
            cmd.extend(["--tt", "wte"])

        if len(cfg.temps) == 3:
            tmin, tmax, tstep = cfg.temps
            cmd.extend(["--tmin", str(tmin), "--tmax", str(tmax), "--tstep", str(tstep)])
        else:
            cmd.extend(["--ts", str(cfg.temps[0])])
        
        print(f"  - Output directory: {self.output_dir}")
        print(f"  - Running command: {shlex.join(cmd)}")
        ret = self._run_command(cmd, cwd=self.output_dir)
        
        if ret == 0:
            print("\n[Done] Phono3py finished successfully.")
            print(f"Check {self.output_dir / self._expected_kappa_name()} for results.")
        else:
            print(f"\n[Error] Phono3py failed with return code {ret}")

    def generate_force_constants(self):
        """Generate FC2 and FC3 files."""
        print("\n[Step 2] Generate Force Constants")
        self.load_force_constant_structure()

        if self.cfg.use_hiphive:
            self._run_timed_stage("HiPhive fitting", self.run_hiphive_fitting)
        else:
            self._run_timed_stage("Finite displacement", self.run_finite_disp_fitting)

    def calculate_kappa(self):
        """Compute thermal conductivity from existing force constants."""
        self._run_timed_stage("Phono3py kappa", self.compute_kappa)

    def run_relax(self):
        """Execute only the structure relaxation stage."""
        self._run_timed_stage("Relax structure", self.relax_structure_stage)
        self._print_timing_summary()

    def run_force_constants(self):
        """Execute only the force-constant generation stage."""
        self.generate_force_constants()
        self._print_timing_summary()

    def run_kappa(self):
        """Execute only the thermal conductivity stage."""
        self.calculate_kappa()
        self._print_timing_summary()

    def run(self):
        """Execute relaxation, force-constant generation, and kappa."""
        self._run_timed_stage("Relax structure", self.relax_structure_stage)
        self.generate_force_constants()
        self.calculate_kappa()
        self._print_timing_summary()
