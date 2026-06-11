# workflow.py
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

def progress_iter(iterable, enabled=True, **kwargs):
    """Return a tqdm iterator when available and enabled."""
    if enabled and tqdm is not None:
        return tqdm(iterable, **kwargs)
    return iterable

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
        self.output_name = getattr(config, "output_name", None)
        self.fc2_path = self.output_dir / "fc2.hdf5"
        self.fc3_path = self.output_dir / "fc3.hdf5"
        self.disp_path = self.output_dir / "phono3py_disp.yaml"
        self.relaxed_poscar_path = self.output_dir / "POSCAR_relaxed"
        self.hiphive_model_path = self.output_dir / "hiphive_model.fcp"
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

    def _expected_kappa_name(self):
        mx, my, mz = self.cfg.mesh
        return f"kappa-m{mx}{my}{mz}.hdf5"

    def _final_kappa_path(self):
        if not self.output_name:
            return self.output_dir / self._expected_kappa_name()
        name = Path(self.output_name).name
        if not name.endswith(".hdf5"):
            name = f"{name}.hdf5"
        return self.output_dir / name

    def _rename_kappa_output(self):
        if not self.output_name:
            return

        default_path = self.output_dir / self._expected_kappa_name()
        final_path = self._final_kappa_path()
        if final_path == default_path:
            return
        if default_path.exists():
            if final_path.exists():
                final_path.unlink()
            default_path.replace(final_path)
            print(f"  - Renamed kappa HDF5 to {final_path}")
            return

        candidates = sorted(
            self.output_dir.glob("kappa-*.hdf5"),
            key=lambda path: path.stat().st_mtime,
            reverse=True,
        )
        if candidates:
            if final_path.exists():
                final_path.unlink()
            candidates[0].replace(final_path)
            print(f"  - Renamed kappa HDF5 to {final_path}")
        else:
            print("  - Warning: no kappa HDF5 file was found to rename.")

    def prepare_structure(self):
        """Step 0: Load and optionally relax the primitive cell."""
        print("\n[Step 0] Prepare Structure")
        print(f"  - Reading from {self.cfg.poscar}")
        self.prim = read(self.cfg.poscar)

        if self.cfg.do_relax:
            print(f"  - Relaxing structure using NEP: {self.cfg.nep_model}")
            self.prim.calc = CPUNEP(self.cfg.nep_model)
            relax_structure(self.prim, fmax=1e-3)
            relax_structure(self.prim, fmax=1e-5)
            write(str(self.relaxed_poscar_path), self.prim)
            print(f"  - Relaxed structure saved to {self.relaxed_poscar_path}")
        else:
            print("  - Skipping relaxation (using input structure as is)")

    def run_hiphive_fitting(self):
        """Step 1 (Path A): Fit Force Constants using HiPhive (Compressive Sensing)."""
        print("\n[Step 1 - HiPhive] Generating Training Data & Fitting")
        cfg = self.cfg
        np.random.seed(42)

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

        # 3. Compute forces with NEP
        print(f"  - Computing forces with NEP for {len(structures)} structures...")
        calc = CPUNEP(cfg.nep_model)
        for i, at in enumerate(progress_iter(
            structures,
            enabled=self.show_progress,
            total=len(structures),
            desc="NEP forces (HiPhive)",
            unit="structure"
        )):
            at.calc = calc
            f = at.get_forces()
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
        print("\n[Step 1 - FiniteDisp] Phono3py Finite Displacement Method")
        cfg = self.cfg
        
        # 1. Initialize Phono3py
        ph3 = Phono3py(
            ase_to_phonopy(self.prim),
            supercell_matrix=cfg.dim,
            phonon_supercell_matrix=cfg.dim, # Assuming same size for FC2 and FC3
            primitive_matrix="auto"
        )
        
        # 2. Generate displacements
        ph3.generate_displacements()
        ph3.save(str(self.disp_path))
        
        fc3_scs = ph3.supercells_with_displacements
        fc2_scs = ph3.phonon_supercells_with_displacements
        print(f"  - Generated {len(fc3_scs)} FC3 supercells and {len(fc2_scs)} FC2 supercells")

        # 3. Compute forces for FC2
        print("  - Computing forces for FC2...")
        calc = CPUNEP(cfg.nep_model)
        forces_fc2 = []
        for sc in progress_iter(
            fc2_scs,
            enabled=self.show_progress,
            total=len(fc2_scs),
            desc="NEP forces (FC2)",
            unit="structure"
        ):
            atoms = phonopy_to_ase(sc)
            atoms.calc = calc
            forces_fc2.append(atoms.get_forces())
        ph3.phonon_forces = np.array(forces_fc2)

        # 4. Compute forces for FC3
        print("  - Computing forces for FC3...")
        forces_fc3 = []
        for sc in progress_iter(
            fc3_scs,
            enabled=self.show_progress,
            total=len(fc3_scs),
            desc="NEP forces (FC3)",
            unit="structure"
        ):
            atoms = phonopy_to_ase(sc)
            atoms.calc = calc
            forces_fc3.append(atoms.get_forces())
        ph3.forces = np.array(forces_fc3)

        # 5. Produce and save FCs
        print("  - Producing FC2 and FC3...")
        ph3.produce_fc2()
        write_fc2_to_hdf5(ph3.fc2, filename=str(self.fc2_path))
        ph3.produce_fc3()
        write_fc3_to_hdf5(ph3.fc3, filename=str(self.fc3_path))
        print(f"  - Generated: {self.fc2_path}, {self.fc3_path}")
        

    def compute_kappa(self):
        """Step 2: Run phono3py CLI to compute thermal conductivity."""
        print("\n[Step 2] Compute Kappa with Phono3py CLI")
        cfg = self.cfg
        
        if cfg.do_relax:
            poscar_name = str(self.relaxed_poscar_path)
        else:
            poscar_name = str(Path(cfg.poscar).resolve())
        
        nx, ny, nz = cfg.dim
        mx, my, mz = cfg.mesh

        missing_fc = [path for path in (self.fc2_path, self.fc3_path) if not path.exists()]
        if missing_fc:
            missing = ", ".join(str(path) for path in missing_fc)
            raise FileNotFoundError(
                f"Missing force constant file(s): {missing}. "
                f"Run with --fc2fc3 true or place fc2.hdf5/fc3.hdf5 in {self.output_dir}."
            )
        
        if cfg.method == 'lbte':
            method_flags = ["--lbte"]
            print("  - Method: LBTE (Linearized Boltzmann Transport Equation)")
        else:
            method_flags = ["--br", "--nu"]
            print("  - Method: RTA (Relaxation Time Approximation)")

        cmd = [
            "phono3py",
            "-c", poscar_name,
            "--dim", str(nx), str(ny), str(nz),
            "--fc2", "--fc3",
            *method_flags,
            "--mesh", str(mx), str(my), str(mz),
        ]

        if cfg.wigner:
            cmd.append("--wigner")

        if len(cfg.temps) == 3:
            tmin, tmax, tstep = cfg.temps
            cmd.extend(["--tmin", str(tmin), "--tmax", str(tmax), "--tstep", str(tstep)])
        else:
            cmd.extend(["--ts", str(cfg.temps[0])])
        
        print(f"  - Output directory: {self.output_dir}")
        print(f"  - Running command: {shlex.join(cmd)}")
        ret = self._run_command(cmd, cwd=self.output_dir)
        
        if ret == 0:
            self._rename_kappa_output()
            print("\n[Done] Phono3py finished successfully.")
            print(f"Check {self._final_kappa_path()} for results.")
        else:
            print(f"\n[Error] Phono3py failed with return code {ret}")

    def run(self):
        """Execute the full workflow."""
        
        if self.cfg.fc2fc3:
            self._run_timed_stage("Prepare structure", self.prepare_structure)
            
            if self.cfg.use_hiphive:
                self._run_timed_stage("HiPhive fitting", self.run_hiphive_fitting)
            else:
                self._run_timed_stage("Finite displacement", self.run_finite_disp_fitting)
            
        self._run_timed_stage("Phono3py kappa", self.compute_kappa)
        self._print_timing_summary()
