# nep_workflow.py
# -*- coding: utf-8 -*-
# ==============================================================================
# Author: Shixian Liu, Fei Yin
# Date: 2025
# Description: Core workflow logic for Thermal Conductivity calculations 
#              using NEP (Neuroevolution Potential), HiPhive (optional), 
#              and Phono3py.
# ==============================================================================

import os
import subprocess
import numpy as np
import h5py

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

# --- Helper Functions ---
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
        self.prim = None 

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
            write("POSCAR_relaxed", self.prim)
            print("  - Relaxed structure saved to POSCAR_relaxed")
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
        for i, at in enumerate(structures):
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
        fcp.write("hiphive_model.fcp")
        print("  - HiPhive model saved to 'hiphive_model.fcp'")

        # 6. Export to Phono3py/ShengBTE formats
        print("  - Exporting FC2 and FC3 from HiPhive model")
        phonopy_obj = Phonopy(ase_to_phonopy(self.prim), supercell_matrix=np.diag(cfg.dim))
        supercell_ase = phonopy_to_ase(phonopy_obj.supercell)
        fcs = fcp.get_force_constants(supercell_ase)
        
        fcs.write_to_phonopy("fc2.hdf5")
        fcs.write_to_phono3py("fc3.hdf5")
        fcs.write_to_shengBTE("FORCE_CONSTANTS_3RD", self.prim)
        fcs.write_to_phonopy("FORCE_CONSTANTS", format="text")
        print("  - Generated: fc2.hdf5, fc3.hdf5, FORCE_CONSTANTS_3RD")

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
        ph3.save("phono3py_disp.yaml")
        
        fc3_scs = ph3.supercells_with_displacements
        fc2_scs = ph3.phonon_supercells_with_displacements
        print(f"  - Generated {len(fc3_scs)} FC3 supercells and {len(fc2_scs)} FC2 supercells")

        # 3. Compute forces for FC2
        print("  - Computing forces for FC2...")
        calc = CPUNEP(cfg.nep_model)
        forces_fc2 = []
        for sc in fc2_scs:
            atoms = phonopy_to_ase(sc)
            atoms.calc = calc
            forces_fc2.append(atoms.get_forces())
        ph3.phonon_forces = np.array(forces_fc2)

        # 4. Compute forces for FC3
        print("  - Computing forces for FC3...")
        forces_fc3 = []
        for sc in fc3_scs:
            atoms = phonopy_to_ase(sc)
            atoms.calc = calc
            forces_fc3.append(atoms.get_forces())
        ph3.forces = np.array(forces_fc3)

        # 5. Produce and save FCs
        print("  - Producing FC2 and FC3...")
        ph3.produce_fc2()
        write_fc2_to_hdf5(ph3.fc2)
        ph3.produce_fc3()
        write_fc3_to_hdf5(ph3.fc3)
        
        # 6. Convert formats
        print("  - Converting to ShengBTE format via HiPhive interface...")
        supercell_ase = phonopy_to_ase(ph3.supercell)
        fcs = ForceConstants.from_arrays(
            supercell=supercell_ase,
            fc2_array=ph3.fc2,
            fc3_array=ph3.fc3
        )
        fcs.write_to_shengBTE("FORCE_CONSTANTS_3RD", self.prim)
        fcs.write_to_phonopy("FORCE_CONSTANTS", format="text")
        print("  - Generated: fc2.hdf5, fc3.hdf5, FORCE_CONSTANTS_3RD, FORCE_CONSTANTS")

    def compute_kappa(self):
        """Step 2: Run phono3py CLI to compute thermal conductivity."""
        print("\n[Step 2] Compute Kappa with Phono3py CLI")
        cfg = self.cfg
        
        if cfg.do_relax:
            poscar_name = "POSCAR_relaxed"
        else:
            poscar_name = cfg.poscar
        
        nx, ny, nz = cfg.dim
        mx, my, mz = cfg.mesh
        tmin, tmax, tstep = cfg.temps
        
        if cfg.method == 'lbte':
            method_flag = "--lbte" 
            print("  - Method: LBTE (Linearized Boltzmann Transport Equation)")
        else:
            method_flag = "--br"
            print("  - Method: RTA (Relaxation Time Approximation)")

        cmd = (
            f"phono3py -c {poscar_name} "
            f"--dim {nx} {ny} {nz} "
            f"--fc2 --fc3 "
            f"{method_flag} " 
            f"--mesh {mx} {my} {mz} "
            f"--tmin {tmin} "
            f"--tmax {tmax} "
            f"--tstep {tstep} "
        )
        
        print(f"  - Running command: {cmd}")
        ret = subprocess.call(cmd, shell=True)
        
        if ret == 0:
            print("\n[Done] Phono3py finished successfully.")
            if cfg.method == 'lbte':
                print(f"Check kappa-m{mx}{my}{mz}.hdf5 for results.")
            else:
                print(f"Check kappa-m{mx}{my}{mz}.hdf5 (RTA) for results.")
        else:
            print(f"\n[Error] Phono3py failed with return code {ret}")

    def run(self):
        """Execute the full workflow."""
        self.prepare_structure()
        
        if self.cfg.use_hiphive:
            self.run_hiphive_fitting()
        else:
            self.run_finite_disp_fitting()
            
        self.compute_kappa()