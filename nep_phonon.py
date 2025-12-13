#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess

import numpy as np
import h5py

from ase.io import read, write
from ase import Atoms

from calorine.calculators import CPUNEP
from calorine.tools import relax_structure

from hiphive.structure_generation import generate_mc_rattled_structures
from hiphive.utilities import prepare_structures
from hiphive import ClusterSpace, StructureContainer, ForceConstantPotential
from hiphive import enforce_rotational_sum_rules

from trainstation import Optimizer

from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

class PhononWorkflow:
    """
    使用 NEP + HiPhive + phono3py 的声子/热导率工作流。

    参数从 argparse.Namespace 传入：
      - poscar, nep_model
      - dim (nx, ny, nz)  同时用于训练超胞和 phonopy 超胞
      - n_structures, rattle_std, min_dist, seed, do_relax
      - cutoffs
      - mesh_x, mesh_y, mesh_z, temps
    """
    def __init__(self, args):
        self.cfg = args   # argparse.Namespace

    # ---------- Step 1: 生成训练结构 ----------
    def build_training_set(self):
        cfg = self.cfg
        np.random.seed(cfg.seed)

        nx, ny, nz = cfg.dim  # nz 这里一般是 1（薄膜）

        print("[Workflow] Step 1: build training set")
        print(f"[Workflow]  - Read primitive cell from {cfg.poscar}")
        prim = read(cfg.poscar)

        # 可选 NEP 松弛
        if cfg.do_relax:
            print(f"[Workflow]  - Relax primitive cell with NEP: {cfg.nep_model}")
            prim.calc = CPUNEP(cfg.nep_model)
            relax_structure(prim, fmax=1e-3)
            relax_structure(prim, fmax=1e-5)

        # 面内扩胞（这里 dim 的 nx, ny 既用作训练超胞，也用作 phonopy 超胞）
        print(f"[Workflow]  - Repeat in-plane: ({nx}, {ny}, 1)")
        atoms_ideal = prim.repeat((nx, ny, 1))
        atoms_ideal.calc = CPUNEP(cfg.nep_model)

        # 生成扰动结构并打力
        print(f"[Workflow]  - Generate {cfg.n_structures} rattled structures")
        structures = generate_mc_rattled_structures(
            atoms_ideal,
            cfg.n_structures,
            cfg.rattle_std,
            cfg.min_dist
        )
        for i, at in enumerate(structures):
            at.calc = CPUNEP(cfg.nep_model)
            _ = at.get_forces()
            if (i + 1) % 50 == 0 or (i + 1) == cfg.n_structures:
                print(f"[Workflow]    forces computed for {i+1}/{cfg.n_structures} structures")

        # 写出训练数据
        prim.wrap()
        write("prim.extxyz", prim)
        write("supercell_ideal.extxyz", atoms_ideal)
        write("supercells_rattled.extxyz", structures)

        print("[Workflow] Step 1 done: prim.extxyz, supercell_ideal.extxyz, supercells_rattled.extxyz")

    # ---------- Step 2: 拟合 FCP ----------
    def fit_force_constants(self):
        cfg = self.cfg
        print("[Workflow] Step 2: fit force constants with HiPhive")

        prim = read("prim.extxyz")
        atoms_ideal = read("supercell_ideal.extxyz")
        rattled = read("supercells_rattled.extxyz", index=":")

        print(f"[Workflow]  - ClusterSpace cutoffs = {cfg.cutoffs}")
        cs = ClusterSpace(prim, cfg.cutoffs)

        sc = StructureContainer(cs)
        for s in prepare_structures(rattled, atoms_ideal):
            sc.add_structure(s)

        print("[Workflow]  - Train Optimizer")
        opt = Optimizer(sc.get_fit_data())
        opt.train()

        parameters = opt.parameters
        parameters_rot = enforce_rotational_sum_rules(cs, parameters, ['Huang', 'Born-Huang'])


        fcp = ForceConstantPotential(cs, parameters_rot)
        fcp.write("fcc-si-film.fcp")

        print("[Workflow] Step 2 done: fcc-si-film.fcp")

    # ---------- Step 3: 导出 FC2/FC3 + phono3py κ ----------
    def compute_kappa(self):
        cfg = self.cfg
        print("[Workflow] Step 3: export FCs and compute kappa with phono3py")

        prim = read("prim.extxyz")

        # Phonopy 超胞
        atoms_ph = PhonopyAtoms(
            symbols=prim.get_chemical_symbols(),
            scaled_positions=prim.get_scaled_positions(),
            cell=prim.cell,
        )

        nx, ny, nz = cfg.dim
        print(f"[Workflow]  - Phonopy supercell dim = ({nx}, {ny}, {nz})")
        phonopy = Phonopy(
            atoms_ph,
            supercell_matrix=np.diag([nx, ny, nz])
        )
        supercell = phonopy.supercell
        supercell_ase = Atoms(
            cell=supercell.cell,
            numbers=supercell.numbers,
            pbc=True,
            scaled_positions=supercell.scaled_positions
        )

        # 从 FCP 导出 FC2/FC3
        print("[Workflow]  - Read fcc-si-film.fcp and export force constants")
        fcp = ForceConstantPotential.read("fcc-si-film.fcp")
        fcs = fcp.get_force_constants(supercell_ase)

        fcs.write_to_phonopy("fc2.hdf5")
        fcs.write_to_phono3py("fc3.hdf5")
        fcs.write_to_phonopy("FORCE_CONSTANTS", format="text")
        fcs.write_to_shengBTE("FORCE_CONSTANTS_3RD", prim)
        write("POSCAR", prim)   # 给 phono3py 用

        print("[Workflow]  - Exported fc2.hdf5, fc3.hdf5, FORCE_CONSTANTS, FORCE_CONSTANTS_3RD")

        # ===== 这里用统一的 mesh 参数 =====
        mx, my, mz = cfg.mesh
        cmd = (
            f'phono3py --dim="{nx} {ny} {nz}" '
            f'--fc2 --fc3 --br '
            f'--mesh="{mx} {my} {mz}" '
            f'--ts="{" ".join(map(str, cfg.temps))}"'
        )
        print("[Workflow]  - Running:", cmd)
        ret = subprocess.call(cmd, shell=True)

        if ret != 0:
            print("[Workflow]  - phono3py exited with code", ret)
            return

        kappa_file = f'kappa-m{mx}{my}{mz}.hdf5'
        if not os.path.isfile(kappa_file):
            print(f"[Workflow]  - kappa file not found: {kappa_file}")
            print("[Workflow]    maybe fc3 is incomplete, or phono3py failed.")
            return

        with h5py.File(kappa_file, "r") as hf:
            temps = hf["temperature"][:]
            print("[Workflow]  - κ(T) temperatures in file:", temps)

        print("[Workflow] Step 3 done, kappa file =", kappa_file)


    # ---------- 一键跑完 ----------
    def run_all(self):
        self.build_training_set()
        self.fit_force_constants()
        self.compute_kappa()
