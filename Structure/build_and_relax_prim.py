#!/usr/bin/env python3
# build_and_relax_prim.py

import os
import argparse
import numpy as np

from ase import Atoms
from ase.build import bulk, diamond100, graphene_nanoribbon
from ase.io import write

from calorine.calculators import CPUNEP
from calorine.tools import relax_structure

NEP_DIR = os.path.join(os.path.dirname(__file__), "../NEP")


# ===================== Helper functions =====================
def make_si_thinfilm(thick_n: int, vac_z: float) -> Atoms:
    """
    Si(001) film with a simple surface reconstruction.
    """
    Nz = thick_n
    slab = diamond100("Si", size=(2, 2, Nz), vacuum=vac_z, orthogonal=True)

    positions = slab.get_positions()
    z_coords = positions[:, 2]

    # find top / bottom layers
    tol = 0.2  # Å
    zmin = z_coords.min()
    zmax = z_coords.max()

    bottom_indices_all = np.where(np.abs(z_coords - zmin) < tol)[0]
    top_indices_all = np.where(np.abs(z_coords - zmax) < tol)[0]

    bottom_indices = list(bottom_indices_all[:4])
    top_indices = list(top_indices_all[:4])

    y_offsets = [0.5, 0.9, -0.9, -0.5]

    # top surface (move inward: z - 0.3)
    for i, idx in enumerate(top_indices):
        positions[idx, 1] += y_offsets[i % 4]
        positions[idx, 2] -= 0.3

    # bottom surface (move inward: z + 0.3)
    for i, idx in enumerate(bottom_indices):
        positions[idx, 1] += y_offsets[i % 4]
        positions[idx, 2] += 0.3

    slab.set_positions(positions)
    return slab


def make_si_wire(thick_n: int, vac_xy: float) -> Atoms:
    """
    1D Si nanowire along z with vacuum in x,y.
    thick_n: repeats in x and y (cross section ~ thick_n * a0).
    """
    uc = bulk("Si", cubic=True)
    a0 = uc.cell.lengths()[0]

    base = uc.repeat((thick_n, thick_n, 1))
    pos = base.get_positions()

    Lx = thick_n * a0 + 2.0 * vac_xy
    Ly = thick_n * a0 + 2.0 * vac_xy

    # center in x,y
    pos[:, 0] += vac_xy
    pos[:, 1] += vac_xy

    wire = Atoms(
        numbers=base.numbers,
        positions=pos,
        cell=[(Lx, 0.0, 0.0), (0.0, Ly, 0.0), (0.0, 0.0, a0)],
        pbc=(True, True, True),  # x,y 有大真空，z 周期 → 纳米线
    )
    return wire


def make_agnr(N: int,
              vacuum_y: float = 10.0,
              vacuum_z: float = 5.0,
              x_len: int = 1) -> Atoms:
    """
    构造 armchair GNR (AGNR)：
      - x 方向为带长方向，周期长度 Lx = 3 * a_cc
      - N 控制宽度（dimer lines 数）
      - y,z 用 vacuum_y / vacuum_z 控制真空

    返回 ASE Atoms，用于后续 NEP relax + 写 POSCAR。
    """
    # 一个长度周期内的 4 个原子（单位：Å）
    a_cc = 1.42               # C-C bond length (Å)
    s = np.sqrt(3.0) * 0.5 * a_cc   # vertical projection

    Lx = 3.0 * a_cc
    base_hex4 = np.array([
        [0.5 * a_cc, 0.0, 0.0],
        [1.0 * a_cc,   s, 0.0],
        [2.0 * a_cc,   s, 0.0],
        [2.5 * a_cc, 0.0, 0.0],
    ], float)

    # 宽度方向叠 N 个单元
    coords = []
    for j in range(N):
        blk = base_hex4.copy()
        blk[:, 1] += 2.0 * j * s   # 每层往 +y 方向平移 2s
        coords.append(blk)
    coords = np.vstack(coords)     # 现在有 4N 个原子（一个长度周期）

    # 再沿 x 方向复制 replicate_len 个周期
    coords_all = []
    for i in range(x_len):
        blk = coords.copy()
        blk[:, 0] += i * Lx
        coords_all.append(blk)
    coords_all = np.vstack(coords_all)

    # 盒子尺寸：x 周期，y/z 真空
    Lx_total = Lx * x_len
    Ly = N*2*s + 2.0 * vacuum_y
    Lz = 2.0 * vacuum_z

    # 把原子整体平移到盒子中央（y,z方向）
    coords_all[:, 1] += 0.5 * Ly
    coords_all[:, 2] += 0.5 * Lz

    natoms = coords_all.shape[0]
    symbols = ["C"] * natoms

    gnr = Atoms(
        symbols=symbols,
        positions=coords_all,
        cell=[(Lx_total, 0.0, 0.0),
              (0.0,       Ly, 0.0),
              (0.0,      0.0, Lz)],
        pbc=(True, True, True),  # x 周期，y/z 大真空
    )
    return gnr

def make_zgnr(N: int,
              vacuum_y: float = 10.0,
              vacuum_z: float = 5.0,
              x_len: int = 2) -> Atoms:
    """
    构造 armchair GNR (ZGNR)：
      - x 方向为带长方向，周期长度 Lx = 2 * s
      - N 控制宽度（dimer lines 数）
      - y,z 用 vacuum_y / vacuum_z 控制真空

    返回 ASE Atoms，用于后续 NEP relax + 写 POSCAR。
    """
    # 一个长度周期内的 4 个原子（单位：Å）
    a_cc = 1.42               # C-C bond length (Å)
    s = np.sqrt(3.0) * 0.5 * a_cc   # vertical projection

    Lx = 2 * s
    base_hex4 = np.array([
        [0, 0.5 * a_cc, 0.0],
        [s, 1.0 * a_cc, 0.0],
        [s, 2.0 * a_cc, 0.0],
        [0, 2.5 * a_cc, 0.0],
    ], float)

    # 宽度方向叠 N 个单元
    coords = []
    for j in range(N):
        blk = base_hex4.copy()
        blk[:, 1] += 3.0 * j * a_cc   # 每层往 +y 方向平移 3a
        coords.append(blk)
    coords = np.vstack(coords)     # 现在有 4N 个原子（一个长度周期）

    # 再沿 x 方向复制 replicate_len 个周期
    coords_all = []
    for i in range(x_len):
        blk = coords.copy()
        blk[:, 0] += i * Lx
        coords_all.append(blk)
    coords_all = np.vstack(coords_all)

    # 盒子尺寸：x 周期，y/z 真空
    Lx_total = Lx * x_len
    Ly = N*2*s + 2.0 * vacuum_y
    Lz = 2.0 * vacuum_z

    # 把原子整体平移到盒子中央（y,z方向）
    coords_all[:, 1] += 0.5 * Ly
    coords_all[:, 2] += 0.5 * Lz

    natoms = coords_all.shape[0]
    symbols = ["C"] * natoms

    gnr = Atoms(
        symbols=symbols,
        positions=coords_all,
        cell=[(Lx_total, 0.0, 0.0),
              (0.0,       Ly, 0.0),
              (0.0,      0.0, Lz)],
        pbc=(True, True, True),  # x 周期，y/z 大真空
    )
    return gnr

# ===================== Main =====================
def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--model",
        type=str,
        choices=["film", "wire", "agnr", "zgnr"],
        required=True,
        help=(
            "Structure type: "
            "'film' = Si(001) slab; "
            "'wire' = 1D Si nanowire; "
            "'agnr' = armchair graphene nanoribbon; "
            "'zgnr' = zigzag graphene nanoribbon."
        ),
    )

    parser.add_argument(
        "--nep",
        type=str,
        default=None,
        help=(
            "NEP potential filename in NEP/. "
            "If not set: film/wire → Si_2025_xuke.txt; "
            "agnr/zgnr → C_2024_liangting.txt."
        ),
    )

    parser.add_argument(
        "--thick",
        type=int,
        default=8,
        help=(
            "Thickness/width parameter (default = 8): "
            "film = number of Si layers along z; "
            "wire = repeats in x,y; "
            "agnr/zgnr = GNR width parameter n."
        ),
    )

    parser.add_argument(
        "--vac",
        type=float,
        default=15.0,
        help="Vacuum thickness in Å (for film: along z; for wire/GNR: transverse directions).",
    )

    args = parser.parse_args()

    model = args.model
    thick = args.thick
    vac = args.vac

    # -------- 根据 model 选择 NEP 默认文件 --------
    if args.nep is None:
        if model in ("film", "wire"):
            nep_name = "Si_2025_xuke.txt"
        else:  # agnr / zgnr
            nep_name = "C_2024_liangting.txt"
    else:
        nep_name = args.nep
    
    nep_file = os.path.join(NEP_DIR, nep_name)

    print(f"==> MODEL : {model}")
    print(f"==> NEP   : {nep_file}")
    print(f"==> THICK : {thick}")
    print(f"==> VAC   : {vac:.3f} Å")

    if model == "film":
        prim = make_si_thinfilm(thick, vac)
    elif model == "wire":
        prim = make_si_wire(thick, vac)
    elif model == "agnr":
        prim = make_agnr(N=thick, vacuum_y=vac, x_len=1)
    elif model == "zgnr":
        prim = make_zgnr(N=thick, vacuum_y=vac, x_len=2)
    else:
        raise ValueError(f"Unknown model: {model}")

    print("==> NEP relaxation ...")
    prim.calc = CPUNEP(nep_file)

    for fmax in [1e-3, 1e-4, 1e-5]:
        print(f"Relaxing structure (fmax = {fmax:.0e}) ...")
        relax_structure(prim, fmax=fmax)
        f = prim.get_forces()
        max_f = np.max(np.abs(f))
        print(f"max|F| = {max_f:8.5f} eV/Å")

    print("Relaxation finished.")

    write("POSCAR", prim, format="vasp")
    print("==> Saved POSCAR")


if __name__ == "__main__":
    main()
