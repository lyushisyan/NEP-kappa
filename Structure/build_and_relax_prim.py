#!/usr/bin/env python3
# build_and_relax_prim.py

import os
import argparse
import numpy as np

from ase import Atoms
from ase.build import bulk, diamond100
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


def make_si_bulk() -> Atoms:
    """
    3D bulk Si primitive cell (no vacuum).
    """
    return bulk("Si")


# ===================== Main =====================
def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--model",
        type=str,
        choices=["film", "wire", "bulk"],
        required=True,
        help=(
            "Structure type: "
            "'film' = Si(001) slab; "
            "'wire' = 1D Si nanowire; "
            "'bulk' = 3D bulk Si."
        ),
    )

    parser.add_argument(
        "--nep",
        type=str,
        default=None,
        help=(
            "NEP potential filename in NEP/. "
            "If not set: film/wire/bulk → Si_2025_Xuke.txt."
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
            "ignored for bulk."
        ),
    )

    parser.add_argument(
        "--vac",
        type=float,
        default=15.0,
        help="Vacuum thickness in Å (for film: along z; for wire: transverse directions; ignored for bulk).",
    )

    args = parser.parse_args()

    model = args.model
    thick = args.thick
    vac = args.vac

    # -------- 根据 model 选择 NEP 默认文件 --------
    if args.nep is None:
        nep_name = "Si_2025_Xuke.txt"
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
    elif model == "bulk":
        print("==> bulk mode: using primitive cell; --thick and --vac are ignored.")
        prim = make_si_bulk()
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

    output_poscar = os.path.join(os.path.dirname(__file__), "POSCAR")
    write(output_poscar, prim, format="vasp")
    print(f"==> Saved POSCAR to {output_poscar}")


if __name__ == "__main__":
    main()
