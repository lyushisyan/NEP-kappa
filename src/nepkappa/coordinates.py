"""Shared unit-cell convention for force-constant generation.

Canonicalize *before* creating supercells or IFCs. Changing only the structure
afterwards changes the lattice-translation labels of those IFCs.
"""

import numpy as np


COORDINATE_TOLERANCE = 1e-10


def canonicalize_atoms(atoms):
    """Return a copy with periodic fractional coordinates in [0, 1).

    Snap floating-point boundary noise to zero, without sorting atoms, changing
    the lattice, or touching non-periodic directions. This must not be applied
    to displaced supercells independently of their force-constant indices.
    """
    result = atoms.copy()
    positions = atoms.get_scaled_positions(wrap=False)
    for axis, periodic in enumerate(atoms.pbc):
        if periodic:
            values = np.remainder(positions[:, axis], 1.0)
            boundary = (values < COORDINATE_TOLERANCE) | (
                1.0 - values < COORDINATE_TOLERANCE
            )
            values[boundary] = 0.0
            positions[:, axis] = values
    result.set_scaled_positions(positions)
    return result
