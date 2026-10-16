"""Shared interpretation of phono3py conductivity datasets."""

import numpy as np


def read_conductivity_components(handle):
    """Read total conductivity and optional verified Wigner components.

    RTA stores the total in ``kappa``. Phono3py 4.1 LBTE stores the particle
    solution there and separately writes ``kappa_intra_exact``. Accept either
    explicit convention, retaining the original HDF5 and rejecting mismatches.
    """
    kappa = np.asarray(handle['kappa'], dtype=float)
    components = {'kappa': kappa}
    wigner_keys = {'kappa_intra', 'kappa_intra_exact', 'kappa_inter'} & set(handle)
    if not wigner_keys:
        return components
    particle_key = 'kappa_intra_exact' if 'kappa_intra_exact' in handle else 'kappa_intra'
    if particle_key not in handle or 'kappa_inter' not in handle:
        raise ValueError('Incomplete Wigner components')
    intra = np.asarray(handle[particle_key], dtype=float)
    inter = np.asarray(handle['kappa_inter'], dtype=float)
    if intra.shape != kappa.shape or inter.shape != kappa.shape:
        raise ValueError('Wigner conductivity components have inconsistent shapes')
    if not all(np.isfinite(values).all() for values in (kappa, intra, inter)):
        raise ValueError('Nonfinite Wigner conductivity components')
    total = intra + inter
    stores_particle = particle_key == 'kappa_intra_exact' and np.allclose(
        kappa, intra, rtol=5e-4, atol=5e-4
    )
    if not stores_particle and not np.allclose(kappa, total, rtol=5e-4, atol=5e-4):
        raise ValueError('Wigner total conductivity differs from particle plus coherence')
    return {'kappa': total, 'kappa_intra': intra, 'kappa_inter': inter}
