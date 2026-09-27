"""Optional compiled block reduction; NumPy remains the reference backend."""
import numpy as np

try:
    from numba import njit
except ImportError:
    njit = None


def _accumulate(n, eta, parent, first, second, weights, coefficients, rhs):
    for event in range(len(parent)):
        p = parent[event]
        eb, ec = 0., 0.
        for j in range(4):
            eb += weights[event, j] * eta[first[event, j]]
            ec += weights[event, j] * eta[second[event, j]]
        if eb <= 0 or ec <= 0:
            raise ValueError('Singular daughter entropy')
        if coefficients[event] == 0:
            continue
        lb = -eb - np.log(-np.expm1(-eb))
        lc = -ec - np.log(-np.expm1(-ec))
        affinity = -eta[p] + eb + ec
        scale = np.log1p(n[p]) + lb + lc + max(affinity, 0.) + np.log(coefficients[event])
        flux = np.exp(scale) * np.sign(affinity) * (-np.expm1(-abs(affinity)))
        rhs[p] -= flux
        for j in range(4):
            value = weights[event, j] * flux
            rhs[first[event, j]] += value
            rhs[second[event, j]] += value


# Do not write JIT caches beside an installed (possibly read-only) package.
# Never enable fastmath: detailed-balance cancellation is part of the model.
accumulate = njit(cache=False)(_accumulate) if njit is not None else None
