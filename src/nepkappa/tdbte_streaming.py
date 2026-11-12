"""Conservative first-order finite-volume advection on a periodic 1D grid."""
import numpy as np


def face_flux(occupation, velocity):
    """Flux at each cell's RIGHT face; the last face wraps to cell zero."""
    n, v = np.asarray(occupation), np.asarray(velocity)
    return np.maximum(v, 0) * n + np.minimum(v, 0) * np.roll(n, -1, axis=0)


def divergence(flux, dx):
    return (flux - np.roll(flux, 1, axis=0)) / dx


def stream(occupation, velocity, dt, dx):
    n, v = np.asarray(occupation), np.asarray(velocity)
    if (n.ndim != 2 or v.shape != (n.shape[1],) or dx <= 0 or dt < 0
            or not np.isfinite([dx, dt]).all() or not np.isfinite(n).all()
            or not np.isfinite(v).all() or np.any(n < 0)):
        raise ValueError('Invalid streaming state, velocity or grid')
    if dt * np.max(abs(v), initial=0.) / dx > 1 + 1e-14:
        raise ValueError('Upwind streaming CFL exceeds one')
    flux = face_flux(n, v)
    # Convex-combination form avoids cancellation near CFL=1. It is exactly
    # equivalent to the conservative divergence formula, without clipping.
    c = dt * abs(v) / dx
    upstream = np.where(v >= 0, np.roll(n, 1, axis=0), np.roll(n, -1, axis=0))
    result = (1 - c) * n + c * upstream
    if np.any(result < 0) or not np.isfinite(result).all():
        raise RuntimeError('Streaming lost nonnegative finite occupations')
    return result, flux
