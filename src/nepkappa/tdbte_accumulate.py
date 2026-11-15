"""Optional compiled block reduction; NumPy remains the reference backend."""
import numpy as np

try:
    from numba import njit, prange, get_num_threads, set_num_threads, config
except ImportError:
    njit = None


def _accumulate(n, eta, parent, first, second, weights, coefficients, rhs):
    for event in range(len(parent)):
        p = parent[event]
        eb, ec = 0., 0.
        for j in range(weights.shape[1]):
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
        for j in range(weights.shape[1]):
            value = weights[event, j] * flux
            rhs[first[event, j]] += value
            rhs[second[event, j]] += value


# Do not write JIT caches beside an installed (possibly read-only) package.
# Never enable fastmath: detailed-balance cancellation is part of the model.
accumulate = njit(cache=False)(_accumulate) if njit is not None else None


def _accumulate_batch(n, eta, parent, first, second, weights, coefficients, rhs):
    for cell in range(len(n)):
        accumulate(n[cell], eta[cell], parent, first, second, weights, coefficients, rhs[cell])


accumulate_batch = njit(cache=False)(_accumulate_batch) if njit is not None else None


def _accumulate_batch_parallel(n, eta, parent, first, second, weights, coefficients, rhs):
    # Each worker owns one RHS row; event deposition within that row stays serial.
    for cell in prange(len(n)):
        accumulate(n[cell], eta[cell], parent, first, second, weights, coefficients, rhs[cell])


accumulate_batch_parallel = njit(parallel=True, cache=False)(_accumulate_batch_parallel) if njit is not None else None


def validate_collision_threads(threads, kernel):
    if isinstance(threads, bool) or not isinstance(threads, (int, np.integer)) or threads < 1:
        raise ValueError('Collision threads must be a positive integer')
    if threads > 1:
        if njit is None or getattr(kernel, 'accumulate', accumulate) is None:
            raise ValueError('Multicore collisions require the Numba backend')
        if threads > config.NUMBA_NUM_THREADS:
            raise ValueError(f'Collision threads exceed Numba limit {config.NUMBA_NUM_THREADS}; '
                             'set NUMBA_NUM_THREADS before launching Python')


def collision_batch(kernel, occupation, batch_size=16, threads=1):
    """Use independent cell workers, restoring the caller's Numba thread count."""
    validate_collision_threads(threads, kernel)
    if threads == 1:
        return _collision_batch(kernel, occupation, batch_size, parallel=False)
    previous = get_num_threads()
    try:
        set_num_threads(threads)
        return _collision_batch(kernel, occupation, batch_size, parallel=True)
    finally:
        set_num_threads(previous)


def _collision_batch(kernel, occupation, batch_size, parallel):
    """Evaluate the unchanged event flux on cells, opening each chunk once.

    Only a cell batch's entropy variables and one event chunk are temporary.
    The result and input have shape (cells, modes); no cell-by-event array or
    dense scattering matrix is formed.
    """
    from contextlib import nullcontext
    n = np.asarray(occupation, dtype=float)
    if (n.ndim != 2 or n.shape[1:] != kernel.frequency.shape or not len(n)
            or not np.isfinite(n).all() or np.any(n[:, ~kernel.zero] <= 0)):
        raise ValueError('Active occupations must be finite and positive; expected (cells, modes)')
    if isinstance(batch_size, bool) or not isinstance(batch_size, (int, np.integer)) or batch_size < 1:
        raise ValueError('Collision batch size must be a positive integer')
    rhs = np.zeros_like(n)
    chunks = getattr(kernel, 'chunks', [None])
    compiled = getattr(kernel, 'accumulate', accumulate) is not None
    for entry in chunks:
        with (kernel._read(entry) if entry is not None else nullcontext(kernel)) as part:
            for start in range(0, len(n), batch_size):
                state = n[start:start + batch_size]
                target = rhs[start:start + batch_size]
                if compiled:
                    eta = np.zeros_like(state)
                    eta[:, ~kernel.zero] = np.log1p(state[:, ~kernel.zero]) - np.log(state[:, ~kernel.zero])
                    reducer = accumulate_batch_parallel if parallel and len(state) > 1 else accumulate_batch
                    reducer(state, eta, part.parent, part.daughter1, part.daughter2,
                                     part.weights, part.coefficients, target)
                else:
                    for cell in range(len(state)):
                        target[cell] += part.collision(state[cell])[0]
    rhs[:, kernel.zero] = 0
    if not np.isfinite(rhs).all():
        raise ValueError('Nonfinite accumulated collision derivative')
    return rhs
