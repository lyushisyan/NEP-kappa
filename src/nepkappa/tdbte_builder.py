"""Build nonlinear energy-shell chunks directly from phono3py force constants.

Native shortest-vector interaction strengths are evaluated one parent q at a
time. There is no global Nq² vertex array, finite-cell Fourier approximation,
Gaussian energy broadening, or dense collision matrix in this route.
"""
import itertools
import json
from pathlib import Path

import numpy as np

from .tdbte import ShellKernel
from .tdbte_quadrature import section_quadrature
from .tdbte_storage import ChunkWriter, atomic_json, sha256

INPUT_FILES = ('phono3py_disp.yaml', 'fc2.hdf5', 'fc3.hdf5')
BUILDER_VERSION = 1
CHUNK_EVENTS = 20000
CUTOFF_THz = 1e-4


def source_hashes(source):
    source = Path(source)
    missing = [name for name in INPUT_FILES if not (source / name).is_file()]
    if missing:
        raise FileNotFoundError(f'TD-BTE force-constant directory {source} lacks: {", ".join(missing)}')
    return {name: sha256(source / name) for name in INPUT_FILES}


def validate_mesh(mesh):
    if (not isinstance(mesh, (list, tuple, np.ndarray)) or len(mesh) != 3
            or any(isinstance(n, (bool, np.bool_)) or not isinstance(n, (int, np.integer))
                   or n < 2 for n in mesh)):
        raise ValueError('tdbte.mesh requires three integers >= 2 (full 3D Gamma-centered grid)')
    return np.asarray(mesh, dtype=int)


def _geometry(addresses, mesh):
    lookup = {tuple(a): i for i, a in enumerate(addresses)}
    if len(lookup) != int(np.prod(mesh)):
        raise ValueError('Expected a full uniform grid without symmetry reduction')

    def indices(a):
        return np.array([lookup[tuple(row % mesh)] for row in a.reshape(-1, 3)],
                        dtype=int).reshape(a.shape[:-1])

    local = []
    for permutation in itertools.permutations(range(3)):
        vertices = np.zeros((4, 3), int)
        for k, axis in enumerate(permutation):
            vertices[k + 1] = vertices[k]
            vertices[k + 1, axis] = 1
        local.append(vertices)
    local = np.asarray(local)
    q1 = indices(addresses[:, None, None, :] + local[None, :, :, :]).reshape(-1, 4)
    return indices, q1, np.tile(local, (len(addresses), 1, 1))


def _prepare_modes(pp, mesh):
    pp.run_phonon_solver()
    grid = pp.bz_grid
    # SNF/generalized grids require a transformed cell measure and different
    # indexing. Reject rather than silently using rectangular-grid geometry.
    if not np.array_equal(grid.D_diag, mesh) or not np.array_equal(grid.Q, np.eye(3, dtype=int)):
        raise ValueError('TD-BTE requires a diagonal, unshifted regular q grid')
    addresses = grid.addresses[grid.grg2bzg] % mesh
    original = pp.get_phonons()[0][grid.grg2bzg].copy()
    if original.ndim != 2 or not np.isfinite(original).all():
        raise ValueError('Invalid native phonon frequencies')
    nb = original.shape[1]
    zero = abs(original) <= CUTOFF_THz
    expected = np.all(addresses == 0, axis=1)[:, None] & (np.arange(nb)[None, :] < 3)
    if nb < 3 or not np.array_equal(zero, expected) or np.any(original[~zero] <= 0):
        raise ValueError('Only three near-zero Gamma translations are allowed; '
                         'check FC2 stability and acoustic sum rules before TD-BTE')
    frequency = original.copy()
    frequency[zero] = 0
    return addresses, frequency, zero, original


def _native_strength(pp, parent, inverse, addresses, mesh, frequency):
    """Map native q0+q1+q2=G to decay parent=-q0; keep ordered daughters."""
    grid = pp.bz_grid
    pp.set_grid_point(int(grid.grg2bzg[inverse[parent]]))
    pp.run()
    triplets, weights = pp.get_triplets_at_q()[:2]
    points = grid.bzg2grg[triplets]
    nq, nb = frequency.shape
    if (len(points) != nq or not np.all(weights == 1)
            or not np.all(points[:, 0] == inverse[parent])
            or len(np.unique(points[:, 1])) != nq):
        raise ValueError('Native backend did not provide every ordered full-grid triplet')
    if np.any((addresses[parent] - addresses[points[:, 1]] - addresses[points[:, 2]]) % mesh):
        raise ValueError('Native triplets violate momentum conservation')
    native = pp.interaction_strength
    if native.shape != (nq, nb, nb, nb) or not np.isfinite(native).all() or np.any(native < 0):
        raise ValueError('Invalid native interaction strengths')
    strength = np.empty_like(native)
    strength[points[:, 1]] = native
    # Match the explicitly recorded rigid-translation acoustic-limit boundary:
    # zero-energy vertex legs have zero strength, but shell interiors remain.
    for _, q1, q2 in points:
        block = strength[q1]
        block[frequency[parent] == 0, :, :] = 0
        block[:, frequency[q1] == 0, :] = 0
        block[:, :, frequency[q2] == 0] = 0
    return strength


def build_kernel(source, mesh, directory, *, chunk_events=CHUNK_EVENTS, excited_branches=None):
    """Build or reuse a complete kernel, checking FC hashes and mesh on reuse.

    Internal chunk_events is exposed to tests, not a required user parameter.
    Partial builds are preserved for diagnosis and never automatically resumed.
    """
    import phono3py
    from phono3py.phonon3.imag_self_energy import ImagSelfEnergy

    mesh = validate_mesh(mesh)
    if not isinstance(chunk_events, int) or chunk_events < 1:
        raise ValueError('chunk_events must be a positive integer')
    source, directory = Path(source).resolve(), Path(directory)
    if (source / 'BORN').exists():
        raise ValueError('NAC/BORN inputs are not yet supported by the TD-BTE builder')
    hashes = source_hashes(source)
    path = directory / 'manifest.json'
    if directory.exists():
        if not path.is_file():
            raise FileExistsError(f'Incomplete TD-BTE kernel at {directory}; preserve it and choose a new result-dir')
        meta = json.loads(path.read_text())
        if (meta.get('builder_version') != BUILDER_VERSION or meta.get('complete') is not True
                or meta.get('input_sha256') != hashes or meta.get('mesh') != mesh.tolist()
                or meta.get('phono3py_version') != phono3py.__version__):
            raise ValueError('Cached TD-BTE kernel differs in inputs, mesh or builder/backend version; use a new result-dir')
        return path  # ChunkedKernel performs checksum and shell validation.

    ph = phono3py.load(source / INPUT_FILES[0], fc2_filename=source / INPUT_FILES[1],
                      fc3_filename=source / INPUT_FILES[2], produce_fc=False,
                      is_mesh_symmetry=False, symmetrize_fc=False, log_level=0)
    if ph.nac_params is not None:
        raise ValueError('Automatic TD-BTE kernel construction does not yet support NAC; '
                         'do not silently discard Born charges for a polar calculation')
    ph.mesh_numbers = mesh
    ph.init_phph_interaction()
    pp = ph.phph_interaction
    if pp.cutoff_frequency != CUTOFF_THz:
        raise ValueError('Native phonon cutoff differs from the TD-BTE acoustic boundary')
    addresses, frequency, zero, original = _prepare_modes(pp, mesh)
    nq, nb = frequency.shape
    if excited_branches is not None and any(b >= nb for b in excited_branches):
        raise ValueError(f'Excited branch index is absent: this structure has {nb} branches')
    indices, q1, geometry = _geometry(addresses, mesh)
    inverse = indices(-addresses)
    if not np.allclose(frequency, frequency[inverse], rtol=1e-10, atol=1e-8):
        raise ValueError('Native phonon frequencies violate time-reversal symmetry')
    # Native strengths include 1/Nq. Ordered daughters need no additional half
    # with this conversion; see tagged-leg linewidth parity regression tests.
    conversion = float(4 * np.pi * ImagSelfEnergy(pp).unit_conversion_factor)
    if not np.isfinite(conversion) or conversion <= 0:
        raise ValueError('Invalid phono3py rate conversion')
    metadata = {
        'builder_version': BUILDER_VERSION, 'phono3py_version': phono3py.__version__,
        'model': 'native shortest-vector vertices; linear shells; entropy interpolation; adjoint deposition',
        'coefficient_time_unit': 'ps', 'coefficient_conversion': conversion,
        'native_cutoff_frequency_THz': float(pp.cutoff_frequency),
        'mesh': mesh.tolist(), 'uniform_mode_weight': 1 / nq,
        'mode_branch_indices': np.tile(np.arange(nb), nq).tolist(),
        'grid_addresses': addresses.tolist(), 'source': str(source), 'input_sha256': hashes,
        'quadrature': 'six Freudenthal tetrahedra/cell; symmetric-degree2',
        'normalization': 'ordered daughters; 4*pi*U*strength*surface weight; strength includes 1/Nq',
        'boundary_assumptions': 'three fixed Gamma translations; eta_Gamma=0; strength_Gamma=0; acoustic limit not off-grid validated',
        'max_rigid_frequency_adjustment_THz': float(abs(original[zero]).max()),
        'source_sum_rule_diagnostics': {
            'fc2_sum_second_atom_max': float(abs(ph.fc2.sum(axis=1)).max()),
            'fc3_sum_second_atom_max': float(abs(ph.fc3.sum(axis=1)).max()),
            'fc3_sum_third_atom_max': float(abs(ph.fc3.sum(axis=2)).max()),
            'note': 'Compact-array partial checks; inputs are not modified or symmetrized'},
        'validated_for_physical_dynamics': False,
        'limitations': 'frequency-sorted branches; acoustic endpoints and mesh/surface quadrature convergence need validation; not the native phono3py tetrahedron operator',
    }
    print(f'TD-BTE kernel: {nq} q points, {nb} branches; native vertex block '
          f'{nq * nb**3 * 8 / 2**20:.1f} MiB; streaming event chunks', flush=True)
    writer = ChunkWriter(directory, frequency.ravel(), zero.ravel(), metadata)
    buffers = [[] for _ in range(5)]
    parents, firsts, seconds, barys, rates = buffers
    max_residual = 0.

    def flush():
        nonlocal max_residual
        if not parents:
            return
        k = ShellKernel(frequency.ravel(), np.asarray(parents, dtype=np.int64),
                        np.asarray(firsts, dtype=np.int64), np.asarray(seconds, dtype=np.int64),
                        np.asarray(barys), np.asarray(rates), zero.ravel())
        writer.append(k)
        residual = k.frequency[k.parent] - (k.weights * (
            k.frequency[k.daughter1] + k.frequency[k.daughter2])).sum(axis=1)
        max_residual = max(max_residual, float(abs(residual).max()))
        for buffer in buffers:
            buffer.clear()

    try:
        for p in range(nq):
            strength = _native_strength(pp, p, inverse, addresses, mesh, frequency)
            q2 = indices(addresses[p] - addresses[q1])
            for a, b, c in itertools.product(range(nb), repeat=3):
                if zero[p, a]:
                    continue
                detuning = frequency[p, a] - frequency[q1, b] - frequency[q2, c]
                hit = (detuning.min(axis=1) <= 0) & (detuning.max(axis=1) >= 0)
                for it in np.flatnonzero(hit):
                    # Reject singular resonant volumes, even if grid vertices
                    # happen to have zero strength. Never publish partial sums.
                    nodes = section_quadrature(geometry[it], detuning[it])
                    vertex_strength = strength[q1[it], a, b, c]
                    for w, surface_weight in nodes:
                        rate = conversion * (w @ vertex_strength) * surface_weight
                        if rate <= 0:
                            continue
                        parents.append(p * nb + a)
                        firsts.append(q1[it] * nb + b)
                        seconds.append(q2[it] * nb + c)
                        barys.append(w)
                        rates.append(rate)
                        if len(parents) >= chunk_events:
                            flush()
            print(f'TD-BTE kernel: parent q {p + 1}/{nq}; {len(writer.chunks)} chunks saved', flush=True)
        flush()
        if source_hashes(source) != hashes:
            raise ValueError('Force-constant inputs changed during kernel construction')
        writer.metadata['max_shell_residual_THz'] = max_residual
        return writer.finish()
    except Exception as exc:
        atomic_json(directory / 'build-failure.json', {
            'complete': False, 'error': str(exc), 'input_sha256': hashes,
            'mesh': mesh.tolist(), 'completed_chunks': len(writer.chunks)})
        raise
