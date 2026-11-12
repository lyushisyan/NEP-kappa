"""Spatial numerical checks with a genuinely nonlinear, energy-conserving toy.

The toy has six modes and both 2<->1+1 and 3<->1+2 scattering. Coefficients
are arbitrary. These tests establish discretization behavior, not material rates.
"""
from contextlib import contextmanager
from itertools import product
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from nepkappa.tdbte import ShellKernel
from nepkappa.tdbte_storage import ChunkWriter, ChunkedKernel
from nepkappa.tdbte_observables import Observables, temperature_harmonics
from nepkappa.tdbte_spatial import (
    spatial_options, initial_occupation, periodic_gaussian, evolve_spatial,
)
from nepkappa.tdbte_streaming import stream, divergence


def toy_kernel(coefficient=.05):
    events = []
    for branch, first, second in ((1, 0, 0), (2, 0, 1)):
        for p, b, c in product(range(2), repeat=3):
            events.append((p*3+branch, b*3+first, c*3+second))
    events = np.asarray(events)
    k = ShellKernel(np.tile([1., 2., 3.], 2), events[:, 0],
        np.repeat(events[:, 1:2], 4, axis=1), np.repeat(events[:, 2:3], 4, axis=1),
        np.full((len(events), 4), .25), np.full(len(events), coefficient), np.zeros(6, bool))
    k.validate()
    return k


def setup(nx=16, length=8., initial=None, collisions=True):
    k = toy_kernel()
    v = np.repeat([1., -1.], 3)
    obs = Observables(k, np.tile(np.arange(3), 2), .5, 1., v)
    cfg = SimpleNamespace(tdbte_spatial={'length-nm': length, 'cells': nx,
        'collisions': collisions, 'initial': initial or {'kind': 'ttg', 'delta-temperature-K': 3}},
        tdbte_temperature=300., tdbte_excitation=.01, tdbte_branches=[2],
        tdbte_duration_ps=.1, tdbte_mesh=[3, 3, 3])
    return k, v, obs, spatial_options(cfg), (np.arange(nx)+.5)*length/nx


@pytest.mark.parametrize('compiled', [False, True])
@pytest.mark.parametrize('batch', [1, 3, 16])
def test_batched_collision_opens_chunks_once_and_matches_scalar(tmp_path, compiled, batch):
    from nepkappa.tdbte_accumulate import accumulate
    if compiled and accumulate is None:
        pytest.skip('Compiled reduction is unavailable in this backend run')
    k = toy_kernel()
    writer = ChunkWriter(tmp_path/'kernel', k.frequency, k.zero, {})
    for start, stop in ((0, 7), (7, 16)):
        writer.append(ShellKernel(k.frequency, k.parent[start:stop], k.daughter1[start:stop],
            k.daughter2[start:stop], k.weights[start:stop], k.coefficients[start:stop], k.zero))
    chunked = ChunkedKernel(writer.finish())
    if not compiled:
        chunked.accumulate = None
    reads = []
    original = chunked._read
    @contextmanager
    def read(entry):
        reads.append(entry)
        with original(entry) as part:
            yield part
    chunked._read = read
    n = k.bose(300) * np.exp(np.random.default_rng(10).uniform(-1, 1, (7, 6)))
    expected = np.array([k.collision(row)[0] for row in n])
    np.testing.assert_allclose(chunked.collision_batch(n, batch), expected, rtol=1e-12, atol=1e-12)
    assert len(reads) == len(chunked.chunks)
    np.testing.assert_allclose(k.collision_batch(n, batch), expected, rtol=1e-12, atol=1e-12)


def test_upwind_is_positive_conservative_and_uniform_stationary():
    n = np.random.default_rng(3).uniform(.1, 3, (31, 6))
    v = np.array([1., -1., 0., .5, -.3, .7])
    updated, flux = stream(n, v, .2, .2)
    assert np.min(updated) > 0
    np.testing.assert_allclose(updated.sum(axis=0), n.sum(axis=0), rtol=1e-14)
    np.testing.assert_allclose(updated-n, -.2*divergence(flux, .2), atol=5e-16)
    uniform = np.tile(n[0], (31, 1))
    np.testing.assert_array_equal(stream(uniform, v, .1, .2)[0], uniform)
    with pytest.raises(ValueError, match='CFL'):
        stream(n, v, .3, .2)


def test_ballistic_translation_and_first_order_spatial_convergence():
    errors = []
    for nx in (32, 64, 128):
        dx = 1/nx
        x = (np.arange(nx)+.5)*dx
        n = np.stack([2 + np.sin(2*np.pi*x), 2 + np.cos(2*np.pi*x)], axis=-1)
        v = np.array([1., -.5])
        time = 0.
        while time < .2:
            dt = min(.4*dx, .2-time)
            n = stream(n, v, dt, dx)[0]
            time += dt
        exact = np.stack([2+np.sin(2*np.pi*(x-.2)), 2+np.cos(2*np.pi*(x+.1))], axis=-1)
        errors.append(np.linalg.norm(n-exact)/np.sqrt(nx))
    assert errors[0]/errors[1] > 1.8
    assert errors[1]/errors[2] > 1.8


def test_single_cell_reproduces_homogeneous_nonlinear_solver():
    k, v, obs, options, x = setup(nx=1, initial={'kind': 'selective', 'profile': 'uniform'})
    n = initial_occupation(obs, x, 8, 300, options['initial'])
    saved = []
    stats = evolve_spatial(k, n, v, 8, [0., .05, .1], .01, options, obs,
                          lambda i, t, n, f: saved.append(n.copy()))
    reference = solve_ivp(lambda t, n: k.collision(n)[0], (0, .1), n[0], t_eval=[0, .05, .1],
                          method='DOP853', rtol=1e-11, atol=1e-14, max_step=.01)
    np.testing.assert_allclose(np.array(saved)[:, 0], reference.y.T, rtol=1e-10, atol=1e-12)
    assert stats['numerical_checks_passed']


@pytest.mark.parametrize('kind', ['ttg', 'gaussian', 'selective'])
def test_initial_fields_branch_power_and_local_continuity(kind):
    initial = {'kind': kind, 'delta-temperature-K': 30, 'sigma-nm': .6, 'amplitude': .2}
    k, v, obs, options, x = setup(initial=initial)
    n = initial_occupation(obs, x, 8, 300, options['initial'])
    rhs = k.collision_batch(n)
    assert abs((rhs @ obs.energy)).max() < 1e-13
    if kind == 'selective':
        assert abs(obs.branch_sum(rhs*obs.energy)).max() > 1e-6
        assert np.max(obs.fields(n, rhs)['departure_from_local_BE']) > .01
    else:
        np.testing.assert_allclose(rhs, 0, atol=1e-12)
    saved = []
    stats = evolve_spatial(k, n, v, .5, [0., .025, .05, .1], .01, options, obs,
                          lambda i, t, n, f: saved.append(f))
    assert stats['numerical_checks_passed']
    assert stats['max_continuity_residual_fraction'] < 1e-12
    assert stats['max_energy_drift_fraction'] < 1e-12
    assert stats['max_collision_energy_residual_fraction'] < 1e-12
    assert stats['minimum_active_occupation'] > 0
    if kind == 'ttg':
        a = 2*temperature_harmonics(saved[0]['temperature_K'], x, 8)[0].real
        np.testing.assert_allclose(a, 30., atol=1e-10)


def test_periodic_pulse_and_energy_temperature_inversion():
    x = np.linspace(-1, 1, 100)
    np.testing.assert_allclose(periodic_gaussian(x, 1, .1, .12), periodic_gaussian(x+1, 1, .1, .12), atol=1e-14)
    k, v, obs, _, _ = setup()
    t = np.array([1., 100., 300., 3000.])
    np.testing.assert_allclose(obs.temperature(obs.bose(t)), t, rtol=1e-12)


@pytest.mark.parametrize('kind', ['selective', 'ttg', 'gaussian'])
def test_split_time_step_convergence_at_fixed_spatial_grid(kind):
    k, v, obs, options, x = setup(initial={'kind': kind, 'sigma-nm': .6, 'amplitude': .8, 'delta-temperature-K': 30})
    initial = initial_occupation(obs, x, 8, 300, options['initial'])
    results = []
    for dt in (.04, .02, .01, .00125):
        saved = []
        evolve_spatial(k, initial, v, .5, [0., .2], dt, options, obs,
                       lambda i, t, n, f: saved.append(n.copy()))
        results.append(saved[-1])
    errors = [np.linalg.norm(n-results[-1]) for n in results[:-1]]
    # Upwind time integration is first order even with Strang collision order.
    assert errors[0]/errors[1] > 1.7
    assert errors[1]/errors[2] > 1.7


@pytest.mark.parametrize('kind', ['ttg', 'gaussian'])
def test_thermal_profile_spatial_convergence_and_analytic_ballistics(kind):
    nonlinear = []
    ballistic_errors = []
    initial_options = {'kind': kind, 'delta-temperature-K': 30, 'sigma-nm': .6}
    for nx in (16, 32, 64, 128):
        k, v, obs, options, x = setup(nx=nx, initial=initial_options)
        n = initial_occupation(obs, x, 8, 300, options['initial'])
        for enabled in (True, False):
            options['collisions'] = enabled
            saved = []
            stats = evolve_spatial(k, n, v, 8/nx, [0., .2], .004, options, obs,
                                  lambda i, t, n, f: saved.append(n.copy()))
            assert stats['numerical_checks_passed']
            if enabled:
                nonlinear.append(saved[-1])
            elif nx < 128:
                exact = np.column_stack([initial_occupation(obs, x-v[j]*.2, 8, 300, options['initial'])[:, j]
                                         for j in range(6)])
                ballistic_errors.append(np.linalg.norm(saved[-1]-exact)/np.sqrt(nx))
    errors = []
    for n in nonlinear[:-1]:
        averaged = nonlinear[-1].reshape(len(n), -1, 6).mean(axis=1)
        errors.append(np.linalg.norm(n-averaged)/np.sqrt(len(n)))
    assert errors[0] > errors[1] > errors[2]
    assert ballistic_errors[0] > ballistic_errors[1] > ballistic_errors[2]


@pytest.mark.parametrize('change', [
    {'cells': True}, {'cfl': 1.1}, {'boundary': 'open'}, {'direction': [0, 0, 0]},
    {'initial': {'kind': 'ttg', 'delta-temperature-K': 300}},
    {'initial': {'kind': 'gaussian'}}, {'initial': {'kind': 'selective', 'branches': [1, 1], 'sigma-nm': 1}},
    {'snapshot-times-ps': [0, .2]}, {'unknown': 2}, {'length-nm': -1},
])
def test_invalid_spatial_inputs(change):
    cfg = SimpleNamespace(tdbte_spatial={'length-nm': 8, 'cells': 16, **change},
        tdbte_temperature=300., tdbte_excitation=.01, tdbte_branches=[2],
        tdbte_duration_ps=.1, tdbte_mesh=[3, 3, 3])
    with pytest.raises(ValueError):
        spatial_options(cfg)


def test_yaml_roundtrip_and_legacy_route(tmp_path):
    from nepkappa.config import parse_workflow_args
    input = tmp_path/'input.yaml'
    input.write_text('''tdbte:
  force-constants: fc
  mesh: [3, 3, 3]
  spatial:
    length-nm: 100
    cells: 16
    direction: [1, 1, 0]
    initial: {kind: ttg, delta-temperature-K: 10}
output: {result-dir: output}
''')
    config = parse_workflow_args(input, command='tdbte')
    assert config.sections.tdbte.spatial['cells'] == 16
    assert config.tdbte_branches is None
    old = parse_workflow_args('examples/tdbte.yaml', command='run')
    assert old.tdbte_spatial is None


def test_native_velocities_match_finite_differences_in_bz_gauge(tmp_path):
    """Real non-Bravais FC2: reciprocal shifts change eigenvector phases."""
    import phono3py
    from phonopy.structure.atoms import PhonopyAtoms
    from phono3py.file_IO import write_fc2_to_hdf5, write_fc3_to_hdf5
    from nepkappa.tdbte_builder import _prepare_modes, prepare_spatial_modes, source_hashes
    ph = phono3py.Phono3py(PhonopyAtoms(symbols=['Si', 'C'], cell=np.eye(3)*3,
        scaled_positions=[[0, 0, 0], [.2, .3, .4]]), supercell_matrix=[3, 3, 3],
        primitive_matrix=np.eye(3), is_symmetry=False, is_mesh_symmetry=False)
    positions = ph.supercell.scaled_positions
    symbols = np.asarray(ph.supercell.symbols)
    fc2 = np.zeros((len(positions), len(positions), 3, 3))
    spring = np.diag([1., 1.7, 2.9])
    for i in np.flatnonzero(symbols == 'Si'):
        for shift in ([0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]):
            target = positions[i] + (np.array([.2, .3, .4])+shift)/3
            delta = positions-target
            delta -= np.rint(delta)
            j = int(np.argmin(np.linalg.norm(delta, axis=1)))
            assert symbols[j] == 'C'
            fc2[i, i] += spring; fc2[j, j] += spring
            fc2[i, j] -= spring; fc2[j, i] -= spring
    ph.fc2 = fc2
    ph.fc3 = np.zeros((2, len(positions), len(positions), 3, 3, 3))
    source = tmp_path/'source'; source.mkdir()
    ph.save(source/'phono3py_disp.yaml')
    write_fc2_to_hdf5(fc2, filename=source/'fc2.hdf5')
    write_fc3_to_hdf5(ph.fc3, filename=source/'fc3.hdf5', p2s_map=ph.primitive.p2s_map)
    ph.mesh_numbers = [3, 3, 3]; ph.init_phph_interaction()
    pp = ph.phph_interaction
    addresses, frequency, zero, _ = _prepare_modes(pp, np.array([3, 3, 3]))
    metadata = {'mesh': [3, 3, 3], 'source': str(source), 'input_sha256': source_hashes(source),
        'uniform_mode_weight': 1/27,
        'phono3py_version': phono3py.__version__, 'grid_addresses': addresses.tolist(),
        'mode_branch_indices': np.tile(np.arange(6), 27).tolist()}
    kernel = SimpleNamespace(frequency=frequency.ravel(), zero=zero.ravel())
    modes = prepare_spatial_modes(metadata, kernel, [1, 0, 0])
    expected = np.zeros((27, 6, 3))
    eps = 1e-6
    for iq, q in enumerate(modes['native_q_fractional']):
        if zero[iq].any():
            continue
        assert np.min(np.diff(frequency[iq])) > 1e-5
        for axis in range(3):
            dq = ph.primitive.cell @ (np.eye(3)[axis]*eps)
            values = []
            for point in (q+dq, q-dq):
                pp.dynamical_matrix.run(point)
                eig = np.linalg.eigvalsh(pp.dynamical_matrix.dynamical_matrix)
                values.append(np.sqrt(eig)*pp.frequency_factor_to_THz)
            expected[iq, :, axis] = .1*(values[0]-values[1])/(2*eps)
    # Gamma optical derivatives vanish here; fixed translations have v=0.
    np.testing.assert_allclose(modes['velocity_nm_ps'].reshape(expected.shape), expected, atol=1e-7, rtol=1e-6)
    assert np.any(modes['native_q_fractional'] < 0)
    metadata['grid_addresses'][1], metadata['grid_addresses'][2] = metadata['grid_addresses'][2], metadata['grid_addresses'][1]
    with pytest.raises(ValueError, match='ordering'):
        prepare_spatial_modes(metadata, kernel, [1, 0, 0])


def test_cli_hdf5_plot_and_report(tmp_path, monkeypatch):
    import json
    import h5py
    import yaml
    from nepkappa.config import parse_workflow_args
    from nepkappa.application import execute_workflow_command
    from nepkappa.report import generate_report
    k = toy_kernel()
    path = tmp_path/'kernel.npz'
    np.savez(path, frequency_THz=k.frequency, zero_modes=k.zero, parent=k.parent,
        daughter1=k.daughter1, daughter2=k.daughter2, weights=k.weights, coefficients=k.coefficients)
    path.with_suffix('.json').write_text(json.dumps({'model': 'algebraic test',
        'mode_branch_indices': [0, 1, 2]*2, 'coefficient_time_unit': 'ps',
        'uniform_mode_weight': .5, 'mesh': [3, 3, 3], 'input_sha256': {}}))
    monkeypatch.setattr('nepkappa.tdbte_builder.prepare_spatial_modes', lambda *a, **kw: {
        'velocity_nm_ps': np.array([[1., 0, 0]]*3 + [[-1., 0, 0]]*3),
        'grid_addresses': np.zeros((2, 3), int), 'native_q_fractional': np.zeros((2, 3)),
        'primitive_volume_nm3': 1., 'primitive_cell_nm': np.eye(3),
        'eigenvectors': np.zeros((2, 3, 3)), 'velocity_basis': 'test fixture',
        'max_degenerate_offdiagonal_velocity_nm_ps': 0.})
    result = tmp_path/'result'
    input = tmp_path/'input.yaml'
    input.write_text(yaml.safe_dump({'tdbte': {'kernel': str(path), 'duration-ps': .02,
        'max-step-ps': .01, 'samples': 3, 'spatial': {'cells': 8, 'length-nm': 8}},
        'output': {'result-dir': str(result)}}))
    config = parse_workflow_args(input, command='tdbte')
    execute_workflow_command('tdbte', config)
    with h5py.File(result/'tdbte/spatial.hdf5') as f:
        assert f.attrs['complete']
        assert f['occupation'].shape == (3, 8, 6)
        assert f['branch_collision_power_eV_nm3_ps'].shape == (3, 8, 3)
        assert f['native_q_fractional'].shape == (2, 3)
    plot = tmp_path/'plot.yaml'
    plot.write_text(yaml.safe_dump({'plot': {'layout': 'both', 'dpi': 50}, 'output': {'result-dir': str(result)}}))
    execute_workflow_command('plot', parse_workflow_args(plot, command='plot'))
    assert (result/'plots/tdbte_spatial.png').is_file()
    _, markdown = generate_report(result)
    assert 'Branch continuity residual' in markdown.read_text()
    with pytest.raises(FileExistsError):
        execute_workflow_command('tdbte', config)


def test_independent_fourier_reference_matches_bose_polynomial_derivative():
    import importlib.util
    from pathlib import Path
    from nepkappa.tdbte import H_THz, KB
    path = Path(__file__).resolve().parents[2]/'tools/validate_spatial_tdbte.py'
    if not path.is_file():
        pytest.skip('Optional local spatial audit tool is unavailable')
    spec = importlib.util.spec_from_file_location('spatial_benchmark', path)
    benchmark = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(benchmark)
    k = toy_kernel()
    n = k.bose(300.)
    J = np.empty((6, 6))
    for i in range(6):
        delta = np.zeros(6); delta[i] = n[i]*1e-5
        J[:, i] = (k.collision(n+delta)[0]-k.collision(n-delta)[0])/(2*delta[i])
    odd = np.vstack([np.eye(3), -np.eye(3)])/np.sqrt(2)
    projected = -odd.T @ J @ odd
    np.testing.assert_allclose(projected, np.diag(np.diag(projected)), atol=1e-10)
    # Independently differentiated nonlinear operator, not fitted to a decay.
    capacity = (H_THz*k.frequency[:3])**2*n[:3]*(1+n[:3])/(KB*300**2)
    alpha = np.sum(capacity/np.diag(projected))/np.sum(capacity)
    np.testing.assert_allclose(benchmark.fourier_diffusivity(), alpha, rtol=1e-8)


def test_complete_synthetic_benchmark_report(tmp_path):
    from nepkappa.tdbte_accumulate import accumulate
    if accumulate is None:
        pytest.skip('Full convergence benchmark runs in the Numba job; reference parity is tested separately')
    import importlib.util
    from pathlib import Path
    path = Path(__file__).resolve().parents[2]/'tools/validate_spatial_tdbte.py'
    if not path.is_file():
        pytest.skip('Optional local spatial audit tool is unavailable')
    spec = importlib.util.spec_from_file_location('spatial_benchmark_report', path)
    benchmark = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(benchmark)
    report = benchmark.synthetic_benchmarks(tmp_path)
    assert report['checks_passed']
    assert {r['kind'] for r in report['thermal_profile_convergence']} == {'ttg', 'gaussian'}


def test_interpolated_batch_and_invalid_occupations():
    k = ShellKernel(np.array([5., 1., 2., 3., 4., 0.]), np.array([0]),
        np.array([[1, 2, 1, 2]]), np.array([[3, 4, 3, 4]]),
        np.full((1, 4), .25), np.array([.2]), np.array([False]*5+[True]))
    n = np.vstack([k.bose(100), k.bose(300)*np.array([1.5, .8, 1.2, 2., .7, 1])])
    np.testing.assert_allclose(k.collision_batch(n), [k.collision(row)[0] for row in n], rtol=1e-12, atol=1e-13)
    n[1, 2] = 0
    with pytest.raises(ValueError, match='positive'):
        k.collision_batch(n)
