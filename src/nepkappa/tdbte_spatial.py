"""Periodic spatial nonlinear phonon TD-BTE: Strang collisions/upwind streaming."""
import json
from pathlib import Path
import time

import h5py
import numpy as np
from scipy.integrate import solve_ivp

from .tdbte_observables import Observables, temperature_harmonics
from .tdbte_streaming import stream, divergence


def _mapping(value, allowed, name):
    if not isinstance(value, dict):
        raise ValueError(f'{name} must be a mapping')
    normalized = {str(key).replace('-', '_'): setting for key, setting in value.items()}
    if len(normalized) != len(value):
        raise ValueError(f'Duplicate normalized key in {name}')
    unknown = normalized.keys() - allowed
    if unknown:
        raise ValueError(f'Unknown {name} option: {sorted(unknown)[0]}')
    return normalized


def spatial_options(config):
    s = _mapping(config.tdbte_spatial, {
        'length_nm', 'cells', 'direction', 'boundary', 'cfl', 'batch_size',
        'collisions', 'rtol', 'atol', 'snapshot_times_ps', 'initial', 'phonons'}, 'tdbte.spatial')
    for key in ('length_nm', 'cells'):
        if key not in s:
            raise ValueError(f'tdbte.spatial.{key.replace("_", "-")} is required')
    s = {'direction': [1., 0., 0.], 'boundary': 'periodic', 'cfl': .8,
         'batch_size': 16, 'collisions': True, 'rtol': 1e-11, 'atol': 1e-14, **s}
    for key in ('cells', 'batch_size'):
        if isinstance(s[key], bool) or not isinstance(s[key], int) or s[key] < 1:
            raise ValueError(f'tdbte.spatial.{key} must be a positive integer')
    for key in ('length_nm', 'cfl', 'rtol', 'atol'):
        if isinstance(s[key], bool) or not isinstance(s[key], (int, float)) or not np.isfinite(s[key]) or s[key] <= 0:
            raise ValueError(f'tdbte.spatial.{key} must be finite and positive')
    if s['cfl'] > 1 or s['rtol'] >= 1 or s['atol'] >= 1:
        raise ValueError('Spatial CFL must be <=1 and collision tolerances <1')
    if s['boundary'] != 'periodic' or not isinstance(s['collisions'], bool):
        raise ValueError('Spatial TD-BTE requires periodic boundaries and boolean collisions')
    direction = np.asarray(s['direction'], dtype=float)
    if direction.shape != (3,) or not np.isfinite(direction).all() or np.linalg.norm(direction) == 0:
        raise ValueError('tdbte.spatial.direction must be a nonzero Cartesian three-vector')
    s['direction'] = (direction / np.linalg.norm(direction)).tolist()
    initial = _mapping(s.get('initial', {}), {'kind', 'delta_temperature_K', 'sigma_nm', 'center_nm',
                                             'profile', 'branches', 'amplitude'}, 'tdbte.spatial.initial')
    initial = {'kind': 'ttg', 'delta_temperature_K': 1., 'profile': 'gaussian',
               'center_nm': s['length_nm']/2, 'amplitude': config.tdbte_excitation, **initial}
    if initial['kind'] not in ('ttg', 'gaussian', 'uniform', 'selective'):
        raise ValueError('Spatial initial.kind must be ttg, gaussian, uniform or selective')
    if initial['kind'] == 'ttg' and s['cells'] < 4:
        raise ValueError('TTG requires at least four spatial cells')
    for key in ('delta_temperature_K', 'center_nm', 'amplitude'):
        if isinstance(initial[key], bool) or not isinstance(initial[key], (int, float)) or not np.isfinite(initial[key]):
            raise ValueError(f'Spatial initial.{key} must be finite')
    if initial['profile'] not in ('gaussian', 'uniform'):
        raise ValueError('Spatial selective profile must be gaussian or uniform')
    if initial['kind'] == 'gaussian' or (initial['kind'] == 'selective' and initial['profile'] == 'gaussian'):
        sigma = initial.get('sigma_nm')
        if isinstance(sigma, bool) or not isinstance(sigma, (int, float)) or not np.isfinite(sigma) or sigma <= 0:
            raise ValueError('Gaussian initial profiles require positive initial.sigma-nm')
    if initial['kind'] == 'ttg' and abs(initial['delta_temperature_K']) >= config.tdbte_temperature:
        raise ValueError('TTG temperature must remain positive throughout the domain')
    if initial['kind'] in ('gaussian', 'uniform') and config.tdbte_temperature + min(initial['delta_temperature_K'], 0) <= 0:
        raise ValueError('Spatial initial temperatures must be positive')
    if initial['kind'] == 'selective':
        initial['branches'] = initial.get('branches', config.tdbte_branches)
        branches = initial['branches']
        if (not isinstance(branches, list) or not branches or any(isinstance(b, bool) or not isinstance(b, int) or b < 0 for b in branches)
                or len(set(branches)) != len(branches) or initial['amplitude'] <= 0):
            raise ValueError('Selective excitation requires unique nonnegative branches and positive amplitude')
    s['initial'] = initial
    if 'phonons' in s and (not isinstance(s['phonons'], str) or not s['phonons'].strip()):
        raise ValueError('tdbte.spatial.phonons must be a nonempty matching force-constant directory path')
    times = s.get('snapshot_times_ps')
    if times is not None:
        t = np.asarray(times, float)
        if (t.ndim != 1 or len(t) < 2 or not np.isfinite(t).all() or np.any(np.diff(t) <= 0)
                or t[0] != 0 or t[-1] != config.tdbte_duration_ps):
            raise ValueError('Spatial snapshot-times-ps must increase from 0 to duration-ps inclusive')
    if getattr(config, 'tdbte_mesh', None) is not None and np.any(np.asarray(config.tdbte_mesh) < 2):
        raise ValueError('Spatial TD-BTE requires a 3D phonon q mesh')
    return s


def periodic_gaussian(x, length, center, sigma):
    """Wrapped Gaussian image sum, normalized to one at its center.

    Unlike a minimum-distance truncated Gaussian, derivatives match at the
    periodic boundary. For sigma << length this is the requested Gaussian.
    """
    distance = (np.asarray(x) - center + length/2) % length - length/2
    images = np.arange(-int(np.ceil(9*sigma/length))-1, int(np.ceil(9*sigma/length))+2)
    profile = np.exp(-.5 * ((distance[..., None] + images*length) / sigma)**2).sum(axis=-1)
    normalization = np.exp(-.5 * (images*length / sigma)**2).sum()
    return profile / normalization


def initial_occupation(obs, x, length, temperature, options):
    kind = options['kind']
    if kind == 'ttg':
        return obs.bose(temperature + options['delta_temperature_K'] * np.cos(2*np.pi*x/length))
    if kind == 'uniform':
        return obs.bose(np.full(len(x), temperature + options['delta_temperature_K']))
    profile = (np.ones(len(x)) if kind == 'selective' and options['profile'] == 'uniform'
               else periodic_gaussian(x, length, options['center_nm'], options['sigma_nm']))
    if kind == 'gaussian':
        return obs.bose(temperature + options['delta_temperature_K'] * profile)
    selected = np.asarray(options['branches'])
    if not set(selected) <= set(obs.labels):
        raise ValueError('Selective excited branch is absent from the collision kernel')
    n = obs.bose(np.full(len(x), temperature))
    n[:, np.isin(obs.branches, selected)] *= 1 + options['amplitude'] * profile[:, None]
    return n


def collision_step(kernel, state, duration, options, stats):
    """DOP853 on the original nonlinear RHS, with explicit positivity retries."""
    if not options['collisions']:
        return state.copy()
    shape = state.shape

    def rhs(t, flat):
        stats['collision_rhs_evaluations'] += 1
        return kernel.collision_batch(flat.reshape(shape), options['batch_size']).ravel()

    maximum = duration
    for attempt in range(9):
        try:
            sol = solve_ivp(rhs, (0., duration), state.ravel(), method='DOP853',
                            t_eval=[duration], rtol=options['rtol'], atol=options['atol'],
                            max_step=maximum, first_step=maximum)
            if not sol.success:
                raise RuntimeError(f'Spatial collision integration failed: {sol.message}')
            result = sol.y[:, -1].reshape(shape)
            if not np.isfinite(result).all() or np.any(result[:, ~kernel.zero] <= 0):
                raise ValueError('Active occupations must be finite and positive')
            return result
        except ValueError as exc:
            if 'occupations must be finite and positive' not in str(exc) or attempt == 8:
                raise
            maximum /= 2
            stats['collision_positivity_retries'] += 1
    raise RuntimeError('Spatial collision positivity retries exhausted')  # pragma: no cover


def evolve_spatial(kernel, initial, velocity, dx, times, maximum_step, options, obs, save):
    """Integrate without retaining internal states; save only requested times."""
    n = np.array(initial, dtype=float, copy=True)
    if n.ndim != 2 or not np.isfinite(n).all() or np.any(n[:, ~kernel.zero] <= 0):
        raise ValueError('Invalid spatial initial occupations')
    speed = np.max(abs(velocity), initial=0.)
    limit = min(maximum_step, options['cfl'] * dx / speed) if speed else maximum_step
    if not np.isfinite(limit) or limit <= 0:
        raise ValueError('Spatial step limit must be positive and finite')
    stats = {'steps': 0, 'collision_rhs_evaluations': 0, 'collision_positivity_retries': 0,
             'max_cfl': 0., 'max_energy_drift_fraction': 0., 'max_continuity_residual_fraction': 0.,
             'max_collision_energy_residual_fraction': 0., 'max_entropy_decrease': 0.,
             'minimum_active_occupation': float(n[:, ~kernel.zero].min())}
    initial_energy = float((n @ obs.energy).sum())
    previous_energy = obs.branch_energy(n)
    branch_shape = previous_energy.shape
    integrated_flux = np.zeros(branch_shape)
    integrated_collision = np.zeros(branch_shape)
    previous_entropy = None
    t = 0.
    progress = time.monotonic()
    for index, target in enumerate(times):
        while t < target:
            dt = min(limit, target - t)
            if t + dt == t:
                raise RuntimeError('Spatial timestep underflow')
            before = n
            first = collision_step(kernel, before, dt/2, options, stats)
            streamed, face = stream(first, velocity, dt, dx)
            n = collision_step(kernel, streamed, dt/2, options, stats)
            integrated_collision += obs.branch_energy(first-before) + obs.branch_energy(n-streamed)
            integrated_flux += dt * obs.branch_sum(face * obs.energy)
            drift = abs(float((n @ obs.energy).sum()) - initial_energy) / initial_energy
            stats['max_energy_drift_fraction'] = max(stats['max_energy_drift_fraction'], drift)
            stats['max_cfl'] = max(stats['max_cfl'], dt * speed / dx)
            stats['minimum_active_occupation'] = min(stats['minimum_active_occupation'], float(n[:, ~kernel.zero].min()))
            stats['steps'] += 1
            t = min(target, t + dt)
            if time.monotonic() - progress >= 30:
                print(f'Spatial TD-BTE: t={t:.6g}/{times[-1]:.6g} ps; {stats["steps"]} split steps', flush=True)
                progress = time.monotonic()
        collision = (kernel.collision_batch(n, options['batch_size']) if options['collisions'] else np.zeros_like(n))
        fields = obs.fields(n, collision)
        residual = fields['branch_energy_eV_nm3'] - previous_energy + divergence(integrated_flux, dx) - integrated_collision
        scale = max(float(abs(fields['energy_eV_nm3']).max()), 1e-300)
        stats['max_continuity_residual_fraction'] = max(stats['max_continuity_residual_fraction'], float(abs(residual).max())/scale)
        stats['max_collision_energy_residual_fraction'] = max(stats['max_collision_energy_residual_fraction'], float(abs(integrated_collision.sum(axis=-1)).max())/scale)
        entropy = float(fields['entropy_kB_nm3'].sum())
        if previous_entropy is not None:
            stats['max_entropy_decrease'] = max(stats['max_entropy_decrease'], (previous_entropy-entropy)/max(abs(previous_entropy), 1e-300))
        previous_entropy = entropy
        fields.update(branch_face_flux_integral_eV_nm2=integrated_flux.copy(),
                      branch_collision_energy_integral_eV_nm3=integrated_collision.copy(),
                      branch_continuity_residual_eV_nm3=residual)
        save(index, target, n, fields)
        previous_energy = fields['branch_energy_eV_nm3'].copy()
        integrated_collision.fill(0)
        integrated_flux.fill(0)
    stats['numerical_checks_passed'] = bool(stats['max_energy_drift_fraction'] < 1e-8
        and stats['max_continuity_residual_fraction'] < 1e-8
        and stats['max_collision_energy_residual_fraction'] < 1e-8
        and stats['minimum_active_occupation'] > 0 and stats['max_entropy_decrease'] < 1e-10)
    return stats


def run_spatial_tdbte(config):
    from .tdbte import load_tdbte_kernel, validate_options, H_THz, KB
    from .tdbte_builder import prepare_spatial_modes
    from .tdbte_storage import atomic_json, sha256
    validate_options(config)
    options = spatial_options(config)
    out = Path(config.result_dir) / 'tdbte'
    if out.exists():
        raise FileExistsError(f'Refusing to overwrite TD-BTE results: {out}')
    kernel, metadata, source, metadata_path, branches, weight = load_tdbte_kernel(config)
    modes = prepare_spatial_modes(metadata, kernel, options['direction'], options.get('phonons'))
    velocity = modes['velocity_nm_ps'] @ options['direction']
    obs = Observables(kernel, branches, weight, modes['primitive_volume_nm3'], velocity)
    nx, length = options['cells'], options['length_nm']
    dx = length / nx
    x = (np.arange(nx) + .5) * dx
    times = np.asarray(options.get('snapshot_times_ps', np.linspace(0., config.tdbte_duration_ps, config.tdbte_samples)), float)
    initial = initial_occupation(obs, x, length, config.tdbte_temperature, options['initial'])
    out.mkdir(parents=True, exist_ok=False)
    report = {'spatial': True, 'numerical_checks_passed': False, 'validated_for_physical_dynamics': False,
              'kernel_path': str(source.resolve()), 'kernel_sha256': sha256(source),
              'metadata_sha256': sha256(metadata_path), 'mesh': metadata['mesh'],
              'input_sha256': metadata['input_sha256'], 'spatial_options': options,
              'collision_backend': getattr(kernel, 'backend', 'numpy'),
              'velocity_basis': modes['velocity_basis'],
              'max_degenerate_offdiagonal_velocity_nm_ps': modes['max_degenerate_offdiagonal_velocity_nm_ps'],
              'integration': 'Strang: half nonlinear DOP853 collision, full upwind streaming, half collision',
              'maximum_split_step_ps': config.tdbte_max_step_ps,
              'limitations': 'Fixed-frequency population model; first-order spatial numerical diffusion; absolute collision rates and physical Fourier limit require independent convergence; degenerate coherences neglected'}
    atomic_json(out/'audit.json', report)
    try:
        with h5py.File(out/'spatial.hdf5', 'x') as handle:
            handle.attrs.update(format='nepkappa-spatial-tdbte-v1', complete=False,
                                metadata_json=json.dumps(report), uniform_mode_weight=weight,
                                primitive_volume_nm3=modes['primitive_volume_nm3'],
                                mode_order='q-major, native frequency-ordered branch minor',
                                energy_zero='thermal phonon energy; zero-point energy omitted')
            handle.attrs['heat_flux_eV_nm2_ps_to_W_m2'] = 1.602176634e11
            handle.attrs['energy_eV_nm3_to_J_m3'] = 1.602176634e8
            handle.attrs['h_eV_per_THz'] = H_THz
            handle.attrs['kB_eV_per_K'] = KB
            handle.attrs['interval_integrals'] = 'Since previous saved snapshot; zero at t=0; face i is right face of cell i'
            for name, value in {'time_ps': times, 'x_nm': x, 'face_x_nm': (np.arange(nx)+1)*dx,
                'frequency_THz': kernel.frequency, 'zero_modes': kernel.zero, 'mode_branch_indices': branches,
                'branch_indices': obs.labels, 'velocity_nm_ps': modes['velocity_nm_ps'],
                'temperature_harmonic_indices': np.arange(1, 5),
                'temperature_harmonic_valid': np.arange(1, 5) <= (nx-1)//2,
                'velocity_x_nm_ps': velocity, 'grid_addresses': modes['grid_addresses'],
                'native_q_fractional': modes['native_q_fractional'],
                'primitive_cell_nm': modes['primitive_cell_nm'], 'native_eigenvectors': modes['eigenvectors']}.items():
                handle.create_dataset(name, data=value)
            occupation = handle.create_dataset('occupation', shape=(len(times), nx, len(kernel.frequency)),
                dtype='f8', chunks=(1, min(nx, options['batch_size']), min(len(kernel.frequency), 1024)),
                compression='gzip', shuffle=True)
            datasets = {}
            def save(index, t, n, fields):
                occupation[index] = n
                harmonics = temperature_harmonics(fields['temperature_K'], x, length)
                fields.update(temperature_harmonics_real_K=harmonics.real,
                              temperature_harmonics_imag_K=harmonics.imag,
                              ttg_cosine_amplitude_K=np.asarray(2*harmonics[0].real))
                for name, value in fields.items():
                    if name not in datasets:
                        datasets[name] = handle.create_dataset(name, shape=(len(times),)+np.shape(value), dtype='f8',
                                                               chunks=True, compression='gzip', shuffle=True)
                    datasets[name][index] = value
                handle.attrs['saved_snapshots'] = index + 1
                handle.flush()
            stats = evolve_spatial(kernel, initial, velocity, dx, times, config.tdbte_max_step_ps, options, obs, save)
            report.update(stats)
            handle.attrs['complete'] = stats['numerical_checks_passed']
            handle.attrs['audit_json'] = json.dumps(report)
        report['spatial_hdf5_sha256'] = sha256(out/'spatial.hdf5')
        atomic_json(out/'audit.json', report)
        if not report['numerical_checks_passed']:
            raise RuntimeError('Spatial TD-BTE numerical audit failed; diagnostic outputs preserved')
    except Exception as exc:
        report['numerical_checks_passed'] = False
        report['error'] = str(exc)
        atomic_json(out/'audit.json', report)
        raise
    print('Spatial TD-BTE: numerical checks passed; absolute physical rates remain unvalidated.', flush=True)
    return report
