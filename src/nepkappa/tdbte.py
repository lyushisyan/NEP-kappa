"""Spatially homogeneous three-phonon occupation dynamics.

Builds an energy-shell kernel from matching FC2/FC3 or reuses an audited
artifact. Passing numerical audits is NOT physical rate validation.
"""
from dataclasses import dataclass
import json
from pathlib import Path
import time

import numpy as np
from scipy.constants import physical_constants
from scipy.integrate import solve_ivp

H_THz = physical_constants['Planck constant in eV/Hz'][0] * 1e12
KB = physical_constants['Boltzmann constant in eV/K'][0]


@dataclass(frozen=True)
class ShellKernel:
    frequency: np.ndarray
    parent: np.ndarray
    daughter1: np.ndarray
    daughter2: np.ndarray
    weights: np.ndarray
    coefficients: np.ndarray
    zero: np.ndarray

    @classmethod
    def load(cls, path):
        """Load numeric arrays only; reject malformed or off-shell artifacts."""
        with np.load(path, allow_pickle=False) as data:
            kernel = cls(*(data[k].copy() for k in (
                'frequency_THz', 'parent', 'daughter1', 'daughter2', 'weights',
                'coefficients', 'zero_modes')))
        kernel.validate()
        return kernel

    def validate(self):
        f, p, w, c = self.frequency, self.parent, self.weights, self.coefficients
        if f.ndim != 1 or len(f) == 0 or p.ndim != 1 or len(p) == 0:
            raise ValueError('Empty or malformed mode/event arrays')
        if self.zero.dtype != np.bool_ or self.zero.shape != f.shape:
            raise ValueError('zero_modes must be a boolean mode mask')
        if w.shape != (len(p), 4) or c.shape != p.shape:
            raise ValueError('Invalid shell weight/coefficient shape')
        for a in (f, w, c):
            if not np.isfinite(a).all() or np.any(a < 0):
                raise ValueError('Frequencies, weights, coefficients must be finite and nonnegative')
        if np.any(f[self.zero] != 0) or np.any(f[~self.zero] <= 0):
            raise ValueError('Only marked rigid modes may have zero frequency')
        if not np.allclose(w.sum(axis=1), 1, rtol=0, atol=1e-12):
            raise ValueError('Barycentric weights must sum to one')
        for index, shape in [(p,p.shape),(self.daughter1,w.shape),(self.daughter2,w.shape)]:
            if index.shape != shape or index.dtype.kind not in 'iu' or np.any(index < 0) or np.any(index >= len(f)):
                raise ValueError('Invalid mode index')
        if np.any(self.zero[p]):
            raise ValueError('Rigid translations cannot be event parents')
        daughters = [(w*f[d]).sum(axis=1) for d in (self.daughter1,self.daughter2)]
        if any(np.any(d <= 0) for d in daughters):
            raise ValueError('Zero-energy shell daughter is singular')
        residual = f[p]-daughters[0]-daughters[1]
        if np.max(abs(residual)) > 1e-10 * max(1.,float(f.max())):
            raise ValueError('Kernel violates energy-shell identity')

    def bose(self, temperature):
        if not np.isfinite(temperature) or temperature <= 0:
            raise ValueError('Bose temperature must be finite and positive')
        n = np.zeros_like(self.frequency, dtype=float)
        x = self.frequency[~self.zero]*H_THz/(KB*temperature)
        n[~self.zero] = np.exp(-x)/(-np.expm1(-x))
        return n

    def collision(self, occupation):
        n = np.asarray(occupation, dtype=float)
        if n.shape != self.frequency.shape or not np.isfinite(n).all() or np.any(n[~self.zero] <= 0):
            raise ValueError('Active occupations must be finite and positive; no clipping is applied')
        eta = np.zeros_like(n)
        eta[~self.zero] = np.log1p(n[~self.zero])-np.log(n[~self.zero])
        a = n[self.parent]
        eb = (self.weights*eta[self.daughter1]).sum(axis=1)
        ec = (self.weights*eta[self.daughter2]).sum(axis=1)
        if np.any(eb <= 0) or np.any(ec <= 0):
            raise ValueError('Zero-entropy shell daughter is singular')
        # Log Bose occupations avoid exp(eta) overflow at low temperatures.
        lb = -eb-np.log(-np.expm1(-eb))
        lc = -ec-np.log(-np.expm1(-ec))
        affinity = -eta[self.parent]+eb+ec
        log_reverse = np.log1p(a)+lb+lc
        positive = self.coefficients > 0
        flux = np.zeros_like(self.coefficients, dtype=float)
        log_scale = log_reverse[positive]+np.maximum(affinity[positive],0)+np.log(self.coefficients[positive])
        flux[positive] = np.exp(log_scale)*np.sign(affinity[positive])*(-np.expm1(-abs(affinity[positive])))
        if not np.isfinite(flux).all():
            raise ValueError('Nonfinite event flux; unsupported numerical regime')
        rhs = np.zeros_like(n)
        np.add.at(rhs,self.parent,-flux)
        for d in (self.daughter1,self.daughter2):
            np.add.at(rhs,d.ravel(),(self.weights*flux[:,None]).ravel())
        rhs[self.zero] = 0
        if not np.isfinite(rhs).all():
            raise ValueError('Nonfinite accumulated collision derivative')
        return rhs, flux


def validate_options(config):
    source = getattr(config, 'tdbte_force_constants', None)
    kernel = getattr(config, 'tdbte_kernel', None)
    mesh = getattr(config, 'tdbte_mesh', None)
    if bool(source) == bool(kernel):
        raise ValueError('Specify exactly one of tdbte.force-constants or tdbte.kernel')
    if source:
        from .tdbte_builder import validate_mesh
        validate_mesh(mesh)
    elif mesh is not None:
        raise ValueError('tdbte.mesh is only used with force-constants; a reused kernel fixes the mesh')
    for key in ('temperature','duration_ps','max_step_ps','excitation'):
        value = getattr(config,'tdbte_'+key)
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f'tdbte.{key} must be finite and positive')
    if not 2 <= config.tdbte_samples <= 100001:
        raise ValueError('tdbte.samples must be between 2 and 100001')
    if not config.tdbte_branches or min(config.tdbte_branches) < 0 or len(set(config.tdbte_branches)) != len(config.tdbte_branches):
        raise ValueError('tdbte.branches requires unique zero-based branch indices')


def run_tdbte(config):
    """Build/load a kernel and evolve it with an independent equilibrium control."""
    from .tdbte_storage import ChunkedKernel, sha256

    validate_options(config)
    out = Path(config.result_dir)/'tdbte'
    if out.exists():
        raise FileExistsError(f'Refusing to overwrite TD-BTE results: {out}')
    if getattr(config, 'tdbte_force_constants', None):
        from .tdbte_builder import build_kernel
        source = build_kernel(config.tdbte_force_constants, config.tdbte_mesh,
                              Path(config.result_dir)/'tdbte-kernel',
                              excited_branches=config.tdbte_branches)
    else:
        source = Path(config.tdbte_kernel)
    if source.suffix.lower() not in ('.json', '.npz'):
        raise ValueError('tdbte.kernel must be an NPZ kernel or a JSON chunk manifest')
    # Sidecar provides mode ordering; never infer a six-branch SiC convention.
    metadata_path = source if source.suffix.lower() == '.json' else source.with_suffix('.json')
    metadata = json.loads(metadata_path.read_text())
    branch_ids = np.asarray(metadata.get('mode_branch_indices', []))
    if source.suffix.lower() == '.json':
        kernel = ChunkedKernel(source)
    else:
        kernel = ShellKernel.load(source)
    if branch_ids.shape != kernel.frequency.shape or branch_ids.dtype.kind not in 'iu' or np.any(branch_ids < 0):
        raise ValueError('Kernel JSON requires explicit mode_branch_indices for every mode')
    if not set(config.tdbte_branches) <= set(branch_ids.tolist()):
        raise ValueError('Excited branch index is absent from kernel')
    if metadata.get('coefficient_time_unit') != 'ps' or not metadata.get('model'):
        raise ValueError('Kernel JSON requires model and coefficient_time_unit: ps')
    weight = metadata.get('uniform_mode_weight')
    if isinstance(weight, bool) or not isinstance(weight, (int, float)) or not np.isfinite(weight) or weight <= 0:
        raise ValueError('Kernel JSON requires a positive uniform_mode_weight; reduced weighted grids are unsupported')
    n0 = kernel.bose(config.tdbte_temperature)
    initial = n0.copy()
    initial[np.isin(branch_ids,config.tdbte_branches)] *= 1+config.tdbte_excitation
    energy = kernel.frequency*H_THz*weight
    injected = float(energy@(initial-n0))
    if injected <= 0:
        raise ValueError('Excitation injects no energy')
    if not np.isfinite(initial).all() or np.any(n0[~kernel.zero] <= 0):
        raise ValueError('Initial occupations underflowed or overflowed; unsupported temperature/excitation')
    out.mkdir(parents=True, exist_ok=False)
    # Record the source before integration. The manifest transitively hashes all
    # chunk arrays; a single-file kernel retains its adjacent metadata hash.
    provenance = {'kernel_path': str(source.resolve()), 'kernel_sha256': sha256(source),
                  'metadata_sha256': sha256(metadata_path),
                  'mesh': metadata.get('mesh'), 'input_sha256': metadata.get('input_sha256'),
                  'kernel_chunks': len(metadata.get('chunks', [])),
                  'collision_backend': getattr(kernel, 'backend', 'numpy'),
                  'kernel_limitations': metadata.get('limitations'),
                  'boundary_assumptions': metadata.get('boundary_assumptions')}
    times = np.linspace(0,config.tdbte_duration_ps,config.tdbte_samples)
    solutions = []
    last_progress = time.monotonic()

    def derivative(t, n):
        nonlocal last_progress
        value = kernel.collision(n)[0]
        if time.monotonic() - last_progress >= 30:
            print(f'TD-BTE integration: t={t:.6g}/{times[-1]:.6g} ps', flush=True)
            last_progress = time.monotonic()
        return value

    for label, state in (('equilibrium control', n0), ('excited state', initial)):
        print(f'TD-BTE: integrating {label}', flush=True)
        try:
            sol = solve_ivp(derivative, (0,times[-1]),state,
                            t_eval=times,method='DOP853',rtol=1e-11,atol=1e-14,
                            max_step=config.tdbte_max_step_ps)
        except Exception as exc:
            (out/'audit.json').write_text(json.dumps({
                'validated_for_physical_dynamics':False,
                'numerical_checks_passed':False,'error':str(exc), **provenance},indent=2)+'\n')
            raise
        if not sol.success or not np.isfinite(sol.y).all() or np.any(sol.y[~kernel.zero] <= 0):
            (out/'audit.json').write_text(json.dumps({
                'validated_for_physical_dynamics':False,
                'numerical_checks_passed':False,'error':sol.message, **provenance},indent=2)+'\n')
            raise RuntimeError('Time integration failed positivity/finite checks')
        solutions.append(sol.y)
    control, excited = solutions
    drift = float(np.max(abs(energy@(excited-initial[:,None])))/injected)
    equilibrium = float(np.max(energy@abs(control-n0[:,None]))/injected)
    active = excited[~kernel.zero]
    entropy = ((1+active)*np.log1p(active)-active*np.log(active)).sum(axis=0)
    entropy_drop = float(max(0.,-np.diff(entropy).min()))
    passed = drift < 1e-8 and equilibrium < 1e-8 and entropy_drop < 1e-10*max(1.,abs(entropy[0]))
    report = {'validated_for_physical_dynamics':False,
              'numerical_checks_passed':bool(passed),'model':metadata['model'],
              **provenance,
              'max_energy_drift_fraction':drift,'max_equilibrium_redistribution_fraction':equilibrium,
              'max_entropy_decrease':entropy_drop,'temperature_K':config.tdbte_temperature,
              'excitation_relative_occupation':config.tdbte_excitation,
              'excited_branch_indices':config.tdbte_branches,
              'injected_energy_eV_per_primitive_cell':injected,
              'uniform_mode_weight':weight,
              'integrator':{'method':'DOP853','rtol':1e-11,'atol':1e-14,
                            'max_step_ps':config.tdbte_max_step_ps,'duration_ps':config.tdbte_duration_ps},
              'limitations':'Homogeneous fixed-frequency three-phonon populations only; no pump coupling, bath, spatial drift or coherences.'}
    (out/'audit.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savez_compressed(out/'trajectories.npz',time_ps=times,frequency_THz=kernel.frequency,
                        equilibrium_control=control,excited=excited,mode_branch_indices=branch_ids,
                        entropy=entropy,uniform_mode_weight=weight)
    branch_labels = np.unique(branch_ids)
    excess = energy[:,None]*(excited-n0[:,None])
    branch_energy = np.array([excess[branch_ids == s].sum(axis=0) for s in branch_labels])
    np.savetxt(out/'branch-energy.csv',np.column_stack([times,branch_energy.T]),
               delimiter=',',comments='',header='time_ps,'+','.join(
                   f'branch_{s}_excess_eV_per_cell' for s in branch_labels))
    if not passed:
        raise RuntimeError('tdbte numerical audit failed; diagnostic artifacts saved, not accepted results')
    print('TD-BTE: numerical checks passed; physical rate validation remains incomplete.')
    return report
