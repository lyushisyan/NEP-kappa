"""Experimental, spatially homogeneous three-phonon occupation dynamics.

Consumes an explicit energy-shell artifact, not arbitrary FC2/FC3. Passing the
audits below is mathematical consistency, NOT physical rate validation.
"""
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

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
    if not config.tdbte_experimental:
        raise ValueError('tdbte requires experimental: true; physical rates are not validated')
    if not config.tdbte_kernel:
        raise ValueError('tdbte.kernel is required (prebuilt energy-shell NPZ)')
    for key in ('temperature','duration_ps','max_step_ps','excitation'):
        value = getattr(config,'tdbte_'+key)
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f'tdbte.{key} must be finite and positive')
    if not 2 <= config.tdbte_samples <= 100001:
        raise ValueError('tdbte.samples must be between 2 and 100001')
    if not config.tdbte_branches or min(config.tdbte_branches) < 0 or len(set(config.tdbte_branches)) != len(config.tdbte_branches):
        raise ValueError('tdbte.branches requires unique zero-based branch indices')


def run_tdbte(config):
    """Audit and evolve a prebuilt artifact, with an independent equilibrium control."""
    validate_options(config)
    source = Path(config.tdbte_kernel)
    # Sidecar provides mode ordering; never infer a six-branch SiC convention.
    metadata = json.loads(source.with_suffix('.json').read_text())
    branch_ids = np.asarray(metadata.get('mode_branch_indices', []))
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
    out = Path(config.result_dir)/'tdbte'
    out.mkdir(parents=True, exist_ok=False)
    times = np.linspace(0,config.tdbte_duration_ps,config.tdbte_samples)
    solutions = []
    for state in (n0,initial):
        try:
            sol = solve_ivp(lambda t,n: kernel.collision(n)[0], (0,times[-1]),state,
                            t_eval=times,method='DOP853',rtol=1e-11,atol=1e-14,
                            max_step=config.tdbte_max_step_ps)
        except Exception as exc:
            (out/'audit.json').write_text(json.dumps({
                'experimental':True,'validated_for_physical_dynamics':False,
                'numerical_checks_passed':False,'error':str(exc)},indent=2)+'\n')
            raise
        if not sol.success or not np.isfinite(sol.y).all() or np.any(sol.y[~kernel.zero] <= 0):
            (out/'audit.json').write_text(json.dumps({
                'experimental':True,'validated_for_physical_dynamics':False,
                'numerical_checks_passed':False,'error':sol.message},indent=2)+'\n')
            raise RuntimeError('Time integration failed positivity/finite checks')
        solutions.append(sol.y)
    control, excited = solutions
    drift = float(np.max(abs(energy@(excited-initial[:,None])))/injected)
    equilibrium = float(np.max(energy@abs(control-n0[:,None]))/injected)
    active = excited[~kernel.zero]
    entropy = ((1+active)*np.log1p(active)-active*np.log(active)).sum(axis=0)
    entropy_drop = float(max(0.,-np.diff(entropy).min()))
    passed = drift < 1e-8 and equilibrium < 1e-8 and entropy_drop < 1e-10*max(1.,abs(entropy[0]))
    report = {'experimental':True,'validated_for_physical_dynamics':False,
              'numerical_checks_passed':bool(passed),'model':metadata['model'],
              'kernel_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
              'metadata_sha256':hashlib.sha256(source.with_suffix('.json').read_bytes()).hexdigest(),
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
    print('EXPERIMENTAL: numerical checks passed; physical rate validation remains incomplete.')
    return report
