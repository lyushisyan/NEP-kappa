"""Energy-equivalent temperature and mode/branch observables in nm, ps, eV."""
import numpy as np
from .tdbte import H_THz, KB


class Observables:
    def __init__(self, kernel, branches, weight, volume_nm3, velocity):
        if not np.isfinite(volume_nm3) or volume_nm3 <= 0:
            raise ValueError('Primitive volume must be finite and positive')
        self.kernel = kernel
        self.branches = np.asarray(branches)
        self.labels = np.unique(branches)
        self.velocity = np.asarray(velocity)
        self.energy = H_THz * kernel.frequency * weight / volume_nm3
        self.weight_density = weight / volume_nm3

    def branch_sum(self, values):
        return np.stack([values[..., self.branches == s].sum(axis=-1) for s in self.labels], axis=-1)

    def branch_energy(self, occupation):
        return self.branch_sum(occupation * self.energy)

    def bose(self, temperatures):
        t = np.asarray(temperatures, float)
        if not np.isfinite(t).all() or np.any(t <= 0):
            raise ValueError('Local temperatures must be positive and finite')
        result = np.zeros(t.shape + self.kernel.frequency.shape)
        x = H_THz * self.kernel.frequency[~self.kernel.zero] / (KB * t[..., None])
        result[..., ~self.kernel.zero] = np.exp(-x) / -np.expm1(-x)
        return result

    def temperature(self, occupation):
        """Invert thermal energy, excluding fixed zero-energy rigid modes."""
        energy = np.asarray(occupation) @ self.energy
        lower = np.zeros_like(energy)
        upper = np.full_like(energy, 300.)
        while np.any(self.bose(upper) @ self.energy < energy):
            upper *= 2
            if not np.isfinite(upper).all():
                raise ValueError('Cannot bracket local energy-equivalent temperature')
        for _ in range(52):
            middle = (lower + upper) / 2
            below = self.bose(middle) @ self.energy < energy
            lower = np.where(below, middle, lower)
            upper = np.where(below, upper, middle)
        return (lower + upper) / 2

    def fields(self, occupation, collision):
        n = np.asarray(occupation)
        t = self.temperature(n)
        equilibrium = self.bose(t)
        active = n[..., ~self.kernel.zero]
        entropy = ((1+active)*np.log1p(active)-active*np.log(active)).sum(axis=-1) * self.weight_density
        energy = self.branch_energy(n)
        return {
            'temperature_K': t,
            'energy_eV_nm3': energy.sum(axis=-1),
            'heat_flux_eV_nm2_ps': (n * self.energy * self.velocity).sum(axis=-1),
            'branch_energy_eV_nm3': energy,
            'branch_heat_flux_eV_nm2_ps': self.branch_sum(n * self.energy * self.velocity),
            'branch_collision_power_eV_nm3_ps': self.branch_sum(collision * self.energy),
            'entropy_kB_nm3': entropy,
            'departure_from_local_BE': (abs(n-equilibrium) @ self.energy) / (n @ self.energy),
        }


def temperature_harmonics(temperature, x, length, count=4):
    """Complex Fourier amplitudes; 2*Re(A[0]) is the TTG cosine amplitude."""
    wave = np.arange(1, count+1)[:, None] * 2 * np.pi * np.asarray(x)[None, :] / length
    result = (np.exp(-1j * wave) @ (temperature - np.mean(temperature))) / len(x)
    result[np.arange(1, count+1) > (len(x)-1)//2] = 0
    return result
