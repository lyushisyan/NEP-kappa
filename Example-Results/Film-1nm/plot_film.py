import os

import h5py
import numpy as np
import matplotlib.pyplot as plt

from ase.io import read
from calorine.calculators import CPUNEP
from calorine.tools import relax_structure, get_force_constants

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# ================= 0. 输入文件 =================
structure_file = os.path.join(BASE_DIR, "POSCAR")
nep_file = os.path.join(BASE_DIR, "Si_2025_Xuke.txt")
kappa_file = os.path.join(BASE_DIR, "kappa-m21211.hdf5")

target_temp = 300.0
supercell_fc2 = [4, 4, 1]
points_per_segment = 50

# ================= 1. 设置期刊出版级风格 =================
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 20,
    'axes.linewidth': 1.5,
    'axes.labelsize': 28,
    'xtick.major.size': 8,
    'xtick.major.width': 1.5,
    'xtick.direction': 'in',
    'ytick.major.size': 8,
    'ytick.major.width': 1.5,
    'ytick.direction': 'in',
    'lines.linewidth': 2,
    'mathtext.fontset': 'stix',
    'legend.fontsize': 18,
    'legend.frameon': False,
})


# ================= 2. 工具函数：高对称路径 =================
def interpolate_k_path(kpoints, points_per_segment):
    """Generate interpolated k-points along a custom path."""
    k_path = []
    for i in range(len(kpoints) - 1):
        start = np.array(kpoints[i], dtype=float)
        end = np.array(kpoints[i + 1], dtype=float)
        segment = [
            (1 - t) * start + t * end
            for t in np.linspace(0, 1, points_per_segment, endpoint=False)
        ]
        k_path.extend(segment)
    k_path.append(np.array(kpoints[-1], dtype=float))
    return np.array(k_path)


def prepare_custom_path(kpoints, points_per_segment, labels=None):
    """Prepare interpolated k-path and cumulative linear distance."""
    k_path = interpolate_k_path(kpoints, points_per_segment)

    linearcoord = [0.0]
    for i in range(1, len(k_path)):
        delta_k = np.linalg.norm(k_path[i] - k_path[i - 1])
        linearcoord.append(linearcoord[-1] + delta_k)
    linearcoord = np.array(linearcoord)

    if labels is None:
        labels = [r'$\Gamma$', r'$\mathrm{X}$', r'$\mathrm{M}$', r'$\Gamma$']

    tick_indices = [i * points_per_segment for i in range(len(kpoints) - 1)] + [len(k_path) - 1]
    tick_positions = linearcoord[tick_indices]

    return {
        'explicit_kpoints_rel': k_path,
        'explicit_kpoints_linearcoord': linearcoord,
        'explicit_kpoints_labels': labels,
        'tick_indices': tick_indices,
        'tick_positions': tick_positions,
    }


def get_tau(g_data, tidx):
    """Convert gamma to relaxation time tau = 1 / (2 gamma)."""
    if g_data is None:
        return None
    g_vals = g_data[tidx, :, :].flatten()
    with np.errstate(divide='ignore', invalid='ignore'):
        t_vals = 1.0 / (2.0 * g_vals)
        t_vals[g_vals < 1e-10] = 0.0
    return t_vals


def filter_points(f, t):
    """Filter valid points for plotting relaxation time."""
    if t is None:
        return np.array([]), np.array([])
    mask = (t > 0) & (t < 1e6) & (f > 0.1)
    return f[mask], t[mask]


# ================= 3. 读取结构 =================
if not os.path.exists(structure_file):
    raise FileNotFoundError(f"找不到结构文件: {structure_file}")
if not os.path.exists(nep_file):
    raise FileNotFoundError(f"找不到 NEP 势文件: {nep_file}")
if not os.path.exists(kappa_file):
    raise FileNotFoundError(f"找不到 HDF5 文件: {kappa_file}")

structure = read(structure_file)


# ================= 4. 读取 HDF5：DOS / 群速度 / 弛豫时间 =================
with h5py.File(kappa_file, 'r') as f:
    frequencies = f['frequency'][:]          # (nq, nb)
    weights = f['weight'][:]                 # (nq,)
    temperatures = f['temperature'][:]       # (nt,)
    gamma_data = f['gamma'][:]               # (nt, nq, nb)
    gv_data = f['group_velocity'][:]         # (nq, nb, 3)

    gamma_N_data = f['gamma_N'][:] if 'gamma_N' in f else None
    gamma_U_data = f['gamma_U'][:] if 'gamma_U' in f else None

t_idx = (np.abs(temperatures - target_temp)).argmin()
print(f"Selected Temperature for spectra: {temperatures[t_idx]:.1f} K")


# ================= 5. 从 HDF5 构造 DOS / v_g / tau =================
freq_min, freq_max = np.min(frequencies) - 0.5, np.max(frequencies) + 0.5
n_qpoints, n_bands = frequencies.shape

flat_freqs = frequencies.flatten()
flat_weights = np.repeat(weights, n_bands)

# DOS
hist, bin_edges = np.histogram(
    flat_freqs,
    bins=200,
    range=(freq_min, freq_max),
    weights=flat_weights
)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
dos_smeared = hist.astype(float)
if np.max(dos_smeared) > 0:
    dos_smeared /= np.max(dos_smeared)

# Group velocity
# 默认 group_velocity 单位是 Å/ps，则乘 0.1 变成 km/s
gv_flat = 0.1 * np.linalg.norm(gv_data, axis=2).flatten()

# Relaxation time
tau_total = get_tau(gamma_data, t_idx)
tau_N = get_tau(gamma_N_data, t_idx)
tau_U = get_tau(gamma_U_data, t_idx)

f_tau_tot, t_tau_tot = filter_points(flat_freqs, tau_total)
f_tau_N, t_tau_N = filter_points(flat_freqs, tau_N)
f_tau_U, t_tau_U = filter_points(flat_freqs, tau_U)


# ================= 6. 重新计算高对称路径色散 =================
custom_kpoints = [
    [0.0, 0.0, 0.0],   # Γ
    [0.5, 0.0, 0.0],   # X
    [0.5, 0.5, 0.0],   # M
    [0.0, 0.0, 0.0],   # Γ
]

custom_path = prepare_custom_path(
    custom_kpoints,
    points_per_segment,
    labels=[r'$\Gamma$', r'$\mathrm{X}$', r'$\mathrm{M}$', r'$\Gamma$']
)

calculator = CPUNEP(nep_file)
structure.calc = calculator

# 结构弛豫
relax_structure(structure, fmax=0.0001)

# 二阶力常数
phonon = get_force_constants(structure, calculator, supercell_fc2)

# 高对称路径色散
phonon.run_band_structure([custom_path['explicit_kpoints_rel']])
band = phonon.get_band_structure_dict()

band_distances = np.array(custom_path['explicit_kpoints_linearcoord'])
band_frequencies = np.array(band['frequencies'][0])   # shape: (n_path, n_bands)

max_freq_hdf5 = np.max(frequencies)
max_freq_band = np.max(band_frequencies)
FREQ_MAX = int(np.ceil(1.1 * max(max_freq_hdf5, max_freq_band)))
print(f"Auto FREQ_MAX = {FREQ_MAX} THz")

# ================= 7. 绘图 (2x2) =================
fig, axes = plt.subplots(2, 2, figsize=(12, 10), constrained_layout=True)
ax_disp = axes[0, 0]
ax_dos  = axes[0, 1]
ax_gv   = axes[1, 0]
ax_tau  = axes[1, 1]

# --- (a) Dispersion ---
for ib in range(band_frequencies.shape[1]):
    ax_disp.plot(
        band_distances,
        band_frequencies[:, ib],
        color='#1f77b4',
        alpha=0.9,
        linewidth=1.5
    )

for x in custom_path['tick_positions']:
    ax_disp.axvline(x=x, color='k', linestyle='-', linewidth=1.0)

ax_disp.set_xticks(custom_path['tick_positions'])
ax_disp.set_xticklabels(custom_path['explicit_kpoints_labels'])
ax_disp.set_xlim(band_distances.min(), band_distances.max())
ax_disp.set_ylim(0, FREQ_MAX)
ax_disp.set_ylabel(r'$\omega~(\mathrm{THz})$')
ax_disp.set_xlabel('')
ax_disp.text(0.03, 0.93, '(a)', transform=ax_disp.transAxes,
             fontweight='bold', fontsize=20)

# --- (b) DOS ---
ax_dos.plot(dos_smeared, bin_centers, 'k-', linewidth=1.5)
ax_dos.fill_betweenx(bin_centers, 0, dos_smeared, color='gray', alpha=0.3)
ax_dos.set_ylim(0, FREQ_MAX)
ax_dos.set_xlim(0, 1.5)
ax_dos.set_xticks([])
ax_dos.set_xlabel('DOS (a.u.)')
ax_dos.set_ylabel(r'$\omega~(\mathrm{THz})$')
ax_dos.text(0.03, 0.93, '(b)', transform=ax_dos.transAxes,
            fontweight='bold', fontsize=20)

# --- (c) Group velocity ---
ax_gv.scatter(
    flat_freqs, gv_flat,
    s=20,
    facecolors='none',
    edgecolors='darkorange',
    linewidths=0.8,
    alpha=0.6
)
ax_gv.set_xlabel(r'$\omega~(\mathrm{THz})$')
ax_gv.set_ylabel(r'$v_g~(\mathrm{km/s})$')
ax_gv.set_xlim(0, FREQ_MAX)
ax_gv.set_ylim(0, 5)
ax_gv.grid(True, linestyle='--', alpha=0.3)
ax_gv.text(0.03, 0.93, '(c)', transform=ax_gv.transAxes,
           fontweight='bold', fontsize=20)

# --- (d) Relaxation time ---
has_nu = (len(f_tau_N) > 0 and len(f_tau_U) > 0)

if has_nu:
    ax_tau.scatter(
        f_tau_N, t_tau_N,
        s=20, facecolors='none', edgecolors='#1f77b4',
        linewidths=0.8, alpha=0.6,
        label=r'$\mathrm{Normal}$'
    )
    ax_tau.scatter(
        f_tau_U, t_tau_U,
        s=20, facecolors='none', edgecolors='crimson',
        linewidths=0.8, alpha=0.6,
        label=r'$\mathrm{Umklapp}$'
    )
    ax_tau.legend(loc='upper right')
else:
    if len(f_tau_tot) > 0:
        ax_tau.scatter(
            f_tau_tot, t_tau_tot,
            s=20, facecolors='none', edgecolors='k',
            linewidths=0.8, alpha=0.6,
            label=r'$\mathrm{Total}$'
        )
        ax_tau.legend(loc='upper right')

ax_tau.set_xlabel(r'$\omega~(\mathrm{THz})$')
ax_tau.set_ylabel(r'$\tau~(\mathrm{ps})$')
ax_tau.set_yscale('log')
ax_tau.set_xlim(0, FREQ_MAX)
ax_tau.set_ylim(1, 1e5)
ax_tau.grid(True, which='both', linestyle='--', alpha=0.3)
ax_tau.text(0.03, 0.93, '(d)', transform=ax_tau.transAxes,
            fontweight='bold', fontsize=20)

# ================= 8. 保存 =================
output_file = os.path.join(BASE_DIR, 'Film_1nm.pdf')
plt.savefig(output_file, dpi=300)
print(f"Done. Saved to {output_file}")