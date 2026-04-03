import h5py
import numpy as np
import matplotlib.pyplot as plt
import warnings
import os

# 忽略不必要的警告
warnings.filterwarnings("ignore")

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# ================= 1. 设置期刊出版级风格 (统一配置) =================
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

# ================= 2. 读取数据 =================
file_path = os.path.join(BASE_DIR, 'kappa-m21211.hdf5')
target_temp = 300.0

try:
    with h5py.File(file_path, 'r') as f:
        frequencies = f['frequency'][:]    
        weights = f['weight'][:]           
        temperatures = f['temperature'][:] 
        gamma_data = f['gamma'][:]         
        gv_data = f['group_velocity'][:]   
        
        gamma_N_data = f['gamma_N'][:] if 'gamma_N' in f else None
        gamma_U_data = f['gamma_U'][:] if 'gamma_U' in f else None

except FileNotFoundError:
    print(f"错误：找不到文件 {file_path}")
    exit()
except Exception as e:
    print(f"读取 HDF5 出错: {e}")
    exit()

t_idx = (np.abs(temperatures - target_temp)).argmin()
print(f"Selected Temperature for spectra: {temperatures[t_idx]:.1f} K")

# ================= 3. 数据处理与路径构建 =================
n_qpoints = frequencies.shape[0]
n_bands = frequencies.shape[1]

nx = int(np.sqrt(n_qpoints))
ny = nx
if nx * ny != n_qpoints:
    nx, ny = 11, 11 

try:
    freq_grid = frequencies.reshape(ny, nx, n_bands)

    path1 = freq_grid[0, :, :]             
    path2 = freq_grid[:, nx-1, :]          
    diag_idx = np.arange(nx)[::-1]
    path3 = freq_grid[diag_idx, diag_idx, :] 

    d1, d2 = 0.5, 0.5
    d3 = np.sqrt(0.5**2 + 0.5**2)
    
    n_points = nx
    x1 = np.linspace(0, d1, n_points)
    x2 = x1[-1] + np.linspace(0, d2, n_points)
    x3 = x2[-1] + np.linspace(0, d3, n_points)
    
    all_x = np.concatenate([x1, x2, x3])
    all_freqs = np.vstack([path1, path2, path3])
    
    plot_dispersion = True
except Exception as e:
    print(f"路径构建失败: {e}")
    plot_dispersion = False

# ================= 4. 物理量计算 =================
freq_min, freq_max = np.min(frequencies) - 0.5, np.max(frequencies) + 0.5
flat_weights = np.repeat(weights, n_bands)
flat_freqs = frequencies.flatten()

# (b) DOS
hist, bin_edges = np.histogram(flat_freqs, bins=200, range=(freq_min, freq_max), weights=flat_weights)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
dos_smeared = hist / np.max(hist)

# (c) 群速度
gv_flat = 0.1 * np.linalg.norm(gv_data, axis=2).flatten()

# 弛豫时间辅助函数
def get_tau(g_data, tidx):
    if g_data is None: return None
    g_vals = g_data[tidx, :, :].flatten()
    with np.errstate(divide='ignore', invalid='ignore'):
        t_vals = 1.0 / (2.0 * g_vals)
        t_vals[g_vals < 1e-10] = 0
    return t_vals

tau_total = get_tau(gamma_data, t_idx)
tau_N = get_tau(gamma_N_data, t_idx)
tau_U = get_tau(gamma_U_data, t_idx)

def filter_points(f, t):
    if t is None: return np.array([]), np.array([])
    mask = (t > 0) & (t < 1e6) & (f > 0.1)
    return f[mask], t[mask]

f_tau_tot, t_tau_tot = filter_points(flat_freqs, tau_total)
f_tau_N, t_tau_N = filter_points(flat_freqs, tau_N)
f_tau_U, t_tau_U = filter_points(flat_freqs, tau_U)

# ================= 5. 绘图 (2x2 布局) =================
fig, axes = plt.subplots(2, 2, figsize=(12, 10), constrained_layout=True)
ax_disp = axes[0, 0]
ax_dos  = axes[0, 1]
ax_gv   = axes[1, 0]
ax_tau  = axes[1, 1]

FREQ_MAX = 16.0 

# --- (a) Dispersion ---
if plot_dispersion:
    for i in range(n_bands):
        ax_disp.plot(all_x, all_freqs[:, i], color='#1f77b4', alpha=0.9, linewidth=1.5)
    
    ax_disp.axvline(x=x1[-1], color='k', linestyle='-', linewidth=1.0)
    ax_disp.axvline(x=x2[-1], color='k', linestyle='-', linewidth=1.0)
    ax_disp.set_xticks([0, x1[-1], x2[-1], x3[-1]])
    ax_disp.set_xticklabels([r'$\Gamma$', r'$\mathrm{X}$', r'$\mathrm{M}$', r'$\Gamma$'])
    ax_disp.set_xlim(0, x3[-1])
else:
    ax_disp.text(0.5, 0.5, "Path Error", ha='center')

ax_disp.set_ylim(0, 16)
ax_disp.set_ylabel(r'$\omega,~\text{THz}$') 
ax_disp.set_xlabel(r'$~$')
ax_disp.text(0.03, 0.93, '(a)', transform=ax_disp.transAxes, fontweight='bold', fontsize=20)

# --- (b) DOS ---
ax_dos.plot(dos_smeared, bin_centers, 'k-', linewidth=1.5)
ax_dos.fill_betweenx(bin_centers, 0, dos_smeared, color='gray', alpha=0.3)
ax_dos.set_ylim(0, 16)
ax_dos.set_xlim(0, 1.5)
ax_dos.set_xticks([])
ax_dos.set_xlabel(r'$\text{DOS, a.u.}$')
ax_dos.set_ylabel(r'$\omega,~\text{THz}$')
ax_dos.text(0.03, 0.93, '(b)', transform=ax_dos.transAxes, fontweight='bold', fontsize=20)

# --- (c) Group Velocity ---
ax_gv.scatter(flat_freqs, gv_flat, s=20, 
              facecolors='none', edgecolors='darkorange', linewidths=0.8, alpha=0.6)
ax_gv.set_xlabel(r'$\omega,~\text{THz}$')
ax_gv.set_ylabel(r'$v_g,~\text{km/s}$')
ax_gv.set_xlim(0, FREQ_MAX)
ax_gv.set_ylim(0, 5)
ax_gv.grid(True, linestyle='--', alpha=0.3)
ax_gv.text(0.03, 0.93, '(c)', transform=ax_gv.transAxes, fontweight='bold', fontsize=20)

# --- (d) Relaxation Time ---
has_nu = (len(f_tau_N) > 0 and len(f_tau_U) > 0)
if has_nu:
    ax_tau.scatter(f_tau_N, t_tau_N, s=20, facecolors='none', edgecolors='#1f77b4', linewidths=0.8, alpha=0.6, label=r'$\mathrm{Normal}$')
    ax_tau.scatter(f_tau_U, t_tau_U, s=20, facecolors='none', edgecolors='crimson', linewidths=0.8, alpha=0.6, label=r'$\mathrm{Umklapp}$')
    ax_tau.legend(loc='upper right')
else:
    if len(f_tau_tot) > 0:
        ax_tau.scatter(f_tau_tot, t_tau_tot, s=20, facecolors='none', edgecolors='k', linewidths=0.8, alpha=0.6, label=r'$\mathrm{Total}$')
        ax_tau.legend(loc='upper right')

ax_tau.set_xlabel(r'$\omega,~\text{THz}$')
ax_tau.set_ylabel(r'$\tau,~\text{ps}$')
ax_tau.set_yscale('log')
ax_tau.set_xlim(0, FREQ_MAX)
ax_tau.set_ylim(1, 1e5) 
ax_tau.grid(True, which='both', linestyle='--', alpha=0.3)
ax_tau.text(0.03, 0.93, '(d)', transform=ax_tau.transAxes, fontweight='bold', fontsize=20)

# 保存
output_file = os.path.join(BASE_DIR, 'Si_film_1nm.pdf')
plt.savefig(output_file, dpi=300)
print(f"Done. Saved to {output_file}")
