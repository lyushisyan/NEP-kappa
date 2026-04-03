import numpy as np
import matplotlib.pyplot as plt
import h5py
import warnings
import os

# --- ASE & Calorine (用于 NEP) ---
from ase.io import read
from calorine.calculators import CPUNEP
from calorine.tools import get_force_constants, relax_structure

# --- Phonopy & Seekpath (用于 DFT 和 路径) ---
from seekpath import get_explicit_k_path
from phonopy import load

# 忽略不必要的警告
warnings.filterwarnings("ignore")

BASE_DIR = os.path.dirname(os.path.abspath(__file__))


def in_example_dir(*parts):
    return os.path.join(BASE_DIR, *parts)

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
    'legend.fontsize': 20,
    'legend.frameon': True,       
})

# ================= 2. 参数配置 =================
DOS_MESH = [50, 50, 50]
FREQ_MAX = 16.0
TARGET_TEMP = 300.0
TRUNCATE_LABEL = 'L'  

# --- NEP/DFT 配置 ---
NEP_CONFIG = {
    'DIR': in_example_dir('NEP'),
    'POSCAR': in_example_dir('NEP', 'POSCAR'),
    'MODEL': in_example_dir('NEP', 'Si.txt'),
    'H5': in_example_dir('NEP', 'kappa-m212121.hdf5'),
    'SUPERCELL': [4, 4, 4],
    'LABEL': r'$\mathrm{NEP}$',
    'COLOR': '#1f77b4',
    'MARKER': 'o',
}
DFT_CONFIG = {
    'DIR': in_example_dir('DFT'),
    'POSCAR': in_example_dir('DFT', 'POSCAR'),
    'FC2': in_example_dir('DFT', 'fc2.hdf5'),
    'H5': in_example_dir('DFT', 'kappa-m212121.hdf5'),
    'SUPERCELL': [3, 3, 3],
    'LABEL': r'$\mathrm{DFT}$',
    'COLOR': '#d62728',
    'MARKER': 'x',
}

# ================= 3. 通用工具函数 =================
def get_any(h5, names):
    for n in names:
        if n in h5: return h5[n][...]
    return None

def calculate_tau(gamma):
    with np.errstate(divide='ignore', invalid='ignore'):
        t = 1.0 / (2.0 * gamma)
        t[gamma < 1e-10] = 0
    return t

def load_hdf5_data(filename, target_temp=300.0):
    if not os.path.exists(filename): return None, None, None
    with h5py.File(filename, 'r') as f:
        freq = get_any(f, ["frequency", "frequencies"])
        if freq is None: return None, None, None
        gv = get_any(f, ["group_velocity", "gv"])
        gv_norm = np.linalg.norm(gv, axis=2) if gv is not None else None
        gamma = get_any(f, ["gamma", "scattering_rate"])
        temperatures = get_any(f, ["temperature"])
        tau = None
        if gamma is not None and temperatures is not None:
            idx = (np.abs(temperatures - target_temp)).argmin()
            tau = calculate_tau(gamma[idx])
        return freq.flatten(), gv_norm.flatten() if gv_norm is not None else None, tau.flatten() if tau is not None else None

def get_kpath(structure, truncate_label=None):
    st_tuple = (structure.cell, structure.get_scaled_positions(), structure.numbers)
    path_data = get_explicit_k_path(st_tuple)
    k_lin = path_data['explicit_kpoints_linearcoord']
    k_lab = path_data['explicit_kpoints_labels']
    k_rel = path_data['explicit_kpoints_rel']
    if truncate_label and truncate_label in k_lab:
        try:
            cut = k_lab.index(truncate_label)
            return k_lin[:cut+1], k_lab[:cut+1], k_rel[:cut+1]
        except ValueError: pass
    return k_lin, k_lab, k_rel

# ================= 4. 数据处理逻辑 =================
def process_nep(conf):
    print(f"\n--- Processing NEP Data ({conf['POSCAR']}) ---")
    try:
        if not os.path.exists(conf['POSCAR']): return None
        structure = read(conf['POSCAR'])
        calc = CPUNEP(conf['MODEL'])
        structure.calc = calc
        relax_structure(structure, fmax=0.0001) 
        phonon = get_force_constants(structure, calc, conf['SUPERCELL'])
        
        k_lin, k_lab, k_rel = get_kpath(structure, TRUNCATE_LABEL)
        phonon.run_band_structure([k_rel], with_eigenvectors=False, labels=k_lab)
        bs = phonon.get_band_structure_dict()
        
        phonon.run_mesh(DOS_MESH)
        phonon.run_total_dos(freq_min=0, freq_max=FREQ_MAX*1.1, freq_pitch=0.05)
        dos = phonon.get_total_dos_dict()
        
        f, g, t = load_hdf5_data(conf['H5'], TARGET_TEMP)
        mask = (f > 0.1) & (f < FREQ_MAX * 1.5) if f is not None else None
        if mask is not None and t is not None: mask &= (t > 0) & (t < 1e5)
        
        return {
            'k_lin': k_lin, 'k_lab': k_lab, 'freq_disp': bs['frequencies'][0],
            'freq_dos': dos['frequency_points'], 'val_dos': dos['total_dos'],
            'f_part': f[mask] if mask is not None else None, 
            'gv': g[mask]*0.1 if g is not None and mask is not None else None, 
            'tau': t[mask] if t is not None and mask is not None else None
        }
    except Exception as e:
        print(f"Error in NEP processing: {e}")
        return None

def process_dft(conf):
    print(f"\n--- Processing DFT Data ({conf['FC2']}) ---")
    try:
        if not os.path.exists(conf['FC2']): return None
        structure = read(conf['POSCAR']) 
        phonon = load(unitcell_filename=conf['POSCAR'], supercell_matrix=conf['SUPERCELL'], force_constants_filename=conf['FC2'])
        
        k_lin, k_lab, k_rel = get_kpath(structure, TRUNCATE_LABEL)
        phonon.run_band_structure([k_rel], with_eigenvectors=False, labels=k_lab)
        bs = phonon.get_band_structure_dict()
        
        phonon.run_mesh(DOS_MESH)
        phonon.run_total_dos(freq_min=0, freq_max=FREQ_MAX*1.1, freq_pitch=0.05)
        dos = phonon.get_total_dos_dict()
        
        f, g, t = load_hdf5_data(conf['H5'], TARGET_TEMP)
        mask = (f > 0.1) & (f < FREQ_MAX * 1.5) if f is not None else None
        if mask is not None and t is not None: mask &= (t > 0) & (t < 1e5)
        
        return {
            'k_lin': k_lin[:bs['frequencies'][0].shape[0]], 'k_lab': k_lab, 'freq_disp': bs['frequencies'][0],
            'freq_dos': dos['frequency_points'], 'val_dos': dos['total_dos'],
            'f_part': f[mask] if mask is not None else None, 
            'gv': g[mask]*0.1 if g is not None and mask is not None else None, 
            'tau': t[mask] if t is not None and mask is not None else None
        }
    except Exception as e:
        print(f"Error in DFT processing: {e}")
        return None

# ================= 5. 主程序与绘图 =================
def main():
    nep_data = process_nep(NEP_CONFIG)
    dft_data = process_dft(DFT_CONFIG)

    print("\n--- Plotting Comparison ---")
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), constrained_layout=True)
    ax_disp, ax_dos = axes[0, 0], axes[0, 1]
    ax_gv, ax_tau = axes[1, 0], axes[1, 1]

    # --- (a) Dispersion ---
    if nep_data:
        for i in range(nep_data['freq_disp'].shape[1]):
            ax_disp.plot(nep_data['k_lin'], nep_data['freq_disp'][:, i], color=NEP_CONFIG['COLOR'], ls='-', lw=1.5, label=NEP_CONFIG['LABEL'] if i==0 else None)
        tick_pos, tick_labs = [], []
        for pos, lab in zip(nep_data['k_lin'], nep_data['k_lab']):
            if lab:
                ax_disp.axvline(x=pos, color='k', linestyle='-', linewidth=1.0)
                tick_pos.append(pos)
                tick_labs.append(r'$\Gamma$' if lab == 'GAMMA' else lab)
        ax_disp.set_xticks(tick_pos)
        ax_disp.set_xticklabels(tick_labs)
        ax_disp.set_xlim(nep_data['k_lin'][0], nep_data['k_lin'][-1])

    if dft_data:
        for i in range(dft_data['freq_disp'].shape[1]):
            ax_disp.plot(dft_data['k_lin'], dft_data['freq_disp'][:, i], color=DFT_CONFIG['COLOR'], ls='--', lw=1.5, label=DFT_CONFIG['LABEL'] if i==0 else None)

    ax_disp.set_ylabel(r'$\omega,~\text{THz}$'); ax_disp.set_ylim(0, 16.5)
    ax_disp.text(0.03, 0.92, '(a)', transform=ax_disp.transAxes, fontweight='bold', fontsize=20)
    ax_disp.legend(loc='upper right', frameon=True)

    # --- (b) DOS ---
    if nep_data:
        ax_dos.plot(nep_data['val_dos'], nep_data['freq_dos'], color=NEP_CONFIG['COLOR'], ls='-', label=NEP_CONFIG['LABEL'])
        ax_dos.fill_betweenx(nep_data['freq_dos'], 0, nep_data['val_dos'], color=NEP_CONFIG['COLOR'], alpha=0.15)
    if dft_data:
        ax_dos.plot(dft_data['val_dos'], dft_data['freq_dos'], color=DFT_CONFIG['COLOR'], ls='--', lw=1.5, label=DFT_CONFIG['LABEL'])

    ax_dos.set_xlim(0, 4); ax_dos.set_ylim(0, 16.5)
    ax_dos.set_xlabel(r'$\text{DOS, a.u.}$'); ax_dos.set_xticks([])
    ax_dos.text(0.03, 0.92, '(b)', transform=ax_dos.transAxes, fontweight='bold', fontsize=20)
    ax_dos.legend(loc='lower right', frameon=True)

    # --- (c) Group Velocity ---
    if nep_data and nep_data.get('gv') is not None:
        ax_gv.scatter(nep_data['f_part'], nep_data['gv'], s=30, facecolors='none', edgecolors=NEP_CONFIG['COLOR'], marker=NEP_CONFIG['MARKER'], alpha=0.6, label=NEP_CONFIG['LABEL'])
    if dft_data and dft_data.get('gv') is not None:
        ax_gv.scatter(dft_data['f_part'], dft_data['gv'], s=30, color=DFT_CONFIG['COLOR'], marker=DFT_CONFIG['MARKER'], alpha=0.6, label=DFT_CONFIG['LABEL'])

    ax_gv.set_xlabel(r'$\omega,~\text{THz}$'); ax_gv.set_ylabel(r'$v_g,~\text{km/s}$')
    ax_gv.set_xlim(0, FREQ_MAX); ax_gv.set_ylim(0, 10)
    ax_gv.grid(True, linestyle='--', alpha=0.3)
    ax_gv.text(0.03, 0.92, '(c)', transform=ax_gv.transAxes, fontweight='bold', fontsize=20)
    handles, labels = ax_gv.get_legend_handles_labels()
    if handles: ax_gv.legend(handles, labels, loc='upper right', frameon=True)

    # --- (d) Relaxation Time ---
    if nep_data and nep_data.get('tau') is not None:
        ax_tau.scatter(nep_data['f_part'], nep_data['tau'], s=30, facecolors='none', edgecolors=NEP_CONFIG['COLOR'], marker=NEP_CONFIG['MARKER'], alpha=0.6, label=NEP_CONFIG['LABEL'])
    if dft_data and dft_data.get('tau') is not None:
        ax_tau.scatter(dft_data['f_part'], dft_data['tau'], s=30, color=DFT_CONFIG['COLOR'], marker=DFT_CONFIG['MARKER'], alpha=0.6, label=DFT_CONFIG['LABEL'])

    ax_tau.set_xlabel(r'$\omega,~\text{THz}$'); ax_tau.set_ylabel(r'$\tau,~\text{ps}$')
    ax_tau.set_yscale('log'); ax_tau.set_xlim(0, FREQ_MAX); ax_tau.set_ylim(1, 1e5)
    ax_tau.grid(True, which='both', linestyle='--', alpha=0.3)
    ax_tau.text(0.03, 0.92, '(d)', transform=ax_tau.transAxes, fontweight='bold', fontsize=20)
    ax_tau.legend(loc='lower left', frameon=True)

    # --- 保存 ---
    outfile = in_example_dir('Comparison_NEP_vs_DFT.pdf')
    plt.savefig(outfile, dpi=300)
    print(f"\nDone. Plot saved to {outfile}")

if __name__ == "__main__":
    main()
