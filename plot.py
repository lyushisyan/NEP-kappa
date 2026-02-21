import argparse
import numpy as np
import matplotlib.pyplot as plt
import h5py
import warnings
import os

# --- ASE & Calorine (用于 NEP) ---
from ase.io import read
from calorine.calculators import CPUNEP
from calorine.tools import get_force_constants, relax_structure

# --- Seekpath (用于路径) ---
from seekpath import get_explicit_k_path

# 忽略不必要的警告
warnings.filterwarnings("ignore")

# ================= 1. 设置期刊出版级风格 =================
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 20,              # 全局字号
    'axes.linewidth': 1.5,
    'axes.labelsize': 30,         # 坐标轴标签字号
    'xtick.major.size': 8,
    'xtick.major.width': 1.5,
    'xtick.direction': 'in',
    'ytick.major.size': 8,
    'ytick.major.width': 1.5,
    'ytick.direction': 'in',
    'lines.linewidth': 2,
    'mathtext.fontset': 'stix',   # 让数学符号更美观
    'legend.fontsize': 20,
    'legend.frameon': False,      # 去掉图例边框
})

# ================= 2. 参数配置 =================

DEFAULT_DOS_MESH = [50, 50, 50]
DEFAULT_FREQ_MAX = 16.0
DEFAULT_TARGET_TEMP = 300.0
DEFAULT_TRUNCATE_LABEL = "L"
DEFAULT_OUTPUT = "NEP_results.pdf"

DEFAULT_NEP_CONFIG = {
    "POSCAR": "NEP/POSCAR",
    "MODEL": "NEP/Si_2025_Xuke.txt",
    "H5": "NEP/kappa-m212121.hdf5",
    "SUPERCELL": [4, 4, 4],
    "LABEL": r"$\mathrm{NEP}$",
    "COLOR": "#1f77b4",
    "MARKER": "o",
}

# ================= 3. 通用工具函数 =================

def get_any(h5, names):
    for n in names:
        if n in h5: return h5[n][...]
    return None

def calculate_tau(gamma):
    # Tau = 1 / (2 * Gamma)
    with np.errstate(divide='ignore', invalid='ignore'):
        t = 1.0 / (2.0 * gamma)
        t[gamma < 1e-10] = 0
    return t

def load_hdf5_data(filename, target_temp=300.0):
    """读取 HDF5 获取 Gv 和 Tau"""
    if not os.path.exists(filename):
        print(f"[Warning] File {filename} not found.")
        return None, None, None

    with h5py.File(filename, 'r') as f:
        freq = get_any(f, ["frequency", "frequencies"])
        if freq is None: return None, None, None
        
        gv = get_any(f, ["group_velocity", "gv"])
        if gv is not None:
            gv_norm = np.linalg.norm(gv, axis=2)
        else:
            gv_norm = None

        gamma = get_any(f, ["gamma", "scattering_rate"])
        temperatures = get_any(f, ["temperature"])
        
        tau = None
        if gamma is not None and temperatures is not None:
            idx = (np.abs(temperatures - target_temp)).argmin()
            gamma_t = gamma[idx]
            tau = calculate_tau(gamma_t)
            
        return freq.flatten(), gv_norm.flatten(), tau.flatten()

def save_data(filename, header, data_tuple):
    """保存数据"""
    try:
        data = np.column_stack(data_tuple)
        np.savetxt(filename, data, header=header, fmt='%.6f')
        print(f"[Export] Saved: {filename}")
    except Exception as e:
        print(f"[Error] Failed to save {filename}: {e}")

def get_kpath(structure, truncate_label=None):
    """统一获取 KPath"""
    st_tuple = (structure.cell, structure.get_scaled_positions(), structure.numbers)
    path_data = get_explicit_k_path(st_tuple)
    k_lin = path_data['explicit_kpoints_linearcoord']
    k_lab = path_data['explicit_kpoints_labels']
    k_rel = path_data['explicit_kpoints_rel']
    
    if truncate_label and truncate_label in k_lab:
        try:
            cut = k_lab.index(truncate_label)
            return k_lin[:cut+1], k_lab[:cut+1], k_rel[:cut+1]
        except ValueError:
            pass
    return k_lin, k_lab, k_rel

# ================= 4. 数据处理逻辑 =================

def process_nep(conf, dos_mesh, freq_max, target_temp, truncate_label):
    print(f"\n--- Processing NEP Data ({conf['POSCAR']}) ---")
    results = {}
    
    try:
        if not os.path.exists(conf['POSCAR']):
            print(f"Error: {conf['POSCAR']} not found.")
            return None

        # 1. 计算 Phonon
        structure = read(conf['POSCAR'])
        calc = CPUNEP(conf['MODEL'])
        structure.calc = calc
        relax_structure(structure, fmax=0.0001) 
        
        phonon = get_force_constants(structure, calc, conf['SUPERCELL'])
        
        # 2. 色散
        k_lin, k_lab, k_rel = get_kpath(structure, truncate_label)
        phonon.run_band_structure([k_rel], with_eigenvectors=False, labels=k_lab)
        bs = phonon.get_band_structure_dict()
        results['freq_disp'] = bs['frequencies'][0]
        results['k_lin'] = k_lin
        results['k_lab'] = k_lab
        
        # 3. DOS
        phonon.run_mesh(dos_mesh)
        phonon.run_total_dos(freq_min=0, freq_max=freq_max * 1.1, freq_pitch=0.05)
        dos = phonon.get_total_dos_dict()
        results['freq_dos'] = dos['frequency_points']
        results['val_dos'] = dos['total_dos']
        
        # 4. Particles (H5)
        f, g, t = load_hdf5_data(conf['H5'], target_temp)
        if f is not None:
            mask = (f > 0.1) & (f < freq_max * 1.5)
            if t is not None: mask &= (t > 0) & (t < 1e5)
            
            results['f_part'] = f[mask]
            results['gv'] = g[mask] * 0.1 if g is not None else None # A/ps -> km/s
            results['tau'] = t[mask] if t is not None else None
            
    except Exception as e:
        print(f"Error in NEP processing: {e}")
        return None
        
    return results

# ================= 5. 参数与主程序 =================

def build_parser():
    parser = argparse.ArgumentParser(description="Plot NEP phonon results.")
    parser.add_argument("--poscar", default=DEFAULT_NEP_CONFIG["POSCAR"], help="POSCAR path")
    parser.add_argument("--model", default=DEFAULT_NEP_CONFIG["MODEL"], help="NEP model file path")
    parser.add_argument("--h5", default=DEFAULT_NEP_CONFIG["H5"], help="kappa hdf5 file path")
    parser.add_argument(
        "--supercell",
        type=int,
        nargs=3,
        default=DEFAULT_NEP_CONFIG["SUPERCELL"],
        metavar=("NX", "NY", "NZ"),
        help="Supercell for force constants",
    )
    parser.add_argument(
        "--dos-mesh",
        type=int,
        nargs=3,
        default=DEFAULT_DOS_MESH,
        metavar=("MX", "MY", "MZ"),
        help="Mesh used for DOS",
    )
    parser.add_argument("--freq-max", type=float, default=DEFAULT_FREQ_MAX, help="Max frequency shown (THz)")
    parser.add_argument("--temp", type=float, default=DEFAULT_TARGET_TEMP, help="Target temperature for tau (K)")
    parser.add_argument(
        "--truncate-label",
        default=DEFAULT_TRUNCATE_LABEL,
        help="Stop band path at this label (set empty string to disable)",
    )
    parser.add_argument("--output", default=DEFAULT_OUTPUT, help="Output figure filename")
    return parser

def main():
    parser = build_parser()
    args = parser.parse_args()

    truncate_label = args.truncate_label if args.truncate_label else None
    nep_config = {
        "POSCAR": args.poscar,
        "MODEL": args.model,
        "H5": args.h5,
        "SUPERCELL": args.supercell,
        "LABEL": DEFAULT_NEP_CONFIG["LABEL"],
        "COLOR": DEFAULT_NEP_CONFIG["COLOR"],
        "MARKER": DEFAULT_NEP_CONFIG["MARKER"],
    }

    # --- 1. 获取数据 ---
    nep_data = process_nep(
        nep_config,
        dos_mesh=args.dos_mesh,
        freq_max=args.freq_max,
        target_temp=args.temp,
        truncate_label=truncate_label,
    )

    # --- 2. 绘图设置 ---
    print("\n--- Plotting NEP Results ---")
    # 使用 13x11 的尺寸，constrained_layout 自动调整布局
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), constrained_layout=True)
    ax_disp, ax_dos = axes[0, 0], axes[0, 1]
    ax_gv, ax_tau = axes[1, 0], axes[1, 1]

    # --- (a) Dispersion ---
    if nep_data:
        k_lin = nep_data['k_lin']
        freqs = nep_data['freq_disp']
        for i in range(freqs.shape[1]):
            label = nep_config['LABEL'] if i == 0 else None
            ax_disp.plot(k_lin, freqs[:, i], color=nep_config['COLOR'],
                         ls='-', alpha=0.9, lw=1.5, label=label)
            
        # 设置 ticks
        tick_pos, tick_labs = [], []
        for pos, lab in zip(nep_data['k_lin'], nep_data['k_lab']):
            if lab:
                ax_disp.axvline(x=pos, color='k', linestyle='-', linewidth=1.0)
                l = r'$\Gamma$' if lab == 'GAMMA' else lab
                tick_pos.append(pos)
                tick_labs.append(l)
        ax_disp.set_xticks(tick_pos)
        ax_disp.set_xticklabels(tick_labs)
        ax_disp.set_xlim(k_lin[0], k_lin[-1])

    # Y轴: Frequency -> \omega
    ax_disp.set_ylabel(r'$\omega,~\text{THz}$')
    ax_disp.set_xlabel(r'$~$')
    ax_disp.set_ylim(0, args.freq_max)
    ax_disp.text(0.03, 0.93, '(a)', transform=ax_disp.transAxes, fontweight='bold', fontsize=18)
    ax_disp.legend(loc='upper right', frameon=True)

    # --- (b) DOS ---
    if nep_data:
        ax_dos.plot(nep_data['val_dos'], nep_data['freq_dos'], 
                    color=nep_config['COLOR'], ls='-', label=nep_config['LABEL'])
        ax_dos.fill_betweenx(nep_data['freq_dos'], 0, nep_data['val_dos'], 
                             color=nep_config['COLOR'], alpha=0.15)

    ax_dos.set_ylim(0, args.freq_max)
    ax_dos.set_xlim(0, 4)
    ax_dos.set_xlabel(r'$\text{DOS, a.u.}$')
    ax_dos.set_xticks([])
    ax_dos.set_ylabel(r'$\omega,~\text{THz}$')
    ax_dos.text(0.03, 0.93, '(b)', transform=ax_dos.transAxes, fontweight='bold', fontsize=18)
    ax_dos.legend(loc='lower right', frameon=True)

    # --- (c) Group Velocity ---
    if nep_data and nep_data.get('gv') is not None:
        ax_gv.scatter(nep_data['f_part'], nep_data['gv'], s=20, 
                      facecolors='none', edgecolors=nep_config['COLOR'],
                      marker=nep_config['MARKER'], alpha=0.6, label=nep_config['LABEL'])
    
    # X轴: Frequency -> \omega, Y轴: Group Velocity -> v_g
    ax_gv.set_xlabel(r'$\omega,~\text{THz}$')
    ax_gv.set_ylabel(r'$v_g,~\text{km/s}$')
    ax_gv.set_xlim(0, args.freq_max)
    ax_gv.set_ylim(0, 10)
    ax_gv.grid(True, linestyle='--', alpha=0.3)
    ax_gv.text(0.03, 0.93, '(c)', transform=ax_gv.transAxes, fontweight='bold', fontsize=18)
    
    handles, labels = ax_gv.get_legend_handles_labels()
    if handles:
        ax_gv.legend(handles, labels, loc='upper right', frameon=True)

    # --- (d) Relaxation Time ---
    if nep_data and nep_data.get('tau') is not None:
        ax_tau.scatter(nep_data['f_part'], nep_data['tau'], s=20, 
                       facecolors='none', edgecolors=nep_config['COLOR'],
                       marker=nep_config['MARKER'], alpha=0.6, label=nep_config['LABEL'])

    # X轴: Frequency -> \omega, Y轴: Relaxation Time -> \tau
    ax_tau.set_xlabel(r'$\omega,~\text{THz}$')
    ax_tau.set_ylabel(r'$\tau,~\text{ps}$')
    ax_tau.set_yscale('log')
    ax_tau.set_xlim(0, args.freq_max)
    ax_tau.set_ylim(1, 1e5)
    ax_tau.grid(True, which='both', linestyle='--', alpha=0.3)
    ax_tau.text(0.03, 0.93, '(d)', transform=ax_tau.transAxes, fontweight='bold', fontsize=18)
    ax_tau.legend(loc='lower left', frameon=True)

    # --- 保存 ---
    outfile = args.output
    plt.savefig(outfile, dpi=300)
    print(f"\nDone. Plot saved to {outfile}")

if __name__ == "__main__":
    main()
