"""Plotting utilities for NEP-kappa results."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

cache_dir = Path(tempfile.gettempdir()) / "nepkappa-matplotlib"
cache_dir.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(cache_dir))
os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))

import h5py
import matplotlib

matplotlib.use("Agg")
matplotlib.set_loglevel("error")

import matplotlib.pyplot as plt
import numpy as np
import seekpath
from phono3py.interface.phono3py_yaml import Phono3pyYaml
from phonopy import Phonopy
from phonopy.file_IO import read_force_constants_hdf5


EV_TO_J = 1.602176634e-19
ANGSTROM3_TO_M3 = 1.0e-30
THZ_ANGSTROM_TO_KM_PER_S = 0.1
DEFAULT_FIGURES = [
    "dispersion",
    "dos",
    "heat_capacity",
    "group_velocity",
    "relaxation_time",
    "kappa",
]


def plot_results(config):
    """Create standard phonon and thermal-transport plots."""
    result_dir = Path(config.result_dir).resolve()
    plot_dir = result_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    disp_path = result_dir / "phono3py_disp.yaml"
    fc2_path = result_dir / "fc2.hdf5"
    kappa_path = find_kappa_file(result_dir, config.mesh)
    missing = [path for path in (disp_path, fc2_path, kappa_path) if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing file(s) required for plotting: "
            + ", ".join(str(path) for path in missing)
        )

    figures = DEFAULT_FIGURES
    layout = getattr(config, "plot_layout", "separate")
    dpi = int(getattr(config, "plot_dpi", 300))

    print("\n[Step 4] Plot Results")
    print(f"  - Reading {disp_path}")
    print(f"  - Reading {fc2_path}")
    print(f"  - Reading {kappa_path}")
    print(f"  - Plot layout: {layout}")
    print(f"  - Band path: {getattr(config, 'plot_path', 'seekpath')}")
    print(f"  - Plot figures: {', '.join(figures)}")

    phono3py_yaml = load_phono3py_yaml(disp_path)
    phonon = make_phonopy(phono3py_yaml, fc2_path)
    plot_data = build_plot_data(phonon, phono3py_yaml.unitcell, kappa_path, config)

    saved = []
    with journal_style():
        if layout in ("separate", "both"):
            saved.extend(write_separate_figures(figures, plot_data, plot_dir, dpi))
        if layout in ("combined", "both"):
            saved.append(write_combined_figure(figures, plot_data, plot_dir, dpi))

    print("  - Generated plots:")
    for path in saved:
        print(f"    {path}")


def journal_style():
    """Return matplotlib style settings suitable for publication figures."""
    return plt.rc_context(
        {
            "font.size": 14,
            "axes.labelsize": 16,
            "axes.titlesize": 16,
            "xtick.labelsize": 13,
            "ytick.labelsize": 13,
            "legend.fontsize": 12,
            "axes.linewidth": 1.2,
            "lines.linewidth": 1.8,
            "xtick.major.width": 1.1,
            "ytick.major.width": 1.1,
            "xtick.minor.width": 1.0,
            "ytick.minor.width": 1.0,
            "xtick.major.size": 5.0,
            "ytick.major.size": 5.0,
            "savefig.dpi": 300,
        }
    )


def find_kappa_file(result_dir, mesh):
    """Find the kappa HDF5 file matching the configured mesh when possible."""
    expected = result_dir / f"kappa-m{int(mesh[0])}{int(mesh[1])}{int(mesh[2])}.hdf5"
    if expected.exists():
        return expected
    candidates = sorted(result_dir.glob("kappa-m*.hdf5"))
    if candidates:
        return candidates[-1]
    return expected


def load_phono3py_yaml(path):
    """Load phono3py_disp.yaml."""
    ph3yml = Phono3pyYaml()
    ph3yml.read(str(path))
    return ph3yml


def make_phonopy(ph3yml, fc2_path):
    """Create a Phonopy object from phono3py metadata and fc2.hdf5."""
    phonon = Phonopy(
        ph3yml.unitcell,
        supercell_matrix=ph3yml.supercell_matrix,
        primitive_matrix=ph3yml.primitive_matrix,
    )
    phonon.force_constants = read_force_constants_hdf5(str(fc2_path))
    return phonon


def build_plot_data(phonon, unitcell, kappa_path, config):
    """Build all data needed by the plotting functions."""
    data = {
        "band": build_band_data(phonon, unitcell, config),
        "dos": build_dos_data(phonon, config.mesh),
        "transport": read_transport_data(kappa_path, unitcell, config),
    }
    return data


def seekpath_structure(unitcell):
    """Return the tuple format expected by seekpath."""
    return (unitcell.cell, unitcell.scaled_positions, unitcell.numbers)


def make_band_paths(unitcell, config, points_per_segment=51):
    """Build high-symmetry q-point paths with seekpath."""
    if getattr(config, "plot_path", "seekpath") == "custom":
        return make_custom_band_paths(config, points_per_segment)

    sp_path = seekpath.get_path(seekpath_structure(unitcell))
    point_coords = sp_path["point_coords"]
    paths = []
    labels = []
    for start_label, end_label in sp_path["path"]:
        start = np.array(point_coords[start_label], dtype=float)
        end = np.array(point_coords[end_label], dtype=float)
        segment = [
            start + (end - start) * i / (points_per_segment - 1)
            for i in range(points_per_segment)
        ]
        paths.append(segment)
        labels.append((format_label(start_label), format_label(end_label)))
    return paths, labels


def make_custom_band_paths(config, points_per_segment):
    """Build high-symmetry q-point paths from YAML."""
    point_coords = getattr(config, "plot_path_points", None)
    path_segments = getattr(config, "plot_path_segments", None)
    if not point_coords or not path_segments:
        raise ValueError(
            "plot.path: custom requires plot.path_points and plot.path_segments."
        )

    paths = []
    labels = []
    for segment in path_segments:
        if not isinstance(segment, (list, tuple)) or len(segment) != 2:
            raise ValueError(
                "Each custom plot.path_segments entry must contain two labels."
            )
        start_label, end_label = str(segment[0]), str(segment[1])
        if start_label not in point_coords or end_label not in point_coords:
            raise ValueError(
                f"Custom path segment {segment} uses undefined q-point label."
            )
        start = np.array(point_coords[start_label], dtype=float)
        end = np.array(point_coords[end_label], dtype=float)
        if start.shape != (3,) or end.shape != (3,):
            raise ValueError("Custom q-points must be three fractional coordinates.")
        qpoints = [
            start + (end - start) * i / (points_per_segment - 1)
            for i in range(points_per_segment)
        ]
        paths.append(qpoints)
        labels.append((format_label(start_label), format_label(end_label)))
    return paths, labels


def format_label(label):
    """Format seekpath labels for matplotlib."""
    if label.upper() in {"G", "GAMMA"}:
        return r"$\Gamma$"
    if "_" in label:
        head, tail = label.split("_", 1)
        return rf"{head}$_{{{tail}}}$"
    return label


def build_band_data(phonon, unitcell, config):
    """Compute band-structure data."""
    paths, labels = make_band_paths(unitcell, config)
    phonon.run_band_structure(paths, with_group_velocities=True)
    band = phonon.get_band_structure_dict()
    band["labels"] = labels
    return band


def build_dos_data(phonon, mesh):
    """Compute total DOS data."""
    phonon.run_mesh(mesh, is_gamma_center=True)
    phonon.run_total_dos()
    return phonon.get_total_dos_dict()


def read_transport_data(kappa_path, unitcell, config):
    """Read kappa HDF5 data and prepare derived transport quantities."""
    with h5py.File(kappa_path, "r") as handle:
        temperature = handle["temperature"][:]
        transport = {
            "temperature": temperature,
            "kappa": handle["kappa"][:],
            "heat_capacity": handle["heat_capacity"][:],
            "group_velocity": handle["group_velocity"][:],
            "frequency": handle["frequency"][:],
            "weight": handle["weight"][:],
            "mesh": handle["mesh"][:],
            "gamma": {
                "total": handle["gamma"][:],
            },
        }
        if "gamma_N" in handle:
            transport["gamma"]["normal"] = handle["gamma_N"][:]
        if "gamma_U" in handle:
            transport["gamma"]["umklapp"] = handle["gamma_U"][:]

    transport["volume_heat_capacity"] = volume_heat_capacity(
        transport["heat_capacity"],
        transport["weight"],
        transport["mesh"],
        unitcell,
    )
    transport["tau_temperature_index"] = int(
        np.argmin(np.abs(temperature - float(getattr(config, "plot_temperature", 300.0))))
    )
    transport["tau_mode"] = getattr(config, "plot_tau", "total")
    transport["kappa_mode"] = getattr(config, "plot_kappa", "all")
    return transport


def volume_heat_capacity(heat_capacity, weight, mesh, unitcell):
    """Return volumetric heat capacity in J m^-3 K^-1."""
    cell_volume_m3 = abs(np.linalg.det(unitcell.cell)) * ANGSTROM3_TO_M3
    mesh_size = int(np.prod(mesh))
    weighted_cv = np.sum(heat_capacity * weight[None, :, None], axis=(1, 2))
    return weighted_cv / mesh_size * EV_TO_J / cell_volume_m3


def write_separate_figures(figures, plot_data, plot_dir, dpi):
    """Write one PNG per requested figure."""
    saved = []
    sizes = {
        "dispersion": (7.4, 4.8),
        "dos": (5.6, 4.6),
        "heat_capacity": (5.8, 4.6),
        "group_velocity": (5.8, 4.6),
        "relaxation_time": (5.8, 4.6),
        "kappa": (5.8, 4.6),
    }
    for name in figures:
        fig, ax = plt.subplots(figsize=sizes[name], constrained_layout=True)
        draw_figure(name, ax, plot_data)
        path = plot_dir / f"{name}.png"
        fig.savefig(path, dpi=dpi)
        plt.close(fig)
        saved.append(path)
    return saved


def write_combined_figure(figures, plot_data, plot_dir, dpi):
    """Write all requested figures into one multi-panel PNG."""
    ncols = 3 if len(figures) == 6 else min(len(figures), 3)
    nrows = int(np.ceil(len(figures) / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(5.8 * ncols, 4.5 * nrows),
        constrained_layout=True,
    )
    axes = np.atleast_1d(axes).ravel()
    for index, name in enumerate(figures):
        draw_figure(name, axes[index], plot_data)
    for ax in axes[len(figures):]:
        ax.axis("off")

    path = plot_dir / "combined.png"
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    return path


def draw_figure(name, ax, plot_data):
    """Draw one named figure on an existing axis."""
    if name == "dispersion":
        draw_dispersion(ax, plot_data["band"])
    elif name == "dos":
        draw_dos(ax, plot_data["dos"])
    elif name == "heat_capacity":
        draw_heat_capacity(ax, plot_data["transport"])
    elif name == "group_velocity":
        draw_group_velocity(ax, plot_data["transport"])
    elif name == "relaxation_time":
        draw_relaxation_time(ax, plot_data["transport"])
    elif name == "kappa":
        draw_kappa(ax, plot_data["transport"])
    else:
        raise ValueError(f"Unknown plot figure: {name}")


def draw_dispersion(ax, band):
    """Draw phonon dispersion along seekpath high-symmetry lines."""
    tick_positions = []
    tick_labels = []
    previous_end_label = None
    for i, (distances, frequencies) in enumerate(
        zip(band["distances"], band["frequencies"])
    ):
        distances = np.asarray(distances)
        frequencies = np.asarray(frequencies)
        start_label, end_label = band["labels"][i]
        ax.plot(distances, frequencies, color="tab:blue", linewidth=1.5)
        if i == 0:
            tick_positions.append(distances[0])
            tick_labels.append(start_label)
        elif previous_end_label != start_label:
            if tick_positions and np.isclose(tick_positions[-1], distances[0]):
                tick_labels[-1] = f"{tick_labels[-1]}|{start_label}"
            else:
                tick_positions.append(distances[0])
                tick_labels.append(start_label)
        tick_positions.append(distances[-1])
        tick_labels.append(end_label)
        ax.axvline(distances[-1], color="0.82", linewidth=1.0)
        previous_end_label = end_label

    ax.set_xlim(tick_positions[0], tick_positions[-1])
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)
    ax.set_ylabel("Frequency (THz)")
    ax.grid(axis="y", color="0.9", linewidth=0.9)


def draw_dos(ax, dos):
    """Draw total phonon DOS."""
    ax.plot(dos["total_dos"], dos["frequency_points"], color="tab:green")
    ax.set_xlabel("DOS")
    ax.set_ylabel("Frequency (THz)")
    ax.grid(color="0.9", linewidth=0.9)


def draw_heat_capacity(ax, transport):
    """Draw volumetric heat capacity."""
    ax.plot(
        transport["temperature"],
        transport["volume_heat_capacity"],
        marker="o",
        color="tab:red",
    )
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel(r"Volumetric heat capacity (J m$^{-3}$ K$^{-1}$)")
    ax.grid(color="0.9", linewidth=0.9)


def draw_group_velocity(ax, transport):
    """Draw group-velocity magnitude versus frequency."""
    frequency = transport["frequency"]
    group_velocity = transport["group_velocity"] * THZ_ANGSTROM_TO_KM_PER_S
    gv_norm = np.linalg.norm(group_velocity, axis=2)
    valid = np.isfinite(frequency) & np.isfinite(gv_norm) & (frequency > 0)

    ax.scatter(frequency[valid], gv_norm[valid], s=12, alpha=0.35, color="tab:purple")
    ax.set_xlabel("Frequency (THz)")
    ax.set_ylabel("Group velocity (km/s)")
    ax.grid(color="0.9", linewidth=0.9)


def draw_relaxation_time(ax, transport):
    """Draw relaxation time from total, normal, and/or Umklapp scattering."""
    frequency = transport["frequency"]
    temp_index = transport["tau_temperature_index"]
    tau_mode = transport["tau_mode"]
    channels = tau_channels(tau_mode, transport["gamma"])
    colors = {
        "total": "tab:orange",
        "normal": "tab:blue",
        "umklapp": "tab:green",
    }
    labels = {
        "total": "total",
        "normal": "N",
        "umklapp": "U",
    }

    for name in channels:
        gamma = transport["gamma"][name][temp_index]
        tau = np.full_like(gamma, np.nan, dtype=float)
        positive = gamma > 0
        tau[positive] = 1.0 / (4.0 * np.pi * gamma[positive])
        valid = np.isfinite(frequency) & np.isfinite(tau) & (frequency > 0)
        ax.scatter(
            frequency[valid],
            tau[valid],
            s=12,
            alpha=0.35,
            color=colors[name],
            label=labels[name],
        )

    ax.set_yscale("log")
    ax.set_xlabel("Frequency (THz)")
    ax.set_ylabel("Relaxation time (ps)")
    if len(channels) > 1:
        ax.legend(frameon=False)
    ax.grid(color="0.9", linewidth=0.9)


def tau_channels(tau_mode, gamma_data):
    """Return relaxation-time channels requested by YAML."""
    if tau_mode == "all":
        return [name for name in ("total", "normal", "umklapp") if name in gamma_data]
    if tau_mode not in gamma_data:
        raise ValueError(
            f"Relaxation-time channel '{tau_mode}' is not available in kappa HDF5."
        )
    return [tau_mode]


def draw_kappa(ax, transport):
    """Draw selected thermal conductivity tensor diagonal components."""
    temperature = transport["temperature"]
    kappa = transport["kappa"]
    kxx, kyy, kzz = kappa[:, 0], kappa[:, 1], kappa[:, 2]
    kavg = (kxx + kyy + kzz) / 3.0
    mode = transport["kappa_mode"]
    components = {
        "x": (kxx, "xx", "o"),
        "y": (kyy, "yy", "s"),
        "z": (kzz, "zz", "^"),
    }

    if mode == "all":
        for _, (values, label, marker) in components.items():
            ax.plot(temperature, values, marker=marker, label=label)
        ax.plot(temperature, kavg, color="black", linewidth=2.4, label="average")
    else:
        values, label, marker = components[mode]
        ax.plot(temperature, values, marker=marker, label=label)

    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel(r"Thermal conductivity (W m$^{-1}$ K$^{-1}$)")
    ax.legend(frameon=False)
    ax.grid(color="0.9", linewidth=0.9)
