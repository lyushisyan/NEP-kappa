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
AXIS_INDEX = {"x": 0, "y": 1, "z": 2}


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
    correction = plot_data["transport"]["geometry_correction"]
    if correction["factor"] != 1.0:
        print(f"  - Effective geometry: {correction['description']}")
        print(f"  - Applying kappa/Cv correction factor: {correction['factor']:.8g}")

    saved = []
    with journal_style():
        if layout in ("separate", "both"):
            saved.extend(write_separate_figures(figures, plot_data, plot_dir, dpi))
        if layout in ("combined", "both"):
            saved.append(write_combined_figure(figures, plot_data, plot_dir, dpi))

    print("  - Generated plots:")
    for path in saved:
        print(f"    {path}")


def compare_results(config):
    """Create DFT-vs-NEP comparison plots."""
    dft_dir = Path(config.dft_dir).resolve()
    nep_dir = Path(config.nep_dir).resolve()
    compare_dir = Path(config.compare_dir).resolve()
    plot_dir = compare_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    figures = DEFAULT_FIGURES
    layout = getattr(config, "plot_layout", "separate")
    dpi = int(getattr(config, "plot_dpi", 300))

    print("\n[Compare] Plot DFT vs NEP Results")
    print(f"  - Reference ({config.reference_label}): {dft_dir}")
    print(f"  - Candidate ({config.candidate_label}): {nep_dir}")
    print(f"  - Compare directory: {compare_dir}")
    print(f"  - Plot layout: {layout}")
    print(f"  - Band path: {getattr(config, 'plot_path', 'seekpath')}")
    print(f"  - Plot figures: {', '.join(figures)}")

    datasets = [
        load_result_plot_data(dft_dir, config, config.reference_label),
        load_result_plot_data(nep_dir, config, config.candidate_label),
    ]
    for dataset in datasets:
        correction = dataset["transport"]["geometry_correction"]
        if correction["factor"] != 1.0:
            print(f"  - {dataset['label']} effective geometry: {correction['description']}")
            print(
                f"  - {dataset['label']} kappa/Cv correction factor: "
                f"{correction['factor']:.8g}"
            )

    saved = []
    with journal_style():
        if layout in ("separate", "both"):
            saved.extend(write_compare_separate_figures(figures, datasets, plot_dir, dpi))
        if layout in ("combined", "both"):
            saved.append(write_compare_combined_figure(figures, datasets, plot_dir, dpi))

    print("  - Generated comparison plots:")
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


def find_kappa_file(result_dir, mesh=None):
    """Find the kappa HDF5 file matching the configured mesh when possible."""
    if mesh is not None:
        expected = result_dir / f"kappa-m{int(mesh[0])}{int(mesh[1])}{int(mesh[2])}.hdf5"
        if expected.exists():
            return expected
    candidates = sorted(result_dir.glob("kappa-m*.hdf5"))
    if candidates:
        return candidates[-1]
    if mesh is None:
        return result_dir / "kappa-m*.hdf5"
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
        "transport": read_transport_data(
            kappa_path,
            phonon.primitive.volume,
            phonon.primitive.cell,
            config,
        ),
    }
    return data


def load_result_plot_data(result_dir, config, label):
    """Load one completed result directory for plotting or comparison."""
    disp_path = result_dir / "phono3py_disp.yaml"
    fc2_path = result_dir / "fc2.hdf5"
    kappa_path = find_kappa_file(result_dir, getattr(config, "mesh", None))
    missing = [path for path in (disp_path, fc2_path, kappa_path) if not path.exists()]
    if missing:
        raise FileNotFoundError(
            f"Missing file(s) required for {label}: "
            + ", ".join(str(path) for path in missing)
        )

    print(f"  - Reading {label}: {disp_path}")
    print(f"  - Reading {label}: {fc2_path}")
    print(f"  - Reading {label}: {kappa_path}")

    phono3py_yaml = load_phono3py_yaml(disp_path)
    phonon = make_phonopy(phono3py_yaml, fc2_path)
    transport = read_transport_data(
        kappa_path,
        phonon.primitive.volume,
        phonon.primitive.cell,
        config,
    )
    data = {
        "label": label,
        "band": build_band_data(phonon, phono3py_yaml.unitcell, config),
        "dos": build_dos_data(phonon, transport["mesh"]),
        "transport": transport,
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


def read_transport_data(kappa_path, primitive_volume, primitive_cell, config):
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

    correction = effective_geometry_correction(config, primitive_cell, primitive_volume)
    effective_volume = primitive_volume / correction["factor"]
    transport["kappa"] = transport["kappa"] * correction["factor"]
    transport["volume_heat_capacity"] = volume_heat_capacity(
        transport["heat_capacity"],
        transport["weight"],
        transport["mesh"],
        effective_volume,
    )
    transport["geometry_correction"] = correction
    transport["tau_temperature_index"] = int(
        np.argmin(np.abs(temperature - float(getattr(config, "plot_temperature", 300.0))))
    )
    transport["tau_mode"] = getattr(config, "plot_tau", "total")
    transport["kappa_mode"] = getattr(config, "plot_kappa", "all")
    return transport


def volume_heat_capacity(heat_capacity, weight, mesh, volume_angstrom3):
    """Return volumetric heat capacity in J m^-3 K^-1."""
    cell_volume_m3 = float(volume_angstrom3) * ANGSTROM3_TO_M3
    mesh_size = int(np.prod(mesh))
    weighted_cv = np.sum(heat_capacity * weight[None, :, None], axis=(1, 2))
    return weighted_cv / mesh_size * EV_TO_J / cell_volume_m3


def effective_geometry_correction(config, primitive_cell, primitive_volume):
    """Return the effective-volume correction for 2D films and 1D wires."""
    dimensionality = int(getattr(config, "dimensionality", 3))
    if dimensionality == 3:
        return {
            "dimensionality": 3,
            "factor": 1.0,
            "description": "3D bulk normalization",
        }

    cell = np.asarray(primitive_cell, dtype=float)
    volume = float(abs(primitive_volume))
    if dimensionality == 2:
        axis_name = getattr(config, "vacuum_axis", "z")
        axis = axis_to_index(axis_name, "vacuum_axis")
        thickness = require_positive(
            getattr(config, "effective_thickness", None),
            "effective_thickness",
        )
        in_plane = [index for index in range(3) if index != axis]
        in_plane_area = np.linalg.norm(np.cross(cell[in_plane[0]], cell[in_plane[1]]))
        cell_thickness = volume / in_plane_area
        factor = float(cell_thickness / thickness)
        return {
            "dimensionality": 2,
            "factor": factor,
            "description": (
                f"2D film, vacuum_axis={axis_name}, "
                f"cell_thickness={cell_thickness:.6g} A, "
                f"effective_thickness={thickness:.6g} A"
            ),
        }

    if dimensionality == 1:
        axis_name = getattr(config, "periodic_axis", "z")
        axis = axis_to_index(axis_name, "periodic_axis")
        effective_area = require_positive(
            getattr(config, "effective_area", None),
            "effective_area",
        )
        periodic_length = np.linalg.norm(cell[axis])
        cell_area = volume / periodic_length
        factor = float(cell_area / effective_area)
        return {
            "dimensionality": 1,
            "factor": factor,
            "description": (
                f"1D nanowire, periodic_axis={axis_name}, "
                f"cell_area={cell_area:.6g} A^2, "
                f"effective_area={effective_area:.6g} A^2"
            ),
        }

    raise ValueError("dimensionality must be 1, 2, or 3.")


def axis_to_index(axis, option_name):
    """Convert an axis label to a lattice-vector index."""
    try:
        return AXIS_INDEX[str(axis).lower()]
    except KeyError as exc:
        raise ValueError(f"{option_name} must be one of x, y, or z.") from exc


def require_positive(value, option_name):
    """Return a positive float or raise a clear plotting error."""
    if value is None:
        raise ValueError(f"{option_name} is required for this dimensionality.")
    value = float(value)
    if value <= 0:
        raise ValueError(f"{option_name} must be positive.")
    return value


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


def write_compare_separate_figures(figures, datasets, plot_dir, dpi):
    """Write one comparison PNG per figure."""
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
        draw_compare_figure(name, ax, datasets)
        path = plot_dir / f"{name}.png"
        fig.savefig(path, dpi=dpi)
        plt.close(fig)
        saved.append(path)
    return saved


def write_compare_combined_figure(figures, datasets, plot_dir, dpi):
    """Write all comparison figures into one multi-panel PNG."""
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
        draw_compare_figure(name, axes[index], datasets)
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


def draw_compare_figure(name, ax, datasets):
    """Draw one named DFT-vs-NEP comparison figure."""
    if name == "dispersion":
        draw_compare_dispersion(ax, datasets)
    elif name == "dos":
        draw_compare_dos(ax, datasets)
    elif name == "heat_capacity":
        draw_compare_heat_capacity(ax, datasets)
    elif name == "group_velocity":
        draw_compare_group_velocity(ax, datasets)
    elif name == "relaxation_time":
        draw_compare_relaxation_time(ax, datasets)
    elif name == "kappa":
        draw_compare_kappa(ax, datasets)
    else:
        raise ValueError(f"Unknown compare figure: {name}")


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


def draw_compare_dispersion(ax, datasets):
    """Draw overlaid phonon dispersions."""
    colors = ["tab:blue", "tab:orange"]
    linestyles = ["-", "--"]
    reference_band = datasets[0]["band"]
    set_band_ticks(ax, reference_band)

    for dataset, color, linestyle in zip(datasets, colors, linestyles):
        band = dataset["band"]
        first_line = True
        for distances, frequencies in zip(band["distances"], band["frequencies"]):
            ax.plot(
                np.asarray(distances),
                np.asarray(frequencies),
                color=color,
                linestyle=linestyle,
                linewidth=1.5,
                label=dataset["label"] if first_line else None,
            )
            first_line = False

    ax.set_ylabel("Frequency (THz)")
    ax.legend(frameon=False)
    ax.grid(axis="y", color="0.9", linewidth=0.9)


def set_band_ticks(ax, band):
    """Set high-symmetry ticks from a band-structure dictionary."""
    tick_positions = []
    tick_labels = []
    previous_end_label = None
    for i, distances in enumerate(band["distances"]):
        distances = np.asarray(distances)
        start_label, end_label = band["labels"][i]
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


def draw_compare_dos(ax, datasets):
    """Draw overlaid phonon DOS curves."""
    colors = ["tab:blue", "tab:orange"]
    linestyles = ["-", "--"]
    for dataset, color, linestyle in zip(datasets, colors, linestyles):
        dos = dataset["dos"]
        ax.plot(
            dos["total_dos"],
            dos["frequency_points"],
            color=color,
            linestyle=linestyle,
            label=dataset["label"],
        )
    ax.set_xlabel("DOS")
    ax.set_ylabel("Frequency (THz)")
    ax.legend(frameon=False)
    ax.grid(color="0.9", linewidth=0.9)


def draw_compare_heat_capacity(ax, datasets):
    """Draw overlaid volumetric heat capacities."""
    colors = ["tab:blue", "tab:orange"]
    markers = ["o", "s"]
    for dataset, color, marker in zip(datasets, colors, markers):
        transport = dataset["transport"]
        ax.plot(
            transport["temperature"],
            transport["volume_heat_capacity"],
            marker=marker,
            color=color,
            label=dataset["label"],
        )
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel(r"Volumetric heat capacity (J m$^{-3}$ K$^{-1}$)")
    ax.legend(frameon=False)
    ax.grid(color="0.9", linewidth=0.9)


def draw_compare_group_velocity(ax, datasets):
    """Draw overlaid group-velocity scatter plots."""
    colors = ["tab:blue", "tab:orange"]
    markers = ["o", "x"]
    for dataset, color, marker in zip(datasets, colors, markers):
        transport = dataset["transport"]
        frequency = transport["frequency"]
        group_velocity = transport["group_velocity"] * THZ_ANGSTROM_TO_KM_PER_S
        gv_norm = np.linalg.norm(group_velocity, axis=2)
        valid = np.isfinite(frequency) & np.isfinite(gv_norm) & (frequency > 0)
        ax.scatter(
            frequency[valid],
            gv_norm[valid],
            s=12,
            alpha=0.35,
            marker=marker,
            color=color,
            label=dataset["label"],
        )
    ax.set_xlabel("Frequency (THz)")
    ax.set_ylabel("Group velocity (km/s)")
    ax.legend(frameon=False)
    ax.grid(color="0.9", linewidth=0.9)


def draw_compare_relaxation_time(ax, datasets):
    """Draw overlaid relaxation-time data."""
    colors = ["tab:blue", "tab:orange"]
    markers = ["o", "x"]
    for dataset, color, marker in zip(datasets, colors, markers):
        transport = dataset["transport"]
        frequency = transport["frequency"]
        temp_index = transport["tau_temperature_index"]
        channels = tau_channels(transport["tau_mode"], transport["gamma"])
        for channel in channels:
            gamma = transport["gamma"][channel][temp_index]
            tau = np.full_like(gamma, np.nan, dtype=float)
            positive = gamma > 0
            tau[positive] = 1.0 / (4.0 * np.pi * gamma[positive])
            valid = np.isfinite(frequency) & np.isfinite(tau) & (frequency > 0)
            label = dataset["label"] if len(channels) == 1 else f"{dataset['label']} {channel}"
            ax.scatter(
                frequency[valid],
                tau[valid],
                s=12,
                alpha=0.35,
                marker=marker,
                color=color,
                label=label,
            )

    ax.set_yscale("log")
    ax.set_xlabel("Frequency (THz)")
    ax.set_ylabel("Relaxation time (ps)")
    ax.legend(frameon=False)
    ax.grid(color="0.9", linewidth=0.9)


def draw_compare_kappa(ax, datasets):
    """Draw overlaid thermal-conductivity components."""
    colors = ["tab:blue", "tab:orange"]
    linestyles = ["-", "--"]
    markers = ["o", "s"]
    component_markers = {"xx": "o", "yy": "s", "zz": "^", "average": "D"}

    for dataset, color, linestyle, marker in zip(datasets, colors, linestyles, markers):
        transport = dataset["transport"]
        temperature = transport["temperature"]
        kappa = transport["kappa"]
        kxx, kyy, kzz = kappa[:, 0], kappa[:, 1], kappa[:, 2]
        kavg = (kxx + kyy + kzz) / 3.0
        mode = transport["kappa_mode"]
        components = {
            "x": (kxx, "xx"),
            "y": (kyy, "yy"),
            "z": (kzz, "zz"),
        }
        if mode == "all":
            for values, component in [
                (kxx, "xx"),
                (kyy, "yy"),
                (kzz, "zz"),
                (kavg, "average"),
            ]:
                ax.plot(
                    temperature,
                    values,
                    color=color,
                    linestyle=linestyle,
                    marker=component_markers[component],
                    label=f"{dataset['label']} {component}",
                )
        else:
            values, component = components[mode]
            ax.plot(
                temperature,
                values,
                color=color,
                linestyle=linestyle,
                marker=marker,
                label=f"{dataset['label']} {component}",
            )

    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel(r"Thermal conductivity (W m$^{-1}$ K$^{-1}$)")
    ax.legend(frameon=False, ncol=2 if datasets[0]["transport"]["kappa_mode"] == "all" else 1)
    ax.grid(color="0.9", linewidth=0.9)
