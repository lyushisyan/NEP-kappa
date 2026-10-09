"""Plotting utilities for NEP-kappa results."""

from __future__ import annotations

import os
import tempfile
import warnings
from pathlib import Path

import yaml

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
from scipy.constants import Avogadro
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
    "scattering_rate",
    "kappa",
]
AXIS_INDEX = {"x": 0, "y": 1, "z": 2}
COMPARE_COLORS = [
    "black",
    "tab:blue",
    "tab:orange",
    "tab:green",
    "tab:red",
    "tab:purple",
    "tab:brown",
    "tab:pink",
    "tab:gray",
    "tab:olive",
    "tab:cyan",
]
COMPARE_LINESTYLES = ["-", "--", "-.", ":"]
COMPARE_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "<", ">", "*"]


def plot_results(config):
    """Create standard phonon and thermal-transport plots."""
    result_dir = Path(config.result_dir).resolve()
    plot_dir = result_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    disp_path = find_phonon_metadata(result_dir)
    fc2_path = result_dir / "fc2.hdf5"
    fourphonon_dir = result_dir / getattr(config, "fp_workdir", "fourphonon")
    has_fourphonon = (fourphonon_dir / "fourphonon-summary.yaml").is_file()
    # A FourPhonon run can share its directory with an older phono3py result
    # on a different mesh. Only use a matching HDF5 as a 3ph reference.
    if has_fourphonon:
        mesh = getattr(config, "mesh", None)
        name = f"kappa-m{int(mesh[0])}{int(mesh[1])}{int(mesh[2])}.hdf5" if mesh is not None else None
        kappa_path = result_dir / name if name else find_kappa_file(result_dir)
    else:
        kappa_path = find_kappa_file(result_dir, config.mesh)
    missing = [path for path in (disp_path, fc2_path) if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing file(s) required for plotting: "
            + ", ".join(str(path) for path in missing)
        )

    layout = getattr(config, "plot_layout", "separate")
    dpi = int(getattr(config, "plot_dpi", 300))

    print("\n[Step 4] Plot Results")
    print(f"  - Reading {disp_path}")
    print(f"  - Reading {fc2_path}")
    if kappa_path.exists():
        print(f"  - Reading {kappa_path}")
    elif has_fourphonon:
        print(f"  - Reading FourPhonon results from {fourphonon_dir}")
    else:
        print("  - No kappa file: plotting harmonic properties from FC2 only")
    print(f"  - Plot layout: {layout}")
    print(f"  - Band path: {getattr(config, 'plot_path', 'seekpath')}")
    phono3py_yaml = load_phono3py_yaml(disp_path)
    phonon = make_phonopy(phono3py_yaml, fc2_path)
    plot_data = build_plot_data(phonon, phono3py_yaml.unitcell, kappa_path, config)
    if has_fourphonon:
        attach_fourphonon_results(plot_data["transport"], fourphonon_dir, config)
    figures = available_figures([plot_data["transport"]], include_fourphonon=True)
    print(f"  - Plot figures: {', '.join(figures)}")
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
    """Create two-result or multi-model comparison plots."""
    configured_datasets = getattr(config, "datasets", None) or [
        {"directory": config.dft_dir, "label": config.reference_label},
        {"directory": config.nep_dir, "label": config.candidate_label},
    ]
    compare_dir = Path(config.compare_dir).resolve()
    plot_dir = compare_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    layout = getattr(config, "plot_layout", "separate")
    dpi = int(getattr(config, "plot_dpi", 300))

    print("\n[Compare] Plot Multiple Phonon/Transport Results")
    for dataset in configured_datasets:
        print(f"  - Dataset ({dataset['label']}): {Path(dataset['directory']).resolve()}")
    print(f"  - Compare directory: {compare_dir}")
    print(f"  - Plot layout: {layout}")
    print(f"  - Band path: {getattr(config, 'plot_path', 'seekpath')}")
    datasets = [
        load_result_plot_data(
            Path(dataset["directory"]).resolve(), config, dataset["label"]
        )
        for dataset in configured_datasets
    ]
    figures = available_figures([dataset["transport"] for dataset in datasets])
    print(f"  - Plot figures: {', '.join(figures)}")
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
        if mesh is not None:
            raise ValueError(f"No conductivity file matches mesh {list(mesh)}. Available: "
                             + ", ".join(p.name for p in candidates))
        if len(candidates) != 1:
            raise ValueError("Multiple conductivity files found; select an explicit mesh")
        return candidates[0]
    if mesh is None:
        return result_dir / "kappa-m*.hdf5"
    return expected


def find_phonon_metadata(result_dir):
    """FC2 still requires its matching cell/supercell metadata."""
    for name in ("phono3py_disp.yaml", "phonopy.yaml", "phonopy_disp.yaml"):
        path = Path(result_dir) / name
        if path.is_file():
            return path
    return Path(result_dir) / "phono3py_disp.yaml"


def load_phono3py_yaml(path):
    """Load phono3py_disp.yaml."""
    ph3yml = Phono3pyYaml()
    ph3yml.read(str(path))
    return ph3yml


def make_phonopy(ph3yml, fc2_path):
    """Create a Phonopy object from phono3py metadata and fc2.hdf5."""
    matrix = getattr(ph3yml, "phonon_supercell_matrix", None)
    if matrix is None:
        matrix = ph3yml.supercell_matrix
    phonon = Phonopy(
        ph3yml.unitcell,
        supercell_matrix=matrix,
        primitive_matrix=ph3yml.primitive_matrix,
    )
    phonon.force_constants = read_force_constants_hdf5(str(fc2_path))
    if getattr(ph3yml, "nac_params", None) is not None:
        phonon.nac_params = ph3yml.nac_params
    return phonon


def build_plot_data(phonon, unitcell, kappa_path, config):
    """Build all data needed by the plotting functions."""
    transport = (
        read_transport_data(kappa_path, phonon.primitive.volume,
                            phonon.primitive.cell, config)
        if kappa_path is not None and Path(kappa_path).exists()
        else build_harmonic_properties(phonon, config)
    )
    data = {
        "band": build_band_data(phonon, unitcell, config),
        "dos": build_dos_data(phonon, transport["mesh"]),
        "transport": transport,
    }
    return data


def load_result_plot_data(result_dir, config, label):
    """Load one completed result directory for plotting or comparison."""
    disp_path = find_phonon_metadata(result_dir)
    fc2_path = result_dir / "fc2.hdf5"
    kappa_path = find_kappa_file(result_dir, getattr(config, "mesh", None))
    missing = [path for path in (disp_path, fc2_path) if not path.exists()]
    if missing:
        raise FileNotFoundError(
            f"Missing file(s) required for {label}: "
            + ", ".join(str(path) for path in missing)
        )

    print(f"  - Reading {label}: {disp_path}")
    print(f"  - Reading {label}: {fc2_path}")
    print(f"  - {label}: {'transport data' if kappa_path.exists() else 'harmonic-only data'}")

    phono3py_yaml = load_phono3py_yaml(disp_path)
    phonon = make_phonopy(phono3py_yaml, fc2_path)
    data = build_plot_data(phonon, phono3py_yaml.unitcell, kappa_path, config)
    data["label"] = label
    return data


def seekpath_structure(unitcell):
    """Return the tuple format expected by seekpath."""
    return (unitcell.cell, unitcell.scaled_positions, unitcell.numbers)


def make_band_paths(unitcell, config, points_per_segment=51):
    """Build high-symmetry q-point paths with seekpath."""
    if getattr(config, "plot_path", "seekpath") == "custom":
        return make_custom_band_paths(config, points_per_segment)

    # Frequencies are evaluated in this actual primitive basis, not seekpath's
    # independently standardized primitive basis.
    sp_path = seekpath.get_path_orig_cell(seekpath_structure(unitcell))
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
        if not np.isfinite(start).all() or not np.isfinite(end).all():
            raise ValueError("Custom q-points must be finite")
        if np.allclose(start, end, rtol=0, atol=1e-12):
            raise ValueError("Custom path segments must have distinct endpoints")
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
    paths, labels = make_band_paths(phonon.primitive, config)
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
        required = {"temperature", "kappa", "heat_capacity", "group_velocity",
                    "frequency", "weight", "mesh"}
        missing = required - set(handle.keys())
        if missing:
            raise ValueError(
                f"Unsupported or incomplete transport file {kappa_path}: missing "
                + ", ".join(sorted(missing))
                + ". Wigner-only datasets and split grid-point files are not standard kappa inputs."
            )
        temperature = handle["temperature"][:]
        transport = {
            "temperature": temperature,
            "kappa": handle["kappa"][:],
            "heat_capacity": handle["heat_capacity"][:],
            "group_velocity": handle["group_velocity"][:],
            "frequency": handle["frequency"][:],
            "weight": handle["weight"][:],
            "mesh": handle["mesh"][:],
            "gamma": {},
        }
        if "gamma" in handle:
            transport["gamma"]["total"] = handle["gamma"][:]
        if "gamma_N" in handle:
            transport["gamma"]["normal"] = handle["gamma_N"][:]
        if "gamma_U" in handle:
            transport["gamma"]["umklapp"] = handle["gamma_U"][:]
        if "mode_kappa" in handle:
            transport["mode_kappa"] = handle["mode_kappa"][:]
        wigner_keys = {"kappa_intra", "kappa_inter"} & set(handle.keys())
        if wigner_keys:
            if wigner_keys != {"kappa_intra", "kappa_inter"}:
                raise ValueError(f"Incomplete Wigner components in {kappa_path}")
            transport["kappa_intra"] = handle["kappa_intra"][:]
            transport["kappa_inter"] = handle["kappa_inter"][:]

    if temperature.ndim != 1 or not len(temperature) or not np.isfinite(temperature).all():
        raise ValueError(f"Invalid temperature array in {kappa_path}")
    frequency = transport["frequency"]
    mesh = transport["mesh"]
    if frequency.ndim != 2 or not np.isfinite(frequency).all():
        raise ValueError("Transport frequencies must be a finite q-by-branch array")
    if mesh.shape != (3,) or np.any(mesh <= 0) or not np.isfinite(mesh).all():
        raise ValueError("Transport mesh must contain three positive dimensions")
    if (transport["heat_capacity"].shape != (len(temperature),)+frequency.shape or
            transport["group_velocity"].shape != frequency.shape+(3,) or
            transport["weight"].shape != (len(frequency),)):
        raise ValueError("Transport heat capacity, group velocity or weights have inconsistent shapes")
    shape = (len(temperature),) + transport["frequency"].shape
    for name, gamma in transport["gamma"].items():
        if gamma.shape != shape:
            raise ValueError(f"Unsupported {name} linewidth shape {gamma.shape}; expected {shape}")
    kappa = transport["kappa"]
    if kappa.ndim != 2 or kappa.shape[0] != len(temperature) or kappa.shape[1] < 3:
        raise ValueError(f"Unsupported conductivity shape {kappa.shape}; select a single transport solution")
    if not np.isfinite(kappa).all():
        raise ValueError("Nonfinite conductivity values; results may be incomplete")
    if "kappa_intra" in transport:
        intra, inter = transport["kappa_intra"], transport["kappa_inter"]
        if intra.shape != kappa.shape or inter.shape != kappa.shape:
            raise ValueError("Wigner conductivity components have inconsistent shapes")
        if not np.isfinite(intra).all() or not np.isfinite(inter).all():
            raise ValueError("Nonfinite Wigner conductivity components")
        if not np.allclose(kappa, intra + inter, rtol=5e-4, atol=5e-4):
            raise ValueError("Wigner total conductivity differs from particle plus coherence")

    correction = effective_geometry_correction(config, primitive_cell, primitive_volume)
    effective_volume = primitive_volume / correction["factor"]
    transport["kappa"] = transport["kappa"] * correction["factor"]
    if "kappa_intra" in transport:
        transport["kappa_scenarios"] = [
            {"label": "particle", "temperature": temperature,
             "kappa": transport["kappa_intra"] * correction["factor"]},
            {"label": "coherence", "temperature": temperature,
             "kappa": transport["kappa_inter"] * correction["factor"]},
            {"label": "total", "temperature": temperature,
             "kappa": transport["kappa"]},
        ]
    if "mode_kappa" in transport:
        transport["mode_kappa"] = mode_kappa_contributions(
            transport["mode_kappa"], transport["mesh"], correction["factor"]
        )
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
    selected = float(temperature[transport["tau_temperature_index"]])
    requested = float(getattr(config, "plot_temperature", 300.0))
    if not np.isclose(selected, requested):
        warnings.warn(f"Requested plot temperature {requested:g} K is unavailable; using {selected:g} K", RuntimeWarning)
    transport["tau_mode"] = getattr(config, "plot_tau", "total")
    transport["kappa_mode"] = getattr(config, "plot_kappa", "all")
    return transport


def read_fourphonon_kappa(path, geometry_factor):
    """Read a normalized FourPhonon T/kxx/kyy/kzz table."""
    if not path.is_file():
        raise FileNotFoundError(f"Missing FourPhonon conductivity table: {path}")
    data = np.atleast_2d(np.loadtxt(path, comments="#"))
    if data.shape[1] != 4 or not np.isfinite(data).all():
        raise ValueError(f"Invalid FourPhonon conductivity table: {path}")
    return {"temperature": data[:, 0], "kappa": data[:, 1:4] * geometry_factor}


def read_fourphonon_rate(path):
    """Read BTE.w_3ph/BTE.w_4ph; FourPhonon writes angular frequency."""
    data = np.atleast_2d(np.loadtxt(path, comments="#"))
    if data.shape[1] != 2 or not np.isfinite(data).all():
        raise ValueError(f"Invalid FourPhonon scattering-rate file: {path}")
    return np.column_stack((data[:, 0] / (2.0 * np.pi), data[:, 1]))


def read_fourphonon_nu(path):
    """Read optional state/frequency/N/U/total output from an NU-capable build."""
    data = np.atleast_2d(np.loadtxt(path, comments="#"))
    if data.shape[1] != 5 or not np.isfinite(data).all():
        raise ValueError(f"Invalid FourPhonon N/U scattering-rate file: {path}")
    return np.column_stack((data[:, 1] / (2.0 * np.pi), data[:, 2:4]))


def attach_fourphonon_results(transport, workdir, config):
    """Add only completed FourPhonon outputs to the standard plotting data."""
    summary_path = workdir / "fourphonon-summary.yaml"
    summary = yaml.safe_load(summary_path.read_text(encoding="utf-8")) or {}
    solutions = summary.get("solutions") or {}
    factor = transport["geometry_correction"]["factor"]
    scenarios = list(transport.get("kappa_scenarios", []))
    if "kappa" in transport and not scenarios:
        scenarios.append({"label": "3ph only", "temperature": transport["temperature"],
                          "kappa": transport["kappa"]})
    labels = {
        "rta": "3ph+4ph RTA",
        "3ph-iterative": "3ph LBTE + 4ph RTA",
        "full-iterative": "3ph+4ph LBTE",
    }
    solver = getattr(config, "fp_solver", "rta")
    for solution in ("rta", "iterative"):
        if solution not in solutions:
            continue
        result = read_fourphonon_kappa(workdir / f"kappa4-{solution}.dat", factor)
        if solution == "rta" and summary.get("wigner", {}).get("components_verified"):
            coherence = read_fourphonon_kappa(workdir / "kappa4-rta-coherence.dat", factor)
            total = read_fourphonon_kappa(workdir / "kappa4-rta-total.dat", factor)
            if not np.array_equal(result["temperature"], coherence["temperature"]) or not np.array_equal(result["temperature"], total["temperature"]):
                raise ValueError("FourPhonon Wigner temperatures do not match")
            if not np.allclose(total["kappa"], result["kappa"] + coherence["kappa"], rtol=5e-4, atol=5e-4):
                raise ValueError("FourPhonon Wigner total differs from particle plus coherence")
            for label, item in (("particle", result), ("coherence", coherence), ("total", total)):
                scenarios.append({"label": f"3ph+4ph {label}", **item})
            primary_data = total
        else:
            label = labels["rta"] if solution == "rta" else labels.get(solver, "3ph+4ph iterative")
            scenarios.append({"label": label, **result})
            primary_data = result
        if solution == summary.get("primary"):
            transport["fourphonon_primary"] = primary_data
    if scenarios:
        transport["kappa_scenarios"] = scenarios
    primary = transport.get("fourphonon_primary")
    if primary is not None:
        transport["heat_capacity_temperature"] = transport["temperature"]
        transport["temperature"] = primary["temperature"]
        transport["kappa"] = primary["kappa"]
    requested = float(getattr(config, "plot_temperature", 300.0))
    dirs = []
    for directory in workdir.glob("T*K"):
        try:
            temperature = float(directory.name[1:-1])
        except ValueError:
            continue
        dirs.append((abs(temperature - requested), temperature, directory))
    if dirs:
        _, selected, directory = min(dirs, key=lambda item: item[0])
        if not np.isclose(selected, requested):
            warnings.warn(f"Requested FourPhonon plot temperature {requested:g} K is unavailable; using {selected:g} K", RuntimeWarning)
        rates = {}
        for channel in ("3ph", "4ph"):
            path = directory / f"BTE.w_{channel}"
            if path.is_file():
                rates[channel] = read_fourphonon_rate(path)
        transport["fourphonon_rates"] = rates
        transport["fourphonon_rate_temperature"] = selected
        for source in ("3ph4ph", "3ph"):
            path = directory / f"BTE.w_{source}_NU"
            if path.is_file():
                transport["fourphonon_nu"] = read_fourphonon_nu(path)
                transport["fourphonon_nu_source"] = source
                break
    transport["kappa_mode"] = getattr(config, "plot_kappa", "all")


def build_harmonic_properties(phonon, config):
    """Compute Cv and group velocities from FC2; never invent scattering data."""
    configured_mesh = getattr(config, "mesh", None)
    mesh = np.asarray([21, 21, 21] if configured_mesh is None else configured_mesh, dtype=int)
    if mesh.shape != (3,) or np.any(mesh <= 0):
        raise ValueError("Harmonic plotting requires three positive mesh dimensions")
    spec = np.asarray(getattr(config, "temps", [100, 1000, 100]), dtype=float)
    if spec.ndim != 1 or not np.isfinite(spec).all():
        raise ValueError("Heat-capacity temperature specification must be finite")
    if spec.size == 1:
        temperatures = spec
    elif spec.size == 3 and spec[2] > 0 and spec[1] >= spec[0]:
        temperatures = np.arange(spec[0], spec[1] + spec[2]*1e-8, spec[2])
    else:
        raise ValueError("Heat-capacity temperatures must be [T] or [start, stop, step]")
    if not np.isfinite(temperatures).all() or np.any(temperatures < 0):
        raise ValueError("Heat-capacity temperatures must be finite and nonnegative")
    phonon.run_mesh(mesh, is_gamma_center=True, with_group_velocities=True)
    modes = phonon.get_mesh_dict()
    frequencies = np.asarray(modes["frequencies"])
    if np.any(frequencies < -1e-4):
        warnings.warn("Imaginary modes present: Cv excludes nonpositive modes and is not a stable-phase thermodynamic prediction", RuntimeWarning)
    phonon.run_thermal_properties(temperatures=temperatures, cutoff_frequency=0.)
    thermal = phonon.get_thermal_properties_dict()
    correction = effective_geometry_correction(config, phonon.primitive.cell,
                                               phonon.primitive.volume)
    # Phonopy returns J/(mol primitive cells K), not eV/K per mode.
    cv = np.asarray(thermal["heat_capacity"])/Avogadro
    cv /= phonon.primitive.volume*ANGSTROM3_TO_M3/correction["factor"]
    return {"temperature":np.asarray(thermal["temperatures"]),
            "volume_heat_capacity":cv, "frequency":frequencies,
            "group_velocity":np.asarray(modes["group_velocities"]),
            "weight":np.asarray(modes["weights"]), "mesh":mesh,
            "geometry_correction":correction, "source":"FC2 harmonic"}


def mode_kappa_contributions(mode_kappa, mesh, geometry_factor=1.0):
    """Return per-mode conductivity whose sum equals the reported tensor.

    phono3py stores irreducible-grid weights inside ``mode_kappa`` and applies
    the full-mesh normalization only when forming ``kappa``.
    """
    return np.asarray(mode_kappa, dtype=float) / int(np.prod(mesh)) * float(
        geometry_factor
    )


def available_figures(transports, include_fourphonon=False):
    """Return standard figures plus analyses supported by every dataset."""
    required = {"heat_capacity": {"temperature", "volume_heat_capacity"},
                "group_velocity": {"frequency", "group_velocity"},
                "relaxation_time": {"frequency", "gamma", "tau_temperature_index"},
                "scattering_rate": {"frequency", "gamma", "tau_temperature_index"},
                "kappa": {"temperature", "kappa"}}
    figures = [name for name in DEFAULT_FIGURES if name not in required or
               (transports and all(required[name] <= set(t) for t in transports))]
    for name in ("relaxation_time", "scattering_rate"):
        if name in figures and not all(t.get("gamma") for t in transports):
            figures.remove(name)
    if transports and all("mode_kappa" in transport for transport in transports):
        figures.append("cumulative_kappa")
    if include_fourphonon and len(transports) == 1:
        rates = transports[0].get("fourphonon_rates", {})
        if "fourphonon_nu" in transports[0]:
            figures.insert(figures.index("kappa") if "kappa" in figures else len(figures),
                           "scattering_rate_nu")
        for channel in ("3ph", "4ph"):
            if channel in rates:
                figures.insert(figures.index("kappa") if "kappa" in figures else len(figures),
                               f"scattering_rate_{channel}")
    return figures


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
        "scattering_rate": (5.8, 4.6),
        "scattering_rate_3ph": (5.8, 4.6),
        "scattering_rate_4ph": (5.8, 4.6),
        "scattering_rate_nu": (5.8, 4.6),
        "cumulative_kappa": (5.8, 4.6),
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
    ncols = 2 if len(figures) == 4 else min(len(figures), 3)
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
        "scattering_rate": (5.8, 4.6),
        "cumulative_kappa": (5.8, 4.6),
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
    ncols = 2 if len(figures) == 4 else min(len(figures), 3)
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
    elif name == "scattering_rate":
        draw_scattering_rate(ax, plot_data["transport"])
    elif name == "scattering_rate_3ph":
        draw_fourphonon_scattering_rate(ax, plot_data["transport"], "3ph")
    elif name == "scattering_rate_4ph":
        draw_fourphonon_scattering_rate(ax, plot_data["transport"], "4ph")
    elif name == "scattering_rate_nu":
        draw_fourphonon_nu(ax, plot_data["transport"])
    elif name == "cumulative_kappa":
        draw_cumulative_kappa(ax, plot_data["transport"])
    elif name == "kappa":
        draw_kappa(ax, plot_data["transport"])
    else:
        raise ValueError(f"Unknown plot figure: {name}")


def draw_compare_figure(name, ax, datasets):
    """Draw one named multi-model comparison figure."""
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
    elif name == "scattering_rate":
        draw_compare_scattering_rate(ax, datasets)
    elif name == "cumulative_kappa":
        draw_compare_cumulative_kappa(ax, datasets)
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
        transport.get("heat_capacity_temperature", transport["temperature"]),
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


def draw_scattering_rate(ax, transport):
    """Draw total, Normal, and/or Umklapp scattering rates in ps^-1."""
    frequency = transport["frequency"]
    temp_index = transport["tau_temperature_index"]
    channels = tau_channels(transport["tau_mode"], transport["gamma"])
    colors = {
        "total": "tab:orange",
        "normal": "tab:blue",
        "umklapp": "tab:green",
    }
    labels = {"total": "total", "normal": "N", "umklapp": "U"}
    for name in channels:
        rate = 4.0 * np.pi * transport["gamma"][name][temp_index]
        valid = np.isfinite(frequency) & np.isfinite(rate) & (frequency > 0) & (rate > 0)
        ax.scatter(
            frequency[valid],
            rate[valid],
            s=12,
            alpha=0.35,
            color=colors[name],
            label=labels[name],
        )
    ax.set_yscale("log")
    ax.set_xlabel("Frequency (THz)")
    ax.set_ylabel(r"Scattering rate (ps$^{-1}$)")
    if len(channels) > 1:
        ax.legend(frameon=False)
    ax.grid(color="0.9", linewidth=0.9)


def draw_fourphonon_scattering_rate(ax, transport, channel):
    """Draw one FourPhonon channel without combining 3ph and 4ph rates."""
    data = transport["fourphonon_rates"][channel]
    valid = (data[:, 0] > 0) & (data[:, 1] > 0)
    color = "tab:blue" if channel == "3ph" else "tab:orange"
    ax.scatter(data[valid, 0], data[valid, 1], s=12, alpha=0.4, color=color)
    ax.set_yscale("log")
    ax.set_xlabel("Frequency (THz)")
    ax.set_ylabel(rf"{channel} scattering rate (ps$^{{-1}}$)")
    ax.grid(color="0.9", linewidth=0.9)


def draw_fourphonon_nu(ax, transport):
    """Draw optional FourPhonon N and U rates together."""
    data = transport["fourphonon_nu"]
    for column, label, color in ((1, "N", "tab:blue"), (2, "U", "tab:green")):
        valid = (data[:, 0] > 0) & (data[:, column] > 0)
        ax.scatter(data[valid, 0], data[valid, column], s=12, alpha=0.4,
                   color=color, label=label)
    ax.set_yscale("log")
    ax.set_xlabel("Frequency (THz)")
    ax.set_ylabel(r"Scattering rate (ps$^{-1}$)")
    ax.legend(frameon=False)
    ax.grid(color="0.9", linewidth=0.9)


def tau_channels(tau_mode, gamma_data):
    """Return relaxation-time channels requested by YAML."""
    if tau_mode == "all":
        return [name for name in ("total", "normal", "umklapp") if name in gamma_data]
    if tau_mode == "nu":
        missing = [name for name in ("normal", "umklapp") if name not in gamma_data]
        if missing:
            datasets = ", ".join(
                "gamma_N" if name == "normal" else "gamma_U" for name in missing
            )
            raise ValueError(
                f"Relaxation-time mode 'nu' requires {datasets} in kappa HDF5."
            )
        return ["normal", "umklapp"]
    if tau_mode not in gamma_data:
        raise ValueError(
            f"Relaxation-time channel '{tau_mode}' is not available in kappa HDF5."
        )
    return [tau_mode]


def draw_kappa(ax, transport):
    """Draw selected thermal conductivity tensor diagonal components."""
    if transport.get("kappa_scenarios"):
        mode = transport["kappa_mode"]
        for scenario in transport["kappa_scenarios"]:
            values = scenario["kappa"][:, :3]
            ordinate = values.mean(axis=1) if mode == "all" else values[:, AXIS_INDEX[mode]]
            ax.plot(scenario["temperature"], ordinate, marker="o", label=scenario["label"])
        ax.set_xlabel("Temperature (K)")
        axis_label = "average" if mode == "all" else mode + mode
        ax.set_ylabel(rf"Thermal conductivity {axis_label} (W m$^{{-1}}$ K$^{{-1}}$)")
        ax.legend(frameon=False)
        ax.grid(color="0.9", linewidth=0.9)
        return
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


def cumulative_kappa_curves(transport):
    """Return frequency-sorted cumulative conductivity curves at plot T."""
    frequency = np.asarray(transport["frequency"], dtype=float).ravel()
    temp_index = transport["tau_temperature_index"]
    mode_kappa = np.asarray(transport["mode_kappa"][temp_index], dtype=float)
    values = {
        "x": mode_kappa[..., 0].ravel(),
        "y": mode_kappa[..., 1].ravel(),
        "z": mode_kappa[..., 2].ravel(),
    }
    values["average"] = (values["x"] + values["y"] + values["z"]) / 3.0
    valid = np.isfinite(frequency) & (frequency > 0)
    order = np.argsort(frequency[valid])
    frequencies = frequency[valid][order]
    return {
        name: (frequencies, np.cumsum(component[valid][order]))
        for name, component in values.items()
    }


def draw_cumulative_kappa(ax, transport):
    """Draw cumulative conductivity against phonon frequency."""
    curves = cumulative_kappa_curves(transport)
    mode = transport["kappa_mode"]
    selected = ["x", "y", "z", "average"] if mode == "all" else [mode]
    labels = {"x": "xx", "y": "yy", "z": "zz", "average": "average"}
    for name in selected:
        frequency, cumulative = curves[name]
        ax.plot(frequency, cumulative, label=labels[name])
    ax.set_xlabel("Frequency (THz)")
    ax.set_ylabel(r"Cumulative thermal conductivity (W m$^{-1}$ K$^{-1}$)")
    ax.legend(frameon=False)
    ax.grid(color="0.9", linewidth=0.9)


def draw_compare_dispersion(ax, datasets):
    """Draw overlaid phonon dispersions."""
    reference_band = datasets[0]["band"]
    set_band_ticks(ax, reference_band)

    for index, dataset in enumerate(datasets):
        color, linestyle, _ = comparison_style(index)
        band = dataset["band"]
        first_line = True
        for segment_index, (distances, frequencies) in enumerate(
            zip(band["distances"], band["frequencies"])
        ):
            # Use the reference dataset's path coordinate so that small
            # lattice-constant differences do not shift high-symmetry points.
            reference_distances = np.asarray(
                reference_band["distances"][segment_index]
            )
            distances = np.asarray(distances)
            plot_distances = (
                reference_distances
                if reference_distances.shape == distances.shape
                else distances
            )
            lines = ax.plot(
                plot_distances,
                np.asarray(frequencies),
                color=color,
                linestyle=linestyle,
                linewidth=1.5,
                label="_nolegend_",
            )
            # A two-dimensional frequency array makes Matplotlib create one
            # Line2D per phonon branch.  Passing the dataset label to plot()
            # assigns it to every branch, so attach it only to the first line.
            if first_line and lines:
                lines[0].set_label(dataset["label"])
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
    for index, dataset in enumerate(datasets):
        color, linestyle, _ = comparison_style(index)
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
    for index, dataset in enumerate(datasets):
        color, _, marker = comparison_style(index)
        transport = dataset["transport"]
        ax.plot(
            transport.get("heat_capacity_temperature", transport["temperature"]),
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
    for index, dataset in enumerate(datasets):
        color, _, marker = comparison_style(index)
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
    for index, dataset in enumerate(datasets):
        color, _, marker = comparison_style(index)
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


def draw_compare_scattering_rate(ax, datasets):
    """Draw overlaid scattering rates for every configured dataset."""
    for index, dataset in enumerate(datasets):
        color, _, marker = comparison_style(index)
        transport = dataset["transport"]
        frequency = transport["frequency"]
        temp_index = transport["tau_temperature_index"]
        channels = tau_channels(transport["tau_mode"], transport["gamma"])
        for channel in channels:
            rate = 4.0 * np.pi * transport["gamma"][channel][temp_index]
            valid = (
                np.isfinite(frequency)
                & np.isfinite(rate)
                & (frequency > 0)
                & (rate > 0)
            )
            label = (
                dataset["label"]
                if len(channels) == 1
                else f"{dataset['label']} {channel}"
            )
            ax.scatter(
                frequency[valid],
                rate[valid],
                s=12,
                alpha=0.35,
                marker=marker,
                color=color,
                label=label,
            )
    ax.set_yscale("log")
    ax.set_xlabel("Frequency (THz)")
    ax.set_ylabel(r"Scattering rate (ps$^{-1}$)")
    ax.legend(frameon=False)
    ax.grid(color="0.9", linewidth=0.9)


def draw_compare_kappa(ax, datasets):
    """Draw overlaid thermal-conductivity components."""
    component_markers = {"xx": "o", "yy": "s", "zz": "^", "average": "D"}

    for index, dataset in enumerate(datasets):
        color, linestyle, marker = comparison_style(index)
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


def draw_compare_cumulative_kappa(ax, datasets):
    """Draw cumulative conductivity for every configured dataset."""
    for index, dataset in enumerate(datasets):
        color, linestyle, _ = comparison_style(index)
        transport = dataset["transport"]
        curves = cumulative_kappa_curves(transport)
        mode = transport["kappa_mode"]
        selected = ["x", "y", "z", "average"] if mode == "all" else [mode]
        labels = {"x": "xx", "y": "yy", "z": "zz", "average": "average"}
        for name in selected:
            frequency, cumulative = curves[name]
            ax.plot(
                frequency,
                cumulative,
                color=color,
                linestyle=linestyle,
                label=f"{dataset['label']} {labels[name]}",
            )
    ax.set_xlabel("Frequency (THz)")
    ax.set_ylabel(r"Cumulative thermal conductivity (W m$^{-1}$ K$^{-1}$)")
    ax.legend(frameon=False, ncol=2)
    ax.grid(color="0.9", linewidth=0.9)


def comparison_style(index):
    """Return a stable color, line style, and marker for any dataset count."""
    return (
        COMPARE_COLORS[index % len(COMPARE_COLORS)],
        COMPARE_LINESTYLES[index % len(COMPARE_LINESTYLES)],
        COMPARE_MARKERS[index % len(COMPARE_MARKERS)],
    )
