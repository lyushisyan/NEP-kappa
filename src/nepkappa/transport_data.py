"""Shared interpretation of phono3py conductivity datasets."""

import warnings

import h5py
import numpy as np

EV_TO_J = 1.602176634e-19
ANGSTROM3_TO_M3 = 1.0e-30
AXIS_INDEX = {"x": 0, "y": 1, "z": 2}


def read_conductivity_components(handle):
    """Read total conductivity and optional verified Wigner components.

    RTA stores the total in ``kappa``. Phono3py 4.1 LBTE stores the particle
    solution there and separately writes ``kappa_intra_exact``. Accept either
    explicit convention, retaining the original HDF5 and rejecting mismatches.
    """
    kappa = np.asarray(handle['kappa'], dtype=float)
    components = {'kappa': kappa}
    wigner_keys = {'kappa_intra', 'kappa_intra_exact', 'kappa_inter'} & set(handle)
    if not wigner_keys:
        return components
    particle_key = 'kappa_intra_exact' if 'kappa_intra_exact' in handle else 'kappa_intra'
    if particle_key not in handle or 'kappa_inter' not in handle:
        raise ValueError('Incomplete Wigner components')
    intra = np.asarray(handle[particle_key], dtype=float)
    inter = np.asarray(handle['kappa_inter'], dtype=float)
    if intra.shape != kappa.shape or inter.shape != kappa.shape:
        raise ValueError('Wigner conductivity components have inconsistent shapes')
    if not all(np.isfinite(values).all() for values in (kappa, intra, inter)):
        raise ValueError('Nonfinite Wigner conductivity components')
    total = intra + inter
    stores_particle = particle_key == 'kappa_intra_exact' and np.allclose(
        kappa, intra, rtol=5e-4, atol=5e-4
    )
    if not stores_particle and not np.allclose(kappa, total, rtol=5e-4, atol=5e-4):
        raise ValueError('Wigner total conductivity differs from particle plus coherence')
    return {'kappa': total, 'kappa_intra': intra, 'kappa_inter': inter}


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
            **read_conductivity_components(handle),
            "heat_capacity": handle["heat_capacity"][:],
            "group_velocity": handle["group_velocity"][:],
            "frequency": handle["frequency"][:],
            "weight": handle["weight"][:],
            "mesh": handle["mesh"][:],
            "gamma": {},
        }
        if "gamma" in handle:
            transport["gamma"]["phph"] = handle["gamma"][:]
        isotope_gamma = handle["gamma_isotope"][:] if "gamma_isotope" in handle else None
        if "gamma_N" in handle:
            transport["gamma"]["normal"] = handle["gamma_N"][:]
        if "gamma_U" in handle:
            transport["gamma"]["umklapp"] = handle["gamma_U"][:]
        if "mode_kappa" in handle:
            transport["mode_kappa"] = handle["mode_kappa"][:]
        # Ordinary LBTE files also store the RTA solution from the same run.
        # Wigner components have separate semantics and retain their own plots.
        if "kappa_RTA" in handle and "kappa_intra" not in transport:
            transport["kappa_RTA"] = np.asarray(handle["kappa_RTA"], dtype=float)

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
    # phono3py stores isotope linewidths separately, normally without a T axis.
    # Use the saved result, not the current YAML switch, when plotting old data.
    if isotope_gamma is not None:
        if isotope_gamma.shape == frequency.shape:
            isotope_gamma = np.broadcast_to(isotope_gamma, shape)
        if isotope_gamma.shape != shape or not np.isfinite(isotope_gamma).all() or np.any(isotope_gamma < 0):
            raise ValueError("Invalid gamma_isotope: expected finite nonnegative q-by-branch or temperature-by-q-by-branch linewidths")
        transport["gamma"]["isotope"] = isotope_gamma
    if "phph" in transport["gamma"]:
        if isotope_gamma is None and getattr(config, "isotope", False) and getattr(config, "plot_tau", "total") in {"total", "all", "isotope"}:
            raise ValueError("Isotope scattering was requested but gamma_isotope is missing in kappa HDF5; use plot.tau: phph for phonon-only rates or supply a result containing gamma_isotope")
        transport["gamma"]["total"] = transport["gamma"]["phph"] + (isotope_gamma if isotope_gamma is not None else 0.0)
    kappa = transport["kappa"]
    if kappa.ndim != 2 or kappa.shape[0] != len(temperature) or kappa.shape[1] < 3:
        raise ValueError(f"Unsupported conductivity shape {kappa.shape}; select a single transport solution")
    if not np.isfinite(kappa).all():
        raise ValueError("Nonfinite conductivity values; results may be incomplete")
    if "kappa_RTA" in transport:
        rta = transport["kappa_RTA"]
        if rta.shape != kappa.shape or not np.isfinite(rta).all():
            raise ValueError("kappa_RTA must match the kappa tensor shape and contain finite values")

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
    elif "kappa_RTA" in transport:
        transport["kappa_RTA"] *= correction["factor"]
        transport["kappa_scenarios"] = [
            {"label": "LBTE", "temperature": temperature,
             "kappa": transport["kappa"], "linestyle": "-"},
            {"label": "RTA", "temperature": temperature,
             "kappa": transport["kappa_RTA"], "linestyle": "--"},
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
    transport["scattering_rate_min"] = getattr(config, "plot_scattering_rate_min", 0.0)
    transport["kappa_mode"] = getattr(config, "plot_kappa", "all")
    return transport


def mode_kappa_contributions(mode_kappa, mesh, geometry_factor=1.0):
    """Return per-mode conductivity whose sum equals the reported tensor.

    phono3py stores irreducible-grid weights inside ``mode_kappa`` and applies
    the full-mesh normalization only when forming ``kappa``.
    """
    return np.asarray(mode_kappa, dtype=float) / int(np.prod(mesh)) * float(
        geometry_factor
    )


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


def scattering_channels(tau_mode, gamma_data):
    """Return scattering-rate channels requested by YAML."""
    if tau_mode in {"total", "all"} and "phph" in gamma_data:
        channels = ["phph"]
        if "isotope" in gamma_data:
            channels.extend(["isotope", "total"])
        if tau_mode == "all":
            channels.extend(name for name in ("normal", "umklapp") if name in gamma_data)
        return channels
    if tau_mode == "all":
        return [name for name in ("total", "normal", "umklapp") if name in gamma_data]
    if tau_mode == "nu":
        missing = [name for name in ("normal", "umklapp") if name not in gamma_data]
        if missing:
            datasets = ", ".join(
                "gamma_N" if name == "normal" else "gamma_U" for name in missing
            )
            raise ValueError(
                f"Scattering-rate mode 'nu' requires {datasets} in kappa HDF5."
            )
        return ["normal", "umklapp"]
    if tau_mode not in gamma_data:
        raise ValueError(
            f"Scattering-rate channel '{tau_mode}' is not available in kappa HDF5."
        )
    return [tau_mode]
