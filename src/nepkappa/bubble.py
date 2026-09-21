"""Diagonal on-shell cubic bubble corrections on an SSCHA harmonic reference.

This is a fixed-input-FC3 approximation, not the ensemble-averaged SSCHA
vertex, the free-energy Hessian, or a frequency-dependent spectral calculation.
It never changes FC2, eigenvectors, transport linewidths, or conductivity files.
"""

from __future__ import annotations

import copy
from pathlib import Path

import h5py
import numpy as np
import phono3py
import phonopy
from phonopy.interface.phonopy_yaml import PhonopyYaml
from phonopy.phonon.grid import get_ir_grid_points
import yaml

from nepkappa.provenance import data_sha256, file_identity
from nepkappa.sscha import temperature_directory_name, temperature_points


SCHEMA_VERSION = 1
# A conservative preflight bound for the dense interaction array alone.
# This is not a guarantee on total process RSS. Large cells need a tiled solver.
MAX_INTERACTION_BYTES = 2 * 1024**3


def bubble_scope():
    """Describe only the correction calculated by this module."""
    return {
        "label": "SSCHA + input-FC3 diagonal on-shell bubble (real part)",
        "dynamic_bubble_calculated": True,
        "bubble_evaluation": "real self-energy at the auxiliary band frequencies",
        "cubic_vertex": "input FC3, not an SSCHA ensemble-averaged third derivative",
        "free_energy_hessian_calculated": False,
        "sscha_cell_optimization": False,
        "bubble_corrected_transport": False,
        "limitations": [
            "Diagonal, one-shot on-shell approximation; no off-diagonal mode mixing.",
            "No frequency-dependent Dyson solution or spectral function is calculated.",
            "No imaginary self-energy is added to existing transport linewidths.",
            "Auxiliary FC2 and existing thermal conductivities are unchanged.",
            "Converge q mesh and principal-value epsilon together; defaults are pilot settings.",
            "Large shifts or nonpositive corrected squared frequencies invalidate a simple quasiparticle interpretation.",
        ],
    }


def on_shell_frequencies(frequencies, deltas, cutoff):
    """Return first-order and on-shell-Dyson estimates in ordinary THz.

    Phono3py returns Delta in ordinary THz, not angular units. The estimates
    are nu + Delta and sqrt(nu**2 + 2*nu*Delta); neither solves Delta(Omega).
    Acoustic/unstable reference modes and nonpositive radicands are marked,
    not clipped to apparently stable modes.
    """
    nu = np.asarray(frequencies, dtype=float)
    delta = np.asarray(deltas, dtype=float)
    if nu.ndim != 1 or delta.ndim != 2 or delta.shape[1] != len(nu):
        raise ValueError("Expected frequency (band,) and delta (epsilon, band).")
    if not np.isfinite(nu).all() or not np.isfinite(delta).all():
        raise ValueError("Non-finite frequencies or bubble shifts.")
    reference_valid = np.broadcast_to(nu > cutoff, delta.shape)
    squared = nu[None, :] ** 2 + 2 * nu[None, :] * delta
    valid = reference_valid & (squared > 0)
    corrected = np.full_like(delta, np.nan)
    corrected[valid] = np.sqrt(squared[valid])
    linear = np.where(reference_valid, nu[None, :] + delta, np.nan)
    ratio = np.full_like(delta, np.nan)
    np.divide(abs(delta), nu[None, :], out=ratio, where=reference_valid)
    return {
        "frequency_linear": linear,
        "frequency_squared_on_shell": squared,
        "frequency_on_shell": corrected,
        "valid_mode": valid,
        "relative_shift": ratio,
        "reference_mode_included": reference_valid,
    }


def _assert_same_cell(first, second, label):
    """Reject mismatched volume, atomic order, masses, or primitive convention."""
    if len(first) != len(second) or list(first.symbols) != list(second.symbols):
        raise ValueError(f"SSCHA/FC3 {label} species or atom ordering mismatch.")
    if not np.allclose(first.cell, second.cell, rtol=0, atol=1e-7):
        raise ValueError(f"SSCHA/FC3 {label} lattice mismatch; use volume-matched FC3.")
    diff = first.scaled_positions - second.scaled_positions
    if not np.allclose(diff - np.rint(diff), 0, rtol=0, atol=1e-7):
        raise ValueError(f"SSCHA/FC3 {label} positions or atom ordering mismatch.")
    if not np.allclose(first.masses, second.masses, rtol=0, atol=1e-7):
        raise ValueError(f"SSCHA/FC3 {label} masses mismatch.")


def _input_paths(config, temperature_dir):
    root = Path(config.result_dir).expanduser().resolve()
    settings = config.sections.scph
    return {
        "sscha_structure": temperature_dir / "phonopy_sscha.yaml",
        "sscha_summary": temperature_dir / "summary.yaml",
        "fc2": temperature_dir / "fc2.hdf5",
        "fc3": (Path(settings.transport_fc3).expanduser().resolve()
                if settings.transport_fc3 else root / "fc3.hdf5"),
        "metadata": (Path(settings.transport_metadata).expanduser().resolve()
                     if settings.transport_metadata else root / "phono3py_disp.yaml"),
    }


def _load_reference(paths, settings):
    # Disable discovery of BORN/FC files from the launch directory. NAC is taken
    # exclusively from the saved SSCHA structure, with the same primitive basis.
    sscha_yaml = PhonopyYaml().read(paths["sscha_structure"])
    sscha = phonopy.load(
        paths["sscha_structure"], force_constants_filename=paths["fc2"],
        produce_fc=False, is_nac=False, symmetrize_fc=False,
    )
    ph = phono3py.load(
        paths["metadata"], fc2_filename=paths["fc2"], fc3_filename=paths["fc3"],
        produce_fc=False, is_nac=False, symmetrize_fc=False, lang="C",
    )
    for label, first, second in (
        ("unit cell", sscha.unitcell, ph.unitcell),
        ("primitive", sscha.primitive, ph.phonon_primitive),
        ("FC2 supercell", sscha.supercell, ph.phonon_supercell),
        ("FC2/FC3 primitive", ph.phonon_primitive, ph.primitive),
    ):
        _assert_same_cell(first, second, label)
    if not np.isclose(sscha.unit_conversion_factor, ph.unit_conversion_factor):
        raise ValueError("SSCHA/FC3 frequency-unit conversion mismatch.")
    if ph.fc2 is None or ph.fc3 is None:
        raise ValueError("bubble requires existing FC2 and FC3; no force fitting is performed.")
    # load() does not expose the interaction cutoff. Reconstruct through the
    # public constructor instead of modifying a private backend attribute.
    loaded = ph
    ph = phono3py.Phono3py(
        loaded.unitcell, supercell_matrix=loaded.supercell_matrix,
        primitive_matrix=loaded.primitive_matrix,
        phonon_supercell_matrix=loaded.phonon_supercell_matrix,
        calculator=loaded.calculator, frequency_factor_to_THz=loaded.unit_conversion_factor,
        cutoff_frequency=settings.cutoff_frequency, lang="C",
    )
    ph.fc2, ph.fc3 = loaded.fc2, loaded.fc3
    ph.nac_params = sscha_yaml.nac_params
    ph.mesh_numbers = list(settings.bubble_mesh)
    ph.init_phph_interaction()
    ph.run_phonon_solver()
    frequencies, _, done = ph.phph_interaction.get_phonons()
    if not np.all(done) or not np.isfinite(frequencies).all():
        raise ValueError("Incomplete or non-finite auxiliary phonons.")
    if frequencies.min() < -1e-4:
        raise ValueError(
            f"Auxiliary phonon mesh has imaginary modes ({frequencies.min():.6g} THz); "
            "a perturbative on-shell bubble is not a stabilization procedure."
        )
    return ph


def _atomic_yaml(path, values):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(yaml.safe_dump(values, sort_keys=False), encoding="utf-8")
    temporary.replace(path)


def _checkpoint(path, fingerprint, frequency, delta):
    temporary = path.with_suffix(".tmp")
    with h5py.File(temporary, "w") as handle:
        handle.attrs["fingerprint"] = fingerprint
        handle["frequency"] = frequency
        handle["delta"] = delta
    temporary.replace(path)


def run_bubble_temperature(config, temperature_dir, temperature):
    """Postprocess one completed SSCHA temperature without altering its inputs."""
    temperature_dir = Path(temperature_dir).expanduser().resolve()
    settings = config.sections.scph
    paths = _input_paths(config, temperature_dir)
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(f"bubble input missing: {path}")
    state = yaml.safe_load(paths["sscha_summary"].read_text(encoding="utf-8")) or {}
    if state.get("status") != "complete" or not np.isclose(
        float(state.get("temperature", np.nan)), temperature, rtol=0, atol=1e-8
    ):
        raise ValueError("bubble requires a completed SSCHA summary at the requested temperature.")
    # Read cell size without allocating FC3 or any interaction workspaces.
    structure = PhonopyYaml().read(paths["sscha_structure"])
    if structure.primitive is None:
        raise ValueError("Saved phonopy_sscha.yaml has no primitive-cell metadata.")
    nb = 3 * len(structure.primitive)
    bound = int(np.prod(settings.bubble_mesh)) * nb**3 * 8
    if bound > MAX_INTERACTION_BYTES:
        raise ValueError(
            f"Dense bubble interaction bound is {bound / 1024**3:.2f} GiB (>2 GiB). "
            "This pilot implementation needs a smaller mesh/cell or a tiled backend. "
            "Do not reduce the physical model without checking its validity."
        )
    inputs = {
        "schema": SCHEMA_VERSION,
        "files": {key: file_identity(path) for key, path in paths.items()},
        "temperature_K": float(temperature),
        "mesh": list(settings.bubble_mesh),
        "epsilons_THz": list(settings.bubble_epsilons),
        "grid_points": (list(settings.bubble_grid_points)
                        if settings.bubble_grid_points is not None else None),
        "cutoff_frequency_THz": settings.cutoff_frequency,
        "phono3py_version": phono3py.__version__,
        "phonopy_version": phonopy.__version__,
        "backend": "C",
    }
    fingerprint = data_sha256(inputs)
    out = temperature_dir / "bubble" / fingerprint[:16]
    out.mkdir(parents=True, exist_ok=True)
    _atomic_yaml(out / "inputs.yaml", inputs)
    ph = _load_reference(paths, settings)
    ir, weights, _ = get_ir_grid_points(ph.grid)
    all_points = ph.grid.grg2bzg[ir]
    if settings.bubble_grid_points is None:
        points = all_points
        selected_weights = weights
        coverage = "full irreducible mesh"
    else:
        points = np.asarray(settings.bubble_grid_points, dtype=int)
        if np.any(points >= len(ph.grid.addresses)):
            raise ValueError("scph.bubble-grid-points contains an out-of-range BZ grid index.")
        # Explicit points may not be IR representatives; do not invent weights.
        selected_weights = None
        coverage = "selected BZ grid points; not a full-BZ integral"
    epsilons = np.asarray(settings.bubble_epsilons)
    frequencies, deltas = [], []
    for index, gp in enumerate(points):
        checkpoint = out / f"gp-{int(gp)}.hdf5"
        if checkpoint.is_file():
            with h5py.File(checkpoint, "r") as handle:
                if handle.attrs.get("fingerprint") != fingerprint:
                    raise ValueError(f"Bubble checkpoint input mismatch: {checkpoint}")
                frequency = handle["frequency"][:]
                delta = handle["delta"][:]
        else:
            # One q and one T per call avoid axis-order ambiguities in older
            # phono3py API documentation and allow exact q-point restarts.
            _, raw = ph.run_real_self_energy(
                grid_points=[int(gp)], temperatures=[float(temperature)],
                frequency_points_at_bands=True, epsilons=list(epsilons),
            )
            expected = (len(epsilons), 1, 1, nb)
            if np.shape(raw) != expected:
                raise RuntimeError(f"Unexpected phono3py bubble shape {np.shape(raw)} != {expected}")
            delta = np.asarray(raw)[:, 0, 0, :]
            frequency = ph.phph_interaction.get_phonons()[0][gp].copy()
            on_shell_frequencies(frequency, delta, settings.cutoff_frequency)
            _checkpoint(checkpoint, fingerprint, frequency, delta)
        if frequency.shape != (nb,) or delta.shape != (len(epsilons), nb):
            raise ValueError(f"Malformed bubble checkpoint: {checkpoint}")
        on_shell_frequencies(frequency, delta, settings.cutoff_frequency)
        frequencies.append(frequency)
        deltas.append(delta)
        print(f"  - Bubble {temperature:g} K: q {index + 1}/{len(points)} (gp={gp})", flush=True)
    frequency = np.asarray(frequencies)
    delta = np.stack(deltas, axis=1)  # epsilon, q, band
    corrections = [on_shell_frequencies(f, d, settings.cutoff_frequency)
                   for f, d in zip(frequencies, deltas)]
    arrays = {key: np.stack([c[key] for c in corrections], axis=1)
              for key in corrections[0]}
    qpoints = ph.grid.addresses[points] @ ph.grid.QDinv.T
    result = out / "bubble.hdf5"
    temporary = result.with_suffix(".tmp")
    with h5py.File(temporary, "w") as handle:
        handle.attrs["schema_version"] = SCHEMA_VERSION
        handle.attrs["fingerprint"] = fingerprint
        handle.attrs["frequency_unit"] = "THz (ordinary frequency, not angular)"
        handle.attrs["squared_frequency_unit"] = "THz^2"
        handle.attrs["delta_axes"] = "epsilon, grid_point, band"
        handle.attrs["method"] = bubble_scope()["label"]
        handle.attrs["transport_updated"] = False
        for key, value in {
            "temperature": [temperature], "epsilon": epsilons,
            "mesh": settings.bubble_mesh, "grid_point": points, "qpoint": qpoints,
            "frequency_auxiliary": frequency, "delta": delta, **arrays,
        }.items():
            handle.create_dataset(key, data=value)
        if selected_weights is not None:
            handle["weight"] = selected_weights
    temporary.replace(result)
    included = arrays["reference_mode_included"]
    ratio = arrays["relative_shift"]
    warnings = ["Mesh, epsilon, and SSCHA statistical convergence are not established."]
    invalid = int(np.count_nonzero(included & ~arrays["valid_mode"]))
    large = int(np.count_nonzero(ratio > 0.1))
    if invalid:
        warnings.append(f"{invalid} epsilon/q/band entries have nonpositive corrected squared frequencies.")
    if large:
        warnings.append(f"{large} epsilon/q/band entries have |Delta|/nu > 0.1; inspect the on-shell approximation.")
    summary = {
        "status": "complete", "temperature_K": float(temperature),
        "approximation": bubble_scope(), "fingerprint": fingerprint,
        "mesh": list(settings.bubble_mesh), "epsilons_THz": epsilons.tolist(),
        "coverage": coverage, "grid_point_count": len(points), "band_count": nb,
        "nac": ph.nac_params is not None,
        "gamma_nac_direction": "none; no directional LO limit at Gamma",
        "maximum_absolute_delta_THz": float(np.max(abs(delta))),
        "maximum_relative_shift": float(np.nanmax(ratio)) if np.any(included) else None,
        "epsilon_spread_max_THz": float(np.max(np.ptp(delta, axis=0))),
        "nonpositive_corrected_entries": invalid,
        "large_shift_entries": large,
        "result": result.name, "warnings": warnings,
    }
    _atomic_yaml(out / "bubble-summary.yaml", summary)
    print(f"  - Bubble output: {result}; auxiliary FC2 and transport unchanged")
    return str(out / "bubble-summary.yaml")


def run_existing_sscha_bubble(config):
    """Reuse existing fixed-volume or QHA-volume SSCHA results, with no ASE work."""
    root = Path(config.result_dir).expanduser().resolve()
    results = []
    for temperature in temperature_points(config.scph_temps):
        label = temperature_directory_name(temperature)
        case = copy.copy(config)
        if config.workflow_preset == "qha-sscha" or config.qha_sscha_enabled:
            case.result_dir = str(root / "qha-sscha" / label)
            # Coupled runs must use the FC3 and metadata at this temperature's
            # QHA volume, never a shared fixed-volume override.
            case.scph_transport_metadata = None
            case.scph_transport_fc3 = None
        directory = Path(case.result_dir) / case.scph_workdir / label
        results.append(run_bubble_temperature(case, directory, float(temperature)))
    return results
