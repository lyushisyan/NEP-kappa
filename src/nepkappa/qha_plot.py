"""QHA overview from completed fits and FC2 recalculated at equilibrium volumes."""

from pathlib import Path
import warnings

import numpy as np
import yaml

from nepkappa.plot import (
    plt, build_dos_data, journal_style, load_phono3py_yaml, make_band_paths,
    make_phonopy,
)
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from phonopy.qha.eos import get_eos


def read_free_energy_curves(path):
    """Read Phonopy's temperature blocks, including the fitted EOS parameters."""
    curves = []
    current = None
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith("# Temperature:"):
            current = {"temperature": float(line.split(":", 1)[1]), "points": []}
            curves.append(current)
        elif line.startswith("# Parameters:") and current is not None:
            current["parameters"] = np.asarray(line.split(":", 1)[1].split(), float)
        elif line.strip() and not line.startswith("#"):
            if current is None:
                raise ValueError("Free-energy data must start with a temperature header.")
            current["points"].append([float(value) for value in line.split()])
    for curve in curves:
        points = np.asarray(curve["points"], float)
        parameters = curve.get("parameters", np.array([]))
        if (points.ndim != 2 or points.shape[1] != 2 or len(points) < 5
                or parameters.shape != (4,) or not np.isfinite(points).all()
                or not np.isfinite(parameters).all() or not np.isfinite(curve["temperature"])
                or np.any(np.diff(points[:, 0]) <= 0)):
            raise ValueError("Invalid QHA free-energy curve or EOS parameters.")
        curve["points"] = points
    if not curves:
        raise ValueError("No QHA free-energy curves found.")
    return curves


def read_temperature_data(path):
    values = np.loadtxt(path, ndmin=2)
    if (values.shape[1] != 2 or not np.isfinite(values).all()
            or np.any(np.diff(values[:, 0]) <= 0)):
        raise ValueError(f"Expected finite, increasing temperature/value pairs in {path}")
    return values


def load_qha_plot_data(config):
    """Read existing data only; never instantiate a force calculator."""
    root = Path(config.result_dir).resolve()
    summary = yaml.safe_load((root / "qha-summary.yaml").read_text())
    manifest_path = root / "qha-phonons.yaml"
    if not manifest_path.is_file():
        raise FileNotFoundError(
            "QHA phonon overview requires FC2 at equilibrium volumes. Set "
            "force-constant.qha.equilibrium-phonons: true and run "
            "nepkappa fc2 on the QHA input before plotting."
        )
    manifest = yaml.safe_load(manifest_path.read_text())
    if summary.get("status") != "complete" or manifest.get("status") != "complete":
        raise ValueError("QHA fitting and equilibrium phonons must both be complete.")
    equilibrium = read_temperature_data(root / "volume-temperature.dat")
    expansion = read_temperature_data(root / "thermal_expansion.dat")
    cases = manifest["cases"]
    case_values = np.array([[case["temperature_K"], case["primitive_volume_A3"]] for case in cases])
    if case_values.shape != equilibrium.shape or not np.allclose(case_values, equilibrium, rtol=0, atol=1e-7):
        raise ValueError("Equilibrium phonons do not match this QHA volume-temperature table.")
    # Keep a dense temperature grid readable while retaining its endpoints.
    indices = np.unique(np.linspace(0, len(cases) - 1, min(10, len(cases))).round().astype(int))
    selected = [cases[i] for i in indices]
    temperatures = np.array([case["temperature_K"] for case in selected])
    all_curves = read_free_energy_curves(root / "helmholtz-volume.dat")
    curves = []
    for temperature in temperatures:
        matches = [curve for curve in all_curves if np.isclose(curve["temperature"], temperature, rtol=0, atol=1e-6)]
        if len(matches) != 1:
            raise ValueError(f"Missing or ambiguous QHA free energy at {temperature:g} K")
        curves.append(matches[0])
    mesh = manifest["mesh"]
    spectra = []
    reference = None
    for point in summary["points"]:
        if np.isclose(point["volume_ratio"], 1.0, rtol=0, atol=1e-8):
            # Resolve under the current result directory, including downloaded results.
            directory = root / "qha-points" / Path(point["directory"]).name
            reference = _load_phonon(directory)
            break
    anchor = reference if reference is not None else _load_phonon(root / selected[0]["directory"])
    paths, labels = make_band_paths(anchor.primitive, config)
    anchor.run_band_structure(paths)
    distances = anchor.get_band_structure_dict()["distances"]
    for case in selected:
        phonon = _load_phonon(root / case["directory"])
        if (not np.array_equal(phonon.primitive.numbers, anchor.primitive.numbers)
                or not np.allclose(phonon.primitive.cell / phonon.primitive.volume ** (1/3),
                                   anchor.primitive.cell / anchor.primitive.volume ** (1/3), atol=1e-6)):
            raise ValueError("QHA phonon primitive bases differ; cannot overlay a common band path.")
        if not np.isclose(phonon.primitive.volume, case["primitive_volume_A3"], rtol=1e-6):
            raise ValueError("QHA phonon cell volume disagrees with the equilibrium manifest.")
        spectra.append(_spectrum(phonon, paths, mesh))
    reference_data = _spectrum(reference, paths, mesh) if reference is not None else None
    return dict(equilibrium=equilibrium, expansion=expansion, temperatures=temperatures,
                curves=curves, eos=summary["eos"], pressure=summary.get("pressure_GPa", 0),
                spectra=spectra, reference=reference_data, distances=distances, labels=labels)


def _load_phonon(directory):
    return make_phonopy(load_phono3py_yaml(directory / "phonopy_params.yaml"),
                        directory / "force_constants.hdf5")


def _spectrum(phonon, paths, mesh):
    phonon.run_band_structure(paths)
    frequencies = phonon.get_band_structure_dict()["frequencies"]
    if min(np.min(values) for values in frequencies) < -0.1:
        warnings.warn("QHA phonon plot contains imaginary frequencies; check structural stability.", UserWarning)
    return {"frequencies": frequencies, "dos": build_dos_data(phonon, mesh)}


def _draw_free_energy(ax, data, colors):
    eos = get_eos(data["eos"])
    minimum_v, minimum_f = [], []
    for curve, color in zip(data["curves"], colors):
        points, parameters = curve["points"], curve["parameters"]
        volumes = np.linspace(points[0, 0], points[-1, 0], 200)
        ax.plot(volumes, eos(volumes, parameters), color=color, lw=1)
        ax.plot(points[:, 0], points[:, 1], ".", color=color, ms=3)
        minimum_v.append(parameters[3])
        minimum_f.append(parameters[0])
    ax.plot(minimum_v, minimum_f, "o--", color="0.3", ms=3, mfc="white", lw=.8, label=r"$V_{\rm eq}(T)$")
    ax.legend(frameon=False, fontsize=8)
    ax.set_xlabel(r"$V$ ($\mathrm{\AA}^3$/primitive cell)")
    label = "$F$" if data["pressure"] == 0 else "$F+PV$"
    ax.set_ylabel(label + " (eV/primitive cell)")


def _draw_volume(ax, data):
    ax.plot(*data["equilibrium"].T, color="#009c9c", lw=1.5)
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel(r"$V$ ($\mathrm{\AA}^3$/primitive cell)")


def _draw_expansion(ax, data):
    values = data["expansion"]
    ax.plot(values[:, 0], values[:, 1] * 1e6, color="#009c9c", lw=1.5)
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel(r"$\alpha_V$ ($10^{-6}$ K$^{-1}$)")


def _draw_phonons(fig, slot, data, colors):
    grid = slot.subgridspec(1, 2, width_ratios=[4, 1], wspace=.05)
    band = fig.add_subplot(grid[0])
    dos = fig.add_subplot(grid[1], sharey=band)
    for spectrum, color in zip(data["spectra"], colors):
        for distances, frequencies in zip(data["distances"], spectrum["frequencies"]):
            band.plot(distances, frequencies, color=color, lw=.9)
        values = spectrum["dos"]
        dos.plot(values["total_dos"], values["frequency_points"], color=color, lw=.9)
    if data["reference"] is not None:
        for distances, frequencies in zip(data["distances"], data["reference"]["frequencies"]):
            band.plot(distances, frequencies, color="0.5", ls="--", lw=.8)
        values = data["reference"]["dos"]
        dos.plot(values["total_dos"], values["frequency_points"], color="0.5", ls="--", lw=.8)
        band.plot([], [], "--", color="0.5", label="Static-volume reference")
        band.legend(frameon=False, fontsize=7, loc="upper right")
    ticks, labels = [], []
    for i, (distances, segment) in enumerate(zip(data["distances"], data["labels"])):
        if i == 0:
            ticks.append(distances[0]); labels.append(segment[0])
        elif labels[-1] != segment[0]:
            labels[-1] += "|" + segment[0]
        ticks.append(distances[-1]); labels.append(segment[1])
        band.axvline(distances[-1], color="0.9", lw=.7, zorder=0)
    band.set_xticks(ticks, labels, fontsize=8)
    band.set_xlim(ticks[0], ticks[-1])
    band.set_ylabel("Frequency (THz)")
    dos.set_xlabel("DOS")
    dos.set_xlim(left=0)
    dos.tick_params(axis="y", labelleft=False)
    dos.tick_params(axis="x", labelbottom=False)
    return band, dos


def plot_qha_results(config):
    """Write QHA (a)--(d) panels with the standard separate/combined/both layout."""
    data = load_qha_plot_data(config)
    temperatures = data["temperatures"]
    norm = Normalize(float(temperatures.min()), float(temperatures.max()))
    mapper = ScalarMappable(norm=norm, cmap="viridis")
    colors = mapper.to_rgba(temperatures)
    output = Path(config.result_dir).resolve() / "plots"
    output.mkdir(exist_ok=True)
    names = ["qha_free_energy", "qha_volume", "qha_thermal_expansion", "qha_phonons"]
    layout = config.plot_layout
    saved = []

    def save(fig, name):
        for suffix in ("png", "pdf"):
            path = output / f"{name}.{suffix}"
            fig.savefig(path, dpi=config.plot_dpi, bbox_inches="tight")
            saved.append(path)
        plt.close(fig)

    def draw(fig, grid):
        axes = [fig.add_subplot(grid[i]) for i in range(3)]
        _draw_free_energy(axes[0], data, colors)
        _draw_volume(axes[1], data)
        _draw_expansion(axes[2], data)
        band, dos = _draw_phonons(fig, grid[3], data, colors)
        axes.append(band)
        for i, ax in enumerate(axes):
            ax.text(-.14, 1.06, f"({chr(97+i)})", transform=ax.transAxes, fontsize=12)
            ax.grid(color="0.9", linewidth=.6)
        color_axis = fig.add_axes([.32, .94, .48, .018])
        fig.colorbar(mapper, cax=color_axis, orientation="horizontal", label="Temperature (K)")

    with journal_style(), plt.rc_context({"font.size": 10, "axes.labelsize": 11,
                                         "xtick.labelsize": 9, "ytick.labelsize": 9}):
        if layout in ("combined", "both"):
            fig = plt.figure(figsize=(10, 7.5))
            draw(fig, fig.add_gridspec(2, 2, left=.09, right=.98, bottom=.09,
                                      top=.84, wspace=.35, hspace=.42))
            save(fig, "qha_combined")
        if layout in ("separate", "both"):
            for i, name in enumerate(names):
                fig = plt.figure(figsize=(6, 4))
                fig.subplots_adjust(left=.16, right=.96, bottom=.16, top=.78 if i in (0, 3) else .94)
                if i == 3:
                    axes = _draw_phonons(fig, fig.add_gridspec(1, 1)[0], data, colors)
                else:
                    ax = fig.add_subplot()
                    if i == 0:
                        _draw_free_energy(ax, data, colors)
                    elif i == 1:
                        _draw_volume(ax, data)
                    else:
                        _draw_expansion(ax, data)
                    ax.grid(color="0.9", linewidth=.6)
                    axes = [ax]
                if i in (0, 3):
                    color_axis = fig.add_axes([.3, .91, .5, .025])
                    fig.colorbar(mapper, cax=color_axis, orientation="horizontal", label="Temperature (K)")
                save(fig, name)
    print(f"  - QHA overview: {len(temperatures)} temperatures; FC2 recalculated at equilibrium volumes")
    print(f"  - Generated {len(saved)} QHA plot files: {output}")
    return saved
