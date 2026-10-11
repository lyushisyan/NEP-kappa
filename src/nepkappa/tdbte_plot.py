"""Plots for completed time-dependent BTE trajectories."""

from __future__ import annotations

import os
from pathlib import Path
import tempfile

import numpy as np


def _read_branch_energy(path: Path):
    if not path.is_file():
        raise FileNotFoundError(f"TD-BTE branch-energy data not found: {path}")
    with path.open(encoding="utf-8") as handle:
        header = handle.readline().strip().split(",")
    if len(header) < 2 or header[0] != "time_ps":
        raise ValueError("TD-BTE branch-energy.csv has an invalid header.")
    data = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    if data.shape[0] < 2 or data.shape[1] != len(header):
        raise ValueError("TD-BTE branch-energy.csv has an invalid data shape.")
    if not np.isfinite(data).all() or not np.all(np.diff(data[:, 0]) > 0):
        raise ValueError("TD-BTE branch-energy.csv needs finite data and increasing time.")
    labels = []
    for name in header[1:]:
        if not name.startswith("branch_") or not name.endswith("_excess_eV_per_cell"):
            raise ValueError(f"Unrecognized TD-BTE branch-energy column: {name}")
        labels.append(name[len("branch_"):-len("_excess_eV_per_cell")])
    return data[:, 0], data[:, 1:], labels


def plot_tdbte_results(config):
    """Draw branch-resolved energy transfer and the total-energy check."""
    result_dir = Path(config.result_dir).resolve()
    times, energies, labels = _read_branch_energy(
        result_dir / "tdbte" / "branch-energy.csv"
    )
    layout = config.plot_layout
    dpi = int(config.plot_dpi)
    if dpi <= 0:
        raise ValueError("plot.dpi must be a positive integer.")

    cache_dir = Path(tempfile.gettempdir()) / "nepkappa-matplotlib"
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plot_dir = result_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    paths = []

    def draw_branches(ax):
        for index, label in enumerate(labels):
            ax.plot(times, energies[:, index], label=f"Branch {label}")
        ax.set(xlabel="Time (ps)", ylabel="Excess energy (eV / primitive cell)",
               title="Energy redistribution")
        ax.legend(fontsize="small", ncol=2)
        ax.grid(alpha=0.3)

    def draw_total(ax):
        total = energies.sum(axis=1)
        ax.plot(times, total, color="tab:purple")
        ax.set(xlabel="Time (ps)", ylabel="Excess energy (eV / primitive cell)",
               title="Total excess energy")
        ax.grid(alpha=0.3)

    if layout in {"separate", "both"}:
        for name, draw in (("tdbte_branch_energy", draw_branches),
                           ("tdbte_total_energy", draw_total)):
            figure, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
            try:
                draw(ax)
                path = plot_dir / f"{name}.png"
                figure.savefig(path, dpi=dpi)
                paths.append(path)
            finally:
                plt.close(figure)
    if layout in {"combined", "both"}:
        figure, axes = plt.subplots(1, 2, figsize=(12, 4.5),
                                   constrained_layout=True)
        try:
            draw_branches(axes[0])
            draw_total(axes[1])
            path = plot_dir / "tdbte_combined.png"
            figure.savefig(path, dpi=dpi)
            paths.append(path)
        finally:
            plt.close(figure)
    print("  - TD-BTE plots:")
    for path in paths:
        print(f"    {path}")
    return paths
