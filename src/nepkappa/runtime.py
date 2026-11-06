"""Small runtime helpers shared by scientific workflow stages."""

from __future__ import annotations

import sys
import threading
import time

from ase import Atoms
from phonopy.structure.atoms import PhonopyAtoms

from nepkappa.output import print_detail, format_activity_elapsed

try:
    from tqdm.auto import tqdm
except ImportError:  # pragma: no cover - tqdm is a declared dependency
    tqdm = None


def progress_iter(iterable, enabled=True, **kwargs):
    """Return a tqdm iterator when available and enabled."""
    if enabled and tqdm is not None:
        return tqdm(iterable, **kwargs)
    return iterable


def run_activity_task(label, func, enabled=True, update_interval=1.0, log_interval=30.0):
    """Run a blocking task while showing elapsed activity feedback."""
    if not enabled:
        return func()

    result = {}

    def worker():
        try:
            result["value"] = func()
        except BaseException as exc:  # re-raised in the caller thread
            result["error"] = exc

    thread = threading.Thread(target=worker)
    thread.start()
    start = time.time()
    progress_stream = getattr(sys.stderr, "terminal_stream", sys.stderr)

    if tqdm is not None:
        with tqdm(
            total=None,
            desc=label,
            unit="s",
            bar_format="{desc}: {elapsed} elapsed",
            file=progress_stream,
            leave=False,
        ) as bar:
            while thread.is_alive():
                thread.join(update_interval)
                bar.update(update_interval)
    else:
        last_log = start
        while thread.is_alive():
            thread.join(update_interval)
            now = time.time()
            if now - last_log >= log_interval:
                print(
                    f"\r{label}: {format_activity_elapsed(now - start)} elapsed",
                    end="",
                    file=progress_stream,
                    flush=True,
                )
                last_log = now

    thread.join()
    if "error" in result:
        raise result["error"]
    print_detail(f"{label}: {format_activity_elapsed(time.time() - start)} elapsed")
    return result.get("value")


def ase_to_phonopy(ase_atoms):
    """Convert an ASE Atoms object to a PhonopyAtoms object."""
    return PhonopyAtoms(
        symbols=ase_atoms.get_chemical_symbols(),
        scaled_positions=ase_atoms.get_scaled_positions(),
        cell=ase_atoms.cell,
    )


def phonopy_to_ase(ph_atoms):
    """Convert a PhonopyAtoms object to an ASE Atoms object."""
    return Atoms(
        symbols=ph_atoms.symbols,
        positions=ph_atoms.positions,
        cell=ph_atoms.cell,
        pbc=True,
    )
