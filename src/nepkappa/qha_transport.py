"""Three-phonon conductivity at temperature-dependent QHA equilibrium volumes."""

from __future__ import annotations

import copy
from pathlib import Path

import numpy as np
from ase.io import read, write
import yaml

from nepkappa.qha_sscha import (
    load_qha_equilibrium_volumes,
    scale_atoms_to_primitive_volume,
    validate_temperature_range,
)
from nepkappa.sscha import temperature_directory_name


def _temperatures(values):
    if len(values) == 1:
        return np.asarray(values, dtype=float)
    minimum, maximum, step = (float(value) for value in values)
    if minimum < 0 or maximum < minimum or step <= 0:
        raise ValueError("kappa.temps must have 0 <= minimum <= maximum and step > 0.")
    count = int(np.floor((maximum - minimum) / step + 1.0e-10))
    return np.asarray([minimum + index * step for index in range(count + 1)])


class QHA3PhWorkflow:
    """Regenerate FC2/FC3 and run phono3py RTA at each expanded volume."""

    def __init__(self, config, *, execution):
        self.cfg = config
        self.execution = execution
        self.output_dir = Path(config.result_dir).expanduser().resolve()
        self.workdir = self.output_dir / "qha-kappa"

    def run(self):
        if not self.cfg.qha_enabled or not self.cfg.qha_volumes:
            raise ValueError("qha-kappa requires QHA and kappa.qha-volumes: true.")
        qha_temperatures, qha_volumes, multiplicity = load_qha_equilibrium_volumes(
            self.output_dir
        )
        requested = _temperatures(self.cfg.temps)
        validate_temperature_range(requested, qha_temperatures)
        volumes = np.interp(requested, qha_temperatures, qha_volumes)
        base = read(self.cfg.poscar)
        self.workdir.mkdir(parents=True, exist_ok=True)
        cases = []
        for temperature, primitive_volume in zip(requested, volumes):
            cases.append(
                self._run_temperature(
                    base, float(temperature), float(primitive_volume), multiplicity
                )
            )
        summary = {
            "status": "complete",
            "method": "phono3py RTA at QHA equilibrium volumes",
            "qha_volume_source": str(self.output_dir / "volume-temperature.dat"),
            "normalization_multiplicity": multiplicity,
            "force_constant_treatment": (
                "FC2 and FC3 regenerated at each temperature's QHA volume; "
                "no explicit anharmonic phonon renormalization"
            ),
            "cases": cases,
        }
        path = self.workdir / "qha-kappa-summary.yaml"
        path.write_text(yaml.safe_dump(summary, sort_keys=False), encoding="utf-8")
        return summary

    def _run_temperature(self, base, temperature, primitive_volume, multiplicity):
        from nepkappa.application import WorkflowStageRunner
        from nepkappa.qha import QHAWorkflow

        case_dir = self.workdir / temperature_directory_name(temperature)
        case_dir.mkdir(parents=True, exist_ok=True)
        atoms, scale = scale_atoms_to_primitive_volume(
            base, primitive_volume, multiplicity
        )
        case = copy.copy(self.cfg)
        case.result_dir = str(case_dir)
        case.poscar = str((case_dir / "POSCAR_qha_volume").resolve())
        case.do_relax = False
        case.workflow_preset = "custom"
        case.workflow_steps = ["fc2fc3", "kappa"]
        case.temps = [temperature]
        runner = WorkflowStageRunner(case, self.execution)
        if self.cfg.qha_relax_internal:
            atoms = QHAWorkflow(case, workflow=runner.core)._relax_internal(
                atoms, case_dir
            )
        write(case.poscar, atoms, format="vasp", direct=True, vasp5=True)
        print(
            f"\n[QHA 3ph] {temperature:g} K: "
            f"V_QHA={primitive_volume:.8f} A^3/primitive cell, "
            f"isotropic scale={scale:.8f}"
        )
        outcome = runner.run("fc2fc3")
        if outcome is not None:
            raise RuntimeError(
                "qha-kappa cannot continue from deferred force calculations; "
                "run inside one allocation with parallel.force-constant.backend: none."
            )
        kappa_result = runner.run("kappa")
        if isinstance(kappa_result, dict) and kappa_result.get("submitted"):
            raise RuntimeError("qha-kappa cannot continue from a deferred kappa job.")
        return {
            "temperature_K": temperature,
            "qha_primitive_volume_A3": primitive_volume,
            "unit_cell_volume_A3": float(atoms.get_volume()),
            "isotropic_scale": scale,
            "directory": str(case_dir),
        }
