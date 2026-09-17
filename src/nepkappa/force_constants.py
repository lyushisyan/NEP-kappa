"""FC2/FC3 strategy selection independent of numerical backend details."""

from __future__ import annotations

import os
from pathlib import Path
from types import SimpleNamespace

from nepkappa.run_state import RunStateStore


class ForceConstantsWorkflow:
    """Route an FC2/FC3 request to finite, HiPhive, thirdorder, or Slurm."""

    def __init__(
        self,
        config,
        *,
        prepare_structure,
        timed_stage,
        submit_slurm,
        run_finite_displacement,
        run_hiphive,
        run_thirdorder,
    ):
        self.config = config
        self.settings = _force_constant_settings(config)
        self.prepare_structure = prepare_structure
        self.timed_stage = timed_stage
        self.submit_slurm = submit_slurm
        self.run_finite_displacement = run_finite_displacement
        self.run_hiphive = run_hiphive
        self.run_thirdorder = run_thirdorder

    def run(self, include_fc3=True):
        """Execute the selected strategy while preserving stage timing labels."""
        print("\n[Step 2] Generate Force Constants")
        self.prepare_structure()

        force_backend = str(
            self.settings.parallel.get("backend", "none")
        ).lower()
        collecting = os.environ.get("NEPKAPPA_FORCE_COLLECT") == "1"
        if force_backend == "slurm" and not collecting:
            return self.timed_stage(
                "Slurm force preparation",
                lambda: self.submit_slurm(include_fc3=include_fc3),
            )

        thirdorder_fc3 = include_fc3 and self.settings.fc3_backend == "thirdorder"
        native_include_fc3 = include_fc3 and not thirdorder_fc3
        if self.settings.use_hiphive:
            self.timed_stage(
                "HiPhive fitting",
                lambda: self.run_hiphive(include_fc3=native_include_fc3),
            )
        else:
            self.timed_stage(
                "Finite displacement",
                lambda: self.run_finite_displacement(
                    include_fc3=native_include_fc3
                ),
            )
        if thirdorder_fc3:
            self.timed_stage("thirdorder FC3", self.run_thirdorder)
        if collecting:
            manifest = (
                Path(getattr(self.config, "result_dir", "result")).resolve()
                / "force-slurm"
                / "submission.yaml"
            )
            if manifest.is_file():
                RunStateStore(manifest).mark_complete()
        return None


def _force_constant_settings(config):
    """Return the typed FC section or a small legacy compatibility view."""
    if hasattr(config, "sections"):
        return config.sections.force_constants
    return SimpleNamespace(
        parallel=dict(getattr(config, "force_parallel", {}) or {}),
        fc3_backend=getattr(config, "fc3_backend", "phono3py"),
        use_hiphive=bool(getattr(config, "use_hiphive", False)),
    )
