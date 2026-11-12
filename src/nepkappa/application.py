"""Application-level command handlers.

This module separates user operations from terminal/logging concerns in
``cli`` and from backend implementation details in the workflow modules.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import os
import time

from nepkappa.output import print_detail, format_duration
from nepkappa.command_registry import calculation_steps


@dataclass
class CommandExecution:
    """Command outcome: complete locally, or prepared/submitted for later work.

    Exceptions propagate to the CLI as failures. Successful scheduler operations
    retain their manifest status instead of claiming scientific completion.
    """

    stage_timings: list[tuple[str, float]] = field(default_factory=list)
    status: str = "complete"

    def run_stage(self, label, operation):
        """Run one operation and record a consistently formatted duration."""
        start = time.time()
        try:
            result = operation()
            if isinstance(result, dict) and result.get("status") in {"prepared", "submitted"}:
                self.status = result["status"]
            return result
        finally:
            elapsed = time.time() - start
            self.stage_timings.append((label, elapsed))
            print_detail(f"  - {label} elapsed time: {format_duration(elapsed)}")


class WorkflowStageRunner:
    """Construct and run scientific stages for one shared workflow context."""

    def __init__(self, config, execution, *, core=None):
        self.config = config
        self.execution = execution
        self._core = core

    @property
    def core(self):
        """Lazily create the structure/calculator context shared by core stages."""
        if self._core is None:
            from nepkappa.workflow import NEPPhononWorkflow

            self._core = NEPPhononWorkflow(self.config)
        return self._core

    def run(self, step):
        """Execute one named plan step through an explicit stage constructor."""
        handlers = {
            "tdbte": self.tdbte,
            "relax": self.relax,
            "fc2": lambda: self.force_constants(include_fc3=False),
            "fc2fc3": lambda: self.force_constants(include_fc3=True),
            "fc4": self.fourth_order,
            "qha": self.qha,
            "qha-kappa": self.qha_kappa,
            "scph": self.scph,
            "bubble": self.bubble,
            "qha-sscha": self.qha_sscha,
            "kappa": self.kappa,
            "kappa4": self.kappa4,
            "plot": self.plot,
            "report": self.report,
        }
        try:
            handler = handlers[step]
        except KeyError as exc:
            raise ValueError(f"Unsupported workflow command: {step}") from exc
        return handler()

    def tdbte(self):
        from nepkappa.tdbte import run_tdbte

        return self.execution.run_stage("Time-dependent BTE", lambda: run_tdbte(self.config))

    def relax(self):
        from nepkappa.stages.structure import StructureRelaxationStage

        return self.execution.run_stage(
            "Relax structure", StructureRelaxationStage(self.core).run
        )

    def force_constants(self, *, include_fc3):
        """Prepare the structure, then fit locally or defer forces to Slurm."""
        from pathlib import Path

        from nepkappa.run_state import RunStateStore
        from nepkappa.stages.force_constants import (
            FiniteDisplacementStage, HiPhiveStage, ThirdOrderStage,
        )

        config, core = self.config, self.core
        print("\n[Step 2] Generate Force Constants")
        core.load_force_constant_structure()
        collecting = os.environ.get("NEPKAPPA_FORCE_COLLECT") == "1"
        parallel = getattr(config, "force_parallel", {}) or {}
        if str(parallel.get("backend", "none")).lower() == "slurm" and not collecting:
            return self.execution.run_stage(
                "Slurm force preparation",
                lambda: core._submit_force_slurm(include_fc3=include_fc3),
            )

        thirdorder = (
            include_fc3 and getattr(config, "fc3_backend", "phono3py") == "thirdorder"
        )
        native_fc3 = include_fc3 and not thirdorder
        if getattr(config, "use_hiphive", False):
            stage = HiPhiveStage(core, include_fc3=native_fc3)
            label = "HiPhive fitting"
        else:
            stage = FiniteDisplacementStage(core, include_fc3=native_fc3)
            label = "Finite displacement"
        self.execution.run_stage(label, stage.run)
        if thirdorder:
            self.execution.run_stage("thirdorder FC3", ThirdOrderStage(core).run)
        if collecting:
            manifest = (
                Path(getattr(config, "result_dir", "result")).resolve()
                / "force-slurm" / "submission.yaml"
            )
            if manifest.is_file():
                RunStateStore(manifest).mark_complete()
        return None

    def fourth_order(self):
        from nepkappa.stages.force_constants import FourthOrderStage

        print("\n[Step 2] Generate Fourth-Order Force Constants")
        self.core.load_force_constant_structure()
        return self.execution.run_stage(
            "FourPhonon FC4", FourthOrderStage(self.core).run
        )

    def qha(self):
        from nepkappa.qha import QHAWorkflow

        return self.execution.run_stage(
            "QHA", QHAWorkflow(self.config, workflow=self.core).run
        )

    def qha_kappa(self):
        from nepkappa.qha_transport import QHA3PhWorkflow

        return self.execution.run_stage(
            "QHA-volume three-phonon kappa",
            QHA3PhWorkflow(
                self.config, execution=self.execution,
                runner_factory=lambda case: WorkflowStageRunner(case, self.execution),
            ).run,
        )

    def scph(self):
        from nepkappa.sscha import PhonopySSCHAWorkflow

        return self.execution.run_stage(
            "Phonopy SSCHA",
            PhonopySSCHAWorkflow(
                self.config, workflow=self.core, execution=self.execution
            ).run,
        )

    def bubble(self):
        from nepkappa.bubble import run_existing_sscha_bubble

        return self.execution.run_stage(
            "SSCHA on-shell bubble", lambda: run_existing_sscha_bubble(self.config)
        )

    def qha_sscha(self):
        from nepkappa.qha_sscha import QHASSCHAWorkflow

        return self.execution.run_stage(
            "QHA-volume SSCHA",
            QHASSCHAWorkflow(
                self.config, execution=self.execution,
                runner_factory=lambda case: WorkflowStageRunner(case, self.execution),
            ).run,
        )

    def kappa(self):
        from nepkappa.transport import Phono3pyTransportWorkflow

        workflow = self.core
        transport = Phono3pyTransportWorkflow(
            self.config,
            workflow.output_dir,
            run_command=workflow._run_command,
            command_resolver=workflow._phono3py_command,
            needs_fc_flags=workflow._phono3py_needs_fc_flags,
        )
        return self.execution.run_stage("Phono3py kappa", transport.run)

    def kappa4(self):
        from nepkappa.fourphonon import FourPhononWorkflow

        return self.execution.run_stage(
            "FourPhonon kappa", FourPhononWorkflow(self.config).run
        )

    def plot(self):
        import copy
        from pathlib import Path

        trajectory = Path(self.config.result_dir) / "tdbte" / "branch-energy.csv"
        spatial = Path(self.config.result_dir) / "tdbte" / "spatial.hdf5"
        if "tdbte" in self.config.workflow_steps or (
            self.config.workflow_steps == ["plot"] and (trajectory.is_file() or spatial.is_file())
        ):
            from nepkappa.tdbte_plot import plot_tdbte_results

            return self.execution.run_stage(
                "Plot TD-BTE results",
                lambda: plot_tdbte_results(self.config),
            )

        from nepkappa.plot import plot_results

        config = self.config
        if config.scph_enabled and not (config.scph_run_transport or config.fp_enabled):
            from nepkappa.sscha_plot import plot_sscha_results

            return self.execution.run_stage(
                "Plot SSCHA results", lambda: plot_sscha_results(config)
            )
        if not config.scph_enabled and not config.qha_volumes and (
            config.qha_enabled or (Path(config.result_dir) / "qha-summary.yaml").is_file()
        ):
            from nepkappa.qha_plot import plot_qha_results

            return self.execution.run_stage(
                "Plot QHA results", lambda: plot_qha_results(config)
            )
        if config.scph_enabled or (config.qha_enabled and config.qha_volumes):
            from nepkappa.temperatures import temperature_directory_name

            config = copy.copy(config)
            temperature_dir = temperature_directory_name(config.plot_temperature)
            result_dir = Path(config.result_dir)
            if config.scph_enabled:
                if config.qha_enabled:
                    result_dir /= Path("qha-sscha") / temperature_dir
                config._plot_harmonic_reference = str(
                    Path(getattr(config, "scph_initial_fc2", None) or result_dir / "fc2.hdf5")
                )
                result_dir /= Path(config.scph_workdir) / temperature_dir
            else:
                result_dir /= Path("qha-kappa") / temperature_dir
            config.result_dir = str(result_dir)
        return self.execution.run_stage(
            "Plot results", lambda: plot_results(config)
        )

    def report(self):
        from nepkappa.report import generate_report

        return self.execution.run_stage(
            "Report results",
            lambda: generate_report(self.config.result_dir),
        )


def execute_workflow_command(command, config) -> CommandExecution:
    """Execute one validated workflow command through its application handler."""
    if command == "run":
        return execute_workflow_plan(config)
    if command == "fc2fc3fc4":
        return execute_workflow_plan(config, steps=calculation_steps(command))

    execution = CommandExecution()
    WorkflowStageRunner(config, execution).run(command)
    _print_timing_summary(execution.stage_timings)
    return execution


def execute_workflow_plan(config, *, steps=None) -> CommandExecution:
    """Run an ordered stage sequence, pausing when scheduler work is deferred."""
    execution = CommandExecution()
    steps = list(config.workflow_steps if steps is None else steps)
    resume_after = os.environ.get("NEPKAPPA_RESUME_AFTER")
    if resume_after:
        if resume_after not in steps:
            raise ValueError(f"Cannot resume after absent stage '{resume_after}'.")
        steps = steps[steps.index(resume_after) + 1 :]
        print(f"[Plan] Resuming after completed stage '{resume_after}'")
    if os.environ.get("NEPKAPPA_FORCE_COLLECT") == "1":
        force_indices = [
            index for index, step in enumerate(steps) if step in {"fc2", "fc2fc3"}
        ]
        if force_indices:
            steps = steps[force_indices[0] :]
            print("[Plan] Resuming from completed Slurm force calculations")
    print(f"[Plan] {' -> '.join(steps)}")

    runner = WorkflowStageRunner(config, execution)

    for step in steps:
        runner.run(step)
        if execution.status != "complete":
            _print_deferred(step)
            break

    _print_timing_summary(execution.stage_timings)
    return execution


def _print_deferred(step):
    print(
        f"  - Plan paused after '{step}' because work was prepared or submitted. "
        "Scientific calculations are not yet complete."
    )


def _print_timing_summary(stage_timings):
    if not stage_timings:
        return
    print_detail("\n[Timing Summary]")
    for label, elapsed in stage_timings:
        print_detail(f"  - {label:<24}: {format_duration(elapsed)}")
