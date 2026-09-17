"""Application-level command handlers.

This module separates user operations from terminal/logging concerns in
``cli`` and from backend implementation details in the workflow modules.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import os
import time


@dataclass
class CommandExecution:
    """Execution metadata returned to the CLI/provenance layer."""

    stage_timings: list[tuple[str, float]] = field(default_factory=list)

    def run_stage(self, label, operation):
        """Run one operation and record a consistently formatted duration."""
        start = time.time()
        try:
            return operation()
        finally:
            elapsed = time.time() - start
            self.stage_timings.append((label, elapsed))
            print(f"  - {label} elapsed time: {_format_duration(elapsed)}")


class WorkflowStageRunner:
    """Construct and run scientific stages for one shared workflow context."""

    def __init__(self, config, execution, *, core=None):
        self.config = config
        self.execution = execution
        self._core = core

    @property
    def core(self):
        """Lazily create the compatibility context shared by all core stages."""
        if self._core is None:
            from nepkappa.workflow import NEPPhononWorkflow

            self._core = NEPPhononWorkflow(self.config)
        return self._core

    def run(self, step):
        """Execute one named plan step through an explicit stage constructor."""
        handlers = {
            "relax": self.relax,
            "fc2": lambda: self.force_constants(include_fc3=False),
            "fc2fc3": lambda: self.force_constants(include_fc3=True),
            "fc4": self.fourth_order,
            "qha": self.qha,
            "scph": self.scph,
            "qha-sscha": self.qha_sscha,
            "kappa": self.kappa,
            "kappa4": self.kappa4,
            "plot": self.plot,
        }
        try:
            handler = handlers[step]
        except KeyError as exc:
            raise ValueError(f"Unsupported workflow command: {step}") from exc
        return handler()

    def relax(self):
        from nepkappa.stages.structure import StructureRelaxationStage

        return self.execution.run_stage(
            "Relax structure", StructureRelaxationStage(self.core).run
        )

    def force_constants(self, *, include_fc3):
        from nepkappa.force_constants import ForceConstantsWorkflow
        from nepkappa.stages.force_constants import (
            FiniteDisplacementStage,
            HiPhiveStage,
            ThirdOrderStage,
        )

        workflow = self.core
        strategy = ForceConstantsWorkflow(
            self.config,
            prepare_structure=workflow.load_force_constant_structure,
            timed_stage=self.execution.run_stage,
            submit_slurm=workflow._submit_force_slurm,
            run_finite_displacement=lambda include_fc3: FiniteDisplacementStage(
                workflow, include_fc3=include_fc3
            ).run(),
            run_hiphive=lambda include_fc3: HiPhiveStage(
                workflow, include_fc3=include_fc3
            ).run(),
            run_thirdorder=lambda: ThirdOrderStage(workflow).run(),
        )
        return strategy.run(include_fc3=include_fc3)

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

    def scph(self):
        from nepkappa.sscha import PhonopySSCHAWorkflow

        return self.execution.run_stage(
            "Phonopy SSCHA",
            PhonopySSCHAWorkflow(self.config, workflow=self.core).run,
        )

    def qha_sscha(self):
        from nepkappa.qha_sscha import QHASSCHAWorkflow

        return self.execution.run_stage(
            "QHA-volume SSCHA",
            QHASSCHAWorkflow(self.config, execution=self.execution).run,
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
        from nepkappa.plot import plot_results

        return self.execution.run_stage(
            "Plot results", lambda: plot_results(self.config)
        )


def execute_workflow_command(command, config) -> CommandExecution:
    """Execute one validated workflow command through its application handler."""
    if command == "run":
        return execute_workflow_plan(config)

    execution = CommandExecution()
    WorkflowStageRunner(config, execution).run(command)
    _print_timing_summary(execution.stage_timings)
    return execution


def execute_workflow_plan(config) -> CommandExecution:
    """Run a preset/custom stage plan behind the single ``nepkappa run`` entry."""
    execution = CommandExecution()
    steps = list(config.sections.workflow.steps)
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
        outcome = runner.run(step)
        if step in {"fc2", "fc2fc3"} and outcome is not None:
            _print_deferred(step)
            break
        if step == "kappa4" and isinstance(outcome, dict) and "submitted" in outcome:
            _print_deferred(step)
            break

    _print_timing_summary(execution.stage_timings)
    return execution


def _print_deferred(step):
    print(
        f"  - Plan paused after '{step}' because work was prepared or submitted. "
        "The dependent collection job will resume the remaining plan."
    )


def _print_timing_summary(stage_timings):
    if not stage_timings:
        return
    print("\n[Timing Summary]")
    for label, elapsed in stage_timings:
        print(f"  - {label:<24}: {_format_duration(elapsed)}")


def _format_duration(seconds):
    """Format a short stage duration without importing scientific backends."""
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    remaining = seconds % 60
    if hours:
        return f"{hours}h {minutes}m {remaining:.2f}s"
    if minutes:
        return f"{minutes}m {remaining:.2f}s"
    return f"{remaining:.2f}s"
