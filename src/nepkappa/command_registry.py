"""Declarative command catalogue shared by the CLI and help tools."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class CommandSpec:
    """User-facing metadata for one NEP-kappa command."""

    name: str
    help: str


CALCULATION_COMMANDS = ("relax", "fc2", "fc2fc3", "qha", "kappa", "tdbte")


COMMAND_SPECS = (
    CommandSpec("info", "Check and show parsed input settings."),
    CommandSpec("run", "Run the workflow defined by an input file."),
    CommandSpec("relax", "Relax the input structure."),
    CommandSpec("fc2", "Generate second-order force constants."),
    CommandSpec("fc2fc3", "Generate second- and third-order force constants."),
    CommandSpec("qha", "Run the configured QHA volume scan."),
    CommandSpec("kappa", "Compute conductivity using the input-selected route."),
    CommandSpec("tdbte", "Run time-dependent BTE from a dynamic input."),
    CommandSpec("plot", "Plot completed results."),
    CommandSpec("report", "Summarize completed results."),
)

PUBLIC_COMMANDS = tuple(spec.name for spec in COMMAND_SPECS)
EXECUTION_COMMANDS = frozenset({"run", *CALCULATION_COMMANDS, "plot"})
VALIDATION_TARGETS = ("run", *CALCULATION_COMMANDS, "plot")


def canonical_command(name: str) -> str:
    """Accept a public calculation or plotting command."""
    if name in EXECUTION_COMMANDS:
        return name
    raise ValueError(f"Unknown workflow command: {name}")


def resolve_calculation_command(name: str, config) -> str:
    """Select the internal operation for a public command from the input plan."""
    if name not in CALCULATION_COMMANDS:
        raise ValueError(f"Unknown calculation command: {name}")
    steps = tuple(config.workflow_steps)
    if name == "qha":
        if "qha" not in steps:
            raise ValueError("The input does not enable QHA; select a QHA input or use 'run'.")
        return "qha"
    if name == "kappa":
        transport_steps = tuple(
            step for step in steps
            if step in {"kappa", "kappa4", "qha-kappa", "qha-sscha"}
        )
        if not transport_steps:
            raise ValueError(
                "The input does not select a thermal-conductivity stage; use 'run' "
                "for its configured workflow."
            )
        if len(transport_steps) > 1:
            raise ValueError(
                "The input selects multiple transport stages; use 'run' to execute "
                "them in the configured order."
            )
        selected = transport_steps[0]
        if selected == "qha-sscha" and not (
            config.qha_sscha_three_phonon or config.qha_sscha_four_phonon
        ):
            raise ValueError(
                "The QHA-SSCHA input does not enable three- or four-phonon "
                "transport; use 'run' for its configured workflow."
            )
        return selected
    return name
