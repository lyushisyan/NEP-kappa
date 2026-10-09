"""Declarative command catalogue shared by the CLI and help tools."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class CommandSpec:
    """User-facing metadata for one NEP-kappa command."""

    name: str
    help: str


STAGE_COMMANDS = (
    "relax", "fc2", "fc2fc3", "fc4", "qha", "qha-kappa", "scph",
    "bubble", "qha-sscha", "kappa", "kappa4", "tdbte",
)


COMMAND_SPECS = (
    CommandSpec("run", "Run the workflow defined by an input file."),
    CommandSpec("stage", "Run one calculation stage (see 'stage --help')."),
    CommandSpec("plot", "Plot completed results."),
    CommandSpec("status", "Show saved and scheduler job states."),
    CommandSpec("report", "Summarize completed results."),
    CommandSpec("info", "Check and show parsed input settings."),
)

PUBLIC_COMMANDS = tuple(spec.name for spec in COMMAND_SPECS)
EXECUTION_COMMANDS = frozenset({"run", "plot"})
VALIDATION_TARGETS = ("run", *STAGE_COMMANDS, "plot")


def canonical_command(name: str) -> str:
    """Accept a public workflow command or a named calculation stage."""
    if name in STAGE_COMMANDS or name in {"run", "plot"}:
        return name
    raise ValueError(f"Unknown workflow command: {name}")
