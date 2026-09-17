"""Declarative command catalogue shared by the CLI and help tools."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class CommandSpec:
    """User-facing metadata for one NEP-kappa command."""

    name: str
    help: str
    config_kind: str = "workflow"
    deprecated_alias_for: Optional[str] = None


COMMAND_SPECS = (
    CommandSpec("run", "Run the selected workflow preset or custom stage plan."),
    CommandSpec("relax", "Relax the input structure only."),
    CommandSpec("fc2", "Generate fc2.hdf5 and phono3py_disp.yaml only."),
    CommandSpec("fc2fc3", "Generate FC2 and FC3 using the configured backend."),
    CommandSpec("fc4", "Generate FORCE_CONSTANTS_4TH using Fourthorder."),
    CommandSpec("fc", "Deprecated alias for fc2fc3.", deprecated_alias_for="fc2fc3"),
    CommandSpec("qha", "Run an isotropic quasi-harmonic approximation workflow."),
    CommandSpec("scph", "Run self-consistent phonons with Phonopy SSCHA."),
    CommandSpec(
        "qha-sscha",
        "Run SSCHA at the temperature-dependent QHA equilibrium volumes.",
    ),
    CommandSpec("kappa", "Compute thermal conductivity from existing force constants."),
    CommandSpec("kappa4", "Compute 3ph+4ph thermal conductivity with FourPhonon."),
    CommandSpec("plot", "Plot phonon and thermal-transport results."),
    CommandSpec("status", "Show stored and live scheduler job states."),
    CommandSpec("report", "Generate a compact report from a result directory.", "report"),
    CommandSpec("init", "Interactively create a validated YAML input.", "initializer"),
    CommandSpec("compare", "Compare DFT and multiple potential-model results.", "compare"),
    CommandSpec("converge", "Prepare, run, and analyze a convergence study.", "convergence"),
    CommandSpec("info", "Print parsed workflow settings without running."),
    CommandSpec("validate", "Validate a configuration without running calculations."),
)

COMMANDS = {spec.name: spec for spec in COMMAND_SPECS}
EXECUTION_COMMANDS = frozenset(
    spec.name
    for spec in COMMAND_SPECS
    if spec.config_kind == "workflow"
    and spec.name not in {"info", "status", "validate"}
)
VALIDATION_TARGETS = tuple(
    spec.name
    for spec in COMMAND_SPECS
    if spec.name
    not in {
        "fc", "info", "init", "status", "report", "validate", "compare", "converge"
    }
)


def canonical_command(name: str) -> str:
    """Resolve a deprecated alias to its canonical workflow command."""
    spec = COMMANDS[name]
    return spec.deprecated_alias_for or spec.name
