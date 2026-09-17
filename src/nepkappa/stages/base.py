"""Common lifecycle contract for independently testable workflow stages."""

from __future__ import annotations

from abc import ABC, abstractmethod


class Stage(ABC):
    """Execute a scientific stage through a stable lifecycle.

    The no-op lifecycle hooks make local stages small today while providing
    the contract required by later scheduler and run-state adapters.
    """

    name = "stage"

    def validate(self):
        """Validate stage-specific prerequisites before creating outputs."""

    def prepare(self):
        """Prepare transient inputs needed by :meth:`execute`."""

    @abstractmethod
    def execute(self):
        """Perform the scientific calculation and return its result."""

    def collect(self, result):
        """Normalize generated artifacts and return the public stage result."""
        return result

    def status(self):
        """Return the local lifecycle state for diagnostics."""
        return "ready"

    def run(self):
        """Run the complete stage lifecycle."""
        self.validate()
        self.prepare()
        return self.collect(self.execute())


class HostStage(Stage):
    """Base class for a stage being incrementally extracted from a host workflow."""

    def __init__(self, host):
        self.host = host
