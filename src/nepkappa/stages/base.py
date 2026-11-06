"""Minimal validation and execution contract for scientific stages."""

from __future__ import annotations

from abc import ABC, abstractmethod


class Stage(ABC):
    """Validate prerequisites, then execute one scientific stage."""

    name = "stage"

    def validate(self):
        """Validate stage-specific prerequisites before creating outputs."""

    @abstractmethod
    def execute(self):
        """Perform the scientific calculation and return its result."""

    def run(self):
        """Validate prerequisites and return the calculation result."""
        self.validate()
        return self.execute()


class HostStage(Stage):
    """A scientific stage using the shared workflow context."""

    def __init__(self, host):
        self.host = host
