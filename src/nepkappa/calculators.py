"""Extensible ASE calculator backends for machine-learning potentials."""

from __future__ import annotations

import importlib
from importlib import metadata
import inspect
from pathlib import Path
import re

import numpy as np
from ase.calculators.calculator import Calculator

from nepkappa.provenance import canonical_data, file_identity


ENTRY_POINT_GROUP = "nepkappa.calculators"


class CalculatorLoadError(RuntimeError):
    """Raised when an external calculator cannot be resolved or constructed."""


def calculator_entry_points():
    """Return installed NEP-kappa calculator entry points across Python versions."""
    discovered = metadata.entry_points()
    if hasattr(discovered, "select"):
        return list(discovered.select(group=ENTRY_POINT_GROUP))
    return list(discovered.get(ENTRY_POINT_GROUP, []))


def import_factory(factory_path):
    """Import ``module:attribute`` (or ``module.attribute``) as a callable."""
    value = str(factory_path).strip()
    if ":" in value:
        module_name, attribute_path = value.split(":", 1)
    else:
        module_name, separator, attribute_path = value.rpartition(".")
        if not separator:
            raise CalculatorLoadError(
                "calculator.factory must use 'module:callable' or "
                "'module.callable' syntax."
            )
    if not module_name or not attribute_path:
        raise CalculatorLoadError(f"Invalid calculator factory: {factory_path!r}")
    try:
        value = importlib.import_module(module_name)
        for attribute in attribute_path.split("."):
            value = getattr(value, attribute)
    except (ImportError, AttributeError) as exc:
        raise CalculatorLoadError(
            f"Could not import calculator factory '{factory_path}': {exc}"
        ) from exc
    if not callable(value):
        raise CalculatorLoadError(
            f"Calculator factory '{factory_path}' is not callable."
        )
    return value


def resolve_factory(name, factory_path=None):
    """Resolve an explicit factory or a named third-party entry point."""
    if factory_path:
        factory = import_factory(factory_path)
        return factory, {
            "source": "factory",
            "factory": str(factory_path),
            "implementation": callable_identity(factory),
        }

    matches = [entry for entry in calculator_entry_points() if entry.name == name]
    if not matches:
        raise CalculatorLoadError(
            f"No calculator plugin named '{name}' is installed. Set "
            "calculator.factory to an ASE Calculator class/factory, or install "
            f"a package exposing the '{ENTRY_POINT_GROUP}' entry-point group."
        )
    if len(matches) > 1:
        providers = ", ".join(
            sorted(str(getattr(entry, "value", entry)) for entry in matches)
        )
        raise CalculatorLoadError(
            f"Multiple calculator plugins are registered as '{name}': {providers}"
        )
    entry = matches[0]
    try:
        factory = entry.load()
    except Exception as exc:
        raise CalculatorLoadError(
            f"Could not load calculator plugin '{name}': {exc}"
        ) from exc
    if not callable(factory):
        raise CalculatorLoadError(f"Calculator plugin '{name}' is not callable.")
    distribution = getattr(entry, "dist", None)
    return factory, {
        "source": "entry_point",
        "entry_point": name,
        "value": getattr(entry, "value", None),
        "distribution": getattr(distribution, "name", None),
        "distribution_version": getattr(distribution, "version", None),
        "implementation": callable_identity(factory),
    }


def callable_identity(factory):
    """Describe and hash the Python implementation behind a calculator factory."""
    source = inspect.getsourcefile(factory)
    return {
        "module": getattr(factory, "__module__", None),
        "qualname": getattr(factory, "__qualname__", getattr(factory, "__name__", None)),
        "source_file": file_identity(source) if source and Path(source).is_file() else None,
    }


def existing_files(value):
    """Find file-valued strings nested inside calculator keyword arguments."""
    paths = []
    if isinstance(value, dict):
        for item in value.values():
            paths.extend(existing_files(item))
    elif isinstance(value, (list, tuple)):
        for item in value:
            paths.extend(existing_files(item))
    elif isinstance(value, (str, Path)):
        path = Path(value).expanduser()
        if path.is_file():
            paths.append(path)
    return paths


class ASECalculatorBackend:
    """Construct and reuse any ASE-compatible calculator."""

    def __init__(self, name, factory=None, kwargs=None, model_files=None):
        self.name = str(name)
        self.factory_path = factory
        self.kwargs = dict(kwargs or {})
        self.model_files = [Path(path).expanduser() for path in (model_files or [])]
        self.factory, self.factory_metadata = resolve_factory(
            self.name, self.factory_path
        )
        self._calculator = None

    def calculator(self):
        """Create the configured calculator once and reuse it across structures."""
        if self._calculator is None:
            try:
                calculator = self.factory(**self.kwargs)
            except Exception as exc:
                raise CalculatorLoadError(
                    f"Could not construct calculator '{self.name}' with the configured "
                    f"kwargs: {exc}"
                ) from exc
            if not isinstance(calculator, Calculator) and not hasattr(
                calculator, "get_forces"
            ):
                raise CalculatorLoadError(
                    f"Factory for '{self.name}' returned {type(calculator).__name__}, "
                    "not an ASE-compatible calculator."
                )
            self._calculator = calculator
        return self._calculator

    def calculate_forces(self, atoms):
        """Attach the calculator and return validated eV/Angstrom forces."""
        atoms.calc = self.calculator()
        forces = np.asarray(atoms.get_forces(), dtype="double")
        expected_shape = (len(atoms), 3)
        if forces.shape != expected_shape:
            raise ValueError(
                f"Calculator '{self.name}' returned forces with shape {forces.shape}; "
                f"expected {expected_shape}."
            )
        if not np.all(np.isfinite(forces)):
            raise ValueError(f"Calculator '{self.name}' returned non-finite forces.")
        return forces

    def cache_inputs(self):
        """Return configuration and model identities for cache/provenance hashing."""
        paths = [*self.model_files, *existing_files(self.kwargs)]
        unique_paths = []
        seen = set()
        for path in paths:
            resolved = str(Path(path).resolve(strict=False))
            if resolved not in seen:
                seen.add(resolved)
                unique_paths.append(Path(path))
        return canonical_data(
            {
                "name": self.name,
                "factory": self.factory_metadata,
                "kwargs": self.kwargs,
                "model_files": [file_identity(path) for path in unique_paths],
            }
        )


class MACECalculatorBackend(ASECalculatorBackend):
    """Native MACE backend for local checkpoints and foundation-model factories."""

    def __init__(
        self,
        *,
        model=None,
        foundation=None,
        device="cpu",
        default_dtype="float64",
        kwargs=None,
        model_files=None,
    ):
        options = dict(kwargs or {})
        hashed_model_files = list(model_files or [])
        if foundation:
            foundation = str(foundation).strip().lower().replace("-", "_")
            if not re.fullmatch(r"[a-z][a-z0-9_]*", foundation):
                raise CalculatorLoadError(
                    "calculator.foundation must be a simple MACE factory suffix, "
                    "for example 'mp'."
                )
            factory = f"mace.calculators:mace_{foundation}"
            if model:
                options.setdefault("model", str(model))
                if Path(str(model)).expanduser().is_file():
                    hashed_model_files.append(model)
        else:
            if not model:
                raise CalculatorLoadError(
                    "A local MACE checkpoint is required when calculator.foundation "
                    "is not set."
                )
            model_path = Path(str(model)).expanduser()
            if not model_path.is_file():
                raise CalculatorLoadError(
                    f"MACE checkpoint not found: {model_path}. Use a trained local "
                    "checkpoint, or set calculator.foundation: mp for a MACE "
                    "foundation model."
                )
            factory = "mace.calculators:MACECalculator"
            options.setdefault("model_paths", [str(model_path)])
            hashed_model_files.append(model_path)

        options.setdefault("device", str(device))
        options.setdefault("default_dtype", str(default_dtype))
        try:
            super().__init__(
                "mace",
                factory=factory,
                kwargs=options,
                model_files=hashed_model_files,
            )
        except CalculatorLoadError as exc:
            raise CalculatorLoadError(
                f"Could not initialize MACE. Install the optional dependency with "
                f"'python -m pip install -e .[mace]'. Details: {exc}"
            ) from exc
        self.model = model
        self.foundation = foundation


def calculator_backend_from_config(config):
    """Construct the native MACE or generic ASE/plugin backend from configuration."""
    getter = config.get if isinstance(config, dict) else lambda key, default=None: getattr(
        config, key, default
    )
    name = str(getter("calculator", "ase")).lower()
    if name == "mace":
        return MACECalculatorBackend(
            model=getter("mace_model"),
            foundation=getter("mace_foundation"),
            device=getter("mace_device", "cpu") or "cpu",
            default_dtype=getter("mace_dtype", "float64") or "float64",
            kwargs=getter("calculator_kwargs", {}),
            model_files=getter("calculator_model_files", []),
        )
    return ASECalculatorBackend(
        name,
        factory=getter("calculator_factory"),
        kwargs=getter("calculator_kwargs", {}),
        model_files=getter("calculator_model_files", []),
    )
