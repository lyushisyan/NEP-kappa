"""Deterministic parameter-sweep preparation and convergence analysis."""

from __future__ import annotations

from dataclasses import dataclass
import csv
from pathlib import Path
import re
from typing import Any, Callable

import h5py
import matplotlib.pyplot as plt
import numpy as np
import yaml


COMPONENTS = {"x": 0, "y": 1, "z": 2}


@dataclass(frozen=True)
class ConvergenceStudyConfig:
    """Validated convergence-study settings."""

    config_path: Path
    base: Path
    parameter: str
    values: tuple[Any, ...]
    directory: Path
    execute: bool
    temperature: float
    component: str
    tolerance: float


@dataclass(frozen=True)
class ConvergenceCase:
    """One generated sweep point."""

    index: int
    value: Any
    label: str
    config_path: Path
    result_dir: Path


def parse_convergence_args(config_path) -> ConvergenceStudyConfig:
    """Read and validate a convergence-study YAML file."""
    path = Path(config_path).resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Convergence input file not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError("Convergence YAML must be a mapping.")
    allowed = {"base", "parameter", "values", "study"}
    unknown = set(data) - allowed
    if unknown:
        raise ValueError(f"Unknown convergence key: {sorted(unknown)[0]}")
    for key in ("base", "parameter", "values"):
        if key not in data:
            raise ValueError(f"Convergence YAML requires '{key}'.")
    if not isinstance(data["values"], list) or len(data["values"]) < 2:
        raise ValueError("Convergence 'values' must contain at least two points.")
    serialized_values = [yaml.safe_dump(value, sort_keys=True) for value in data["values"]]
    if len(serialized_values) != len(set(serialized_values)):
        raise ValueError("Convergence 'values' must not contain duplicates.")
    study = data.get("study", {}) or {}
    if not isinstance(study, dict):
        raise ValueError("Convergence 'study' must be a mapping.")
    allowed_study = {"directory", "execute", "temperature", "component", "tolerance"}
    unknown_study = set(study) - allowed_study
    if unknown_study:
        raise ValueError(f"Unknown convergence study key: {sorted(unknown_study)[0]}")

    base = _relative_to(path.parent, data["base"])
    if not base.is_file():
        raise FileNotFoundError(f"Base workflow input file not found: {base}")
    parameter = str(data["parameter"]).strip()
    if "." not in parameter:
        raise ValueError("Convergence parameter must be a dotted YAML path, e.g. kappa.mesh.")
    directory = _relative_to(path.parent, study.get("directory", "convergence-study"))
    component = str(study.get("component", "average")).lower()
    if component not in {*COMPONENTS, "average"}:
        raise ValueError("Convergence component must be x, y, z, or average.")
    tolerance = float(study.get("tolerance", 0.02))
    if not 0 < tolerance < 1:
        raise ValueError("Convergence tolerance must be between 0 and 1.")
    execute = study.get("execute", False)
    if not isinstance(execute, bool):
        raise ValueError("Convergence study.execute must be true or false.")
    return ConvergenceStudyConfig(
        config_path=path,
        base=base,
        parameter=parameter,
        values=tuple(data["values"]),
        directory=directory,
        execute=execute,
        temperature=float(study.get("temperature", 300.0)),
        component=component,
        tolerance=tolerance,
    )


def run_convergence_study(
    config: ConvergenceStudyConfig,
    run_case: Callable[[Path], int] | None = None,
):
    """Prepare all cases, optionally execute them, and analyze completed results."""
    cases = prepare_convergence_cases(config)
    print(f"[Convergence] Parameter: {config.parameter}")
    print(f"  - Study directory: {config.directory}")
    print(f"  - Cases: {len(cases)}")
    for case in cases:
        print(f"    {case.label}: {case.config_path}")

    failed = []
    if config.execute:
        if run_case is None:
            raise ValueError("An execution callback is required when study.execute is true.")
        for case in cases:
            print(f"\n[Convergence] Running {case.label}")
            if run_case(case.config_path) != 0:
                failed.append(case.label)
                print(f"  - Case failed: {case.label}")
    else:
        print("  - Execution disabled; set study.execute: true to run prepared cases.")

    summary = analyze_convergence(config, cases)
    summary["failed_cases"] = failed
    _write_yaml(config.directory / "convergence-summary.yaml", summary)
    return summary


def prepare_convergence_cases(config: ConvergenceStudyConfig):
    """Generate one validated workflow YAML candidate per sweep value."""
    with config.base.open("r", encoding="utf-8") as handle:
        base_data = yaml.safe_load(handle) or {}
    if not isinstance(base_data, dict):
        raise ValueError("Base workflow YAML must be a mapping.")

    cases = []
    inputs_dir = config.directory / "inputs"
    runs_dir = config.directory / "runs"
    inputs_dir.mkdir(parents=True, exist_ok=True)
    runs_dir.mkdir(parents=True, exist_ok=True)
    used_labels = set()
    for index, value in enumerate(config.values, start=1):
        label = _unique_label(value, used_labels)
        used_labels.add(label)
        case_data = yaml.safe_load(yaml.safe_dump(base_data, sort_keys=False))
        _set_dotted_value(case_data, config.parameter, value)
        result_dir = (runs_dir / label).resolve()
        output = _mapping_section(case_data, "output")
        _set_normalized_key(output, "result-dir", str(result_dir))
        case_path = (inputs_dir / f"{index:02d}-{label}.yaml").resolve()
        _write_yaml(case_path, case_data)
        cases.append(ConvergenceCase(index, value, label, case_path, result_dir))
    return cases


def analyze_convergence(config: ConvergenceStudyConfig, cases):
    """Extract target-temperature kappa and write CSV/PNG convergence outputs."""
    records = []
    for case in cases:
        record = {
            "index": case.index,
            "label": case.label,
            "value": case.value,
            "result_dir": str(case.result_dir),
            "status": "missing",
            "temperature": None,
            "kappa": None,
            "relative_error": None,
        }
        kappa_files = sorted(case.result_dir.glob("kappa-m*.hdf5"))
        if kappa_files:
            with h5py.File(kappa_files[-1], "r") as handle:
                temperatures = np.asarray(handle["temperature"][:], dtype=float)
                index = int(np.argmin(np.abs(temperatures - config.temperature)))
                tensor = np.asarray(handle["kappa"][index], dtype=float)
            value = (
                float(np.mean(tensor[:3]))
                if config.component == "average"
                else float(tensor[COMPONENTS[config.component]])
            )
            record.update(
                status="complete",
                temperature=float(temperatures[index]),
                kappa=value,
            )
        records.append(record)

    complete = [record for record in records if record["status"] == "complete"]
    converged_label = None
    analysis_complete = len(complete) == len(records)
    if complete:
        reference = complete[-1]["kappa"]
        denominator = max(abs(reference), np.finfo(float).eps)
        for record in complete:
            record["relative_error"] = abs(record["kappa"] - reference) / denominator
        if analysis_complete:
            for start in range(len(complete)):
                if all(
                    record["relative_error"] <= config.tolerance
                    for record in complete[start:]
                ):
                    converged_label = complete[start]["label"]
                    break
        _write_convergence_plot(config, complete)

    _write_convergence_csv(config.directory / "convergence.csv", records)
    print(f"  - Completed results: {len(complete)}/{len(records)}")
    if complete:
        print(f"  - Provisional reference case: {complete[-1]['label']}")
        if analysis_complete:
            print(f"  - First converged case: {converged_label}")
        else:
            print("  - Convergence decision deferred until every case is complete.")
        print(f"  - Table: {config.directory / 'convergence.csv'}")
        print(f"  - Plot: {config.directory / 'convergence.png'}")
    return {
        "parameter": config.parameter,
        "temperature_requested": config.temperature,
        "component": config.component,
        "tolerance": config.tolerance,
        "analysis_complete": analysis_complete,
        "reference_case": complete[-1]["label"] if complete else None,
        "first_converged_case": converged_label,
        "records": records,
    }


def _write_convergence_plot(config, records):
    labels = [record["label"] for record in records]
    values = [record["kappa"] for record in records]
    figure, axis = plt.subplots(figsize=(6.4, 4.8), constrained_layout=True)
    axis.plot(range(len(values)), values, marker="o")
    axis.set_xticks(range(len(labels)), labels, rotation=30, ha="right")
    axis.set_xlabel(config.parameter)
    axis.set_ylabel(r"Thermal conductivity (W m$^{-1}$ K$^{-1}$)")
    axis.grid(color="0.9", linewidth=0.9)
    figure.savefig(config.directory / "convergence.png", dpi=300)
    plt.close(figure)


def _write_convergence_csv(path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=records[0].keys())
        writer.writeheader()
        writer.writerows(records)


def _set_dotted_value(data, dotted_path, value):
    parts = dotted_path.split(".")
    current = data
    for part in parts[:-1]:
        key = _matching_key(current, part)
        if key is None or not isinstance(current[key], dict):
            raise ValueError(
                f"Convergence parameter path '{dotted_path}' does not exist "
                "as a mapping in the base YAML."
            )
        current = current[key]
    key = _matching_key(current, parts[-1])
    if key is None:
        raise ValueError(
            f"Convergence parameter path '{dotted_path}' does not exist in the base YAML."
        )
    current[key] = value


def _mapping_section(data, requested):
    key = _matching_key(data, requested)
    if key is None:
        key = requested
        data[key] = {}
    if not isinstance(data[key], dict):
        raise ValueError(f"Convergence path '{requested}' is not a YAML mapping.")
    return data[key]


def _set_normalized_key(data, requested, value):
    key = _matching_key(data, requested) or requested
    data[key] = value


def _matching_key(data, requested):
    normalized = str(requested).replace("_", "-")
    return next(
        (key for key in data if str(key).replace("_", "-") == normalized), None
    )


def _unique_label(value, used):
    raw = "x".join(str(item) for item in value) if isinstance(value, list) else str(value)
    label = re.sub(r"[^A-Za-z0-9_.-]+", "-", raw).strip("-.") or "value"
    if label not in used:
        return label
    suffix = 2
    while f"{label}-{suffix}" in used:
        suffix += 1
    return f"{label}-{suffix}"


def _relative_to(parent, value):
    path = Path(str(value)).expanduser()
    return path.resolve() if path.is_absolute() else (parent / path).resolve()


def _write_yaml(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(data, handle, sort_keys=False, allow_unicode=True)
