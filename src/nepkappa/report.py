"""Compact, reproducible summaries for completed NEP-kappa result trees."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import yaml


ARTIFACT_NAMES = {
    "fc2.hdf5",
    "fc3.hdf5",
    "FORCE_CONSTANTS_2ND",
    "FORCE_CONSTANTS_3RD",
    "FORCE_CONSTANTS_4TH",
    "phono3py_disp.yaml",
    "qha-summary.yaml",
    "sscha-summary.yaml",
    "qha-sscha-summary.yaml",
    "fourphonon-summary.yaml",
}


def resolve_result_directory(source):
    """Resolve a result directory directly or from a workflow YAML file."""
    source = Path(source).expanduser()
    if source.is_dir():
        return source.resolve()
    if not source.is_file():
        raise FileNotFoundError(f"Result directory or YAML input not found: {source}")
    if source.suffix.lower() not in {".yaml", ".yml"}:
        raise ValueError("Report input must be a result directory or a YAML file.")
    data = yaml.safe_load(source.read_text(encoding="utf-8")) or {}
    output = data.get("output", {})
    result = output.get("result-dir", output.get("result_dir"))
    if not result:
        raise ValueError(f"No output.result-dir is defined in {source}")
    return Path(result).expanduser().resolve()


def generate_report(result_dir):
    """Write ``report.yaml`` and ``report.md`` below one result directory."""
    root = Path(result_dir).expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"Result directory not found: {root}")
    warnings = []
    report = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "result_directory": str(root),
        "provenance": _read_yaml_if_present(root / "provenance.yaml", warnings),
        "qha": _read_qha(root, warnings),
        "sscha": _read_sscha(root, warnings),
        "transport": _read_transport(root, warnings),
        "artifacts": _find_artifacts(root),
        "figures": _find_figures(root),
        "warnings": warnings,
    }
    yaml_path = root / "report.yaml"
    markdown_path = root / "report.md"
    yaml_path.write_text(
        yaml.safe_dump(report, sort_keys=False, allow_unicode=True), encoding="utf-8"
    )
    markdown_path.write_text(_render_markdown(report), encoding="utf-8")
    return (yaml_path, markdown_path)


def _read_yaml_if_present(path, warnings):
    if not path.is_file():
        return None
    try:
        return yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except Exception as exc:
        warnings.append(f"Could not read {path.name}: {exc}")
        return None


def _read_qha(root, warnings):
    path = root / "volume-temperature.dat"
    if not path.is_file():
        return None
    try:
        values = np.loadtxt(path, comments="#", ndmin=2)
        if values.shape[1] < 2:
            raise ValueError("expected at least two columns")
        records = [
            {"temperature_K": float(row[0]), "volume_A3": float(row[1])}
            for row in values
        ]
        return {
            "volume_temperature_file": str(path.relative_to(root)),
            "points": records,
        }
    except Exception as exc:
        warnings.append(f"Could not summarize QHA volumes: {exc}")
        return None


def _read_sscha(root, warnings):
    records = []
    for path in sorted(root.rglob("summary.yaml")):
        data = _read_yaml_if_present(path, warnings)
        if not isinstance(data, dict) or "temperature" not in data:
            continue
        records.append(
            {
                "file": str(path.relative_to(root)),
                "status": data.get("status"),
                "temperature_K": _finite_or_none(data.get("temperature")),
                "free_energy_eV_per_primitive_cell": _finite_or_none(
                    data.get("free_energy_ev_per_primitive_cell")
                ),
                "max_departure_sigma": _finite_or_none(
                    data.get("maximum_kept_free_energy_departure_sigma")
                ),
                "max_translational_drift_eV_per_A2": _finite_or_none(
                    data.get("maximum_translational_drift_ev_per_angstrom2")
                ),
            }
        )
    return records


def _read_transport(root, warnings):
    records = []
    for path in sorted(root.rglob("kappa-m*.hdf5")):
        try:
            with h5py.File(path, "r") as handle:
                if "temperature" not in handle or "kappa" not in handle:
                    continue
                temperatures = np.atleast_1d(np.asarray(handle["temperature"], dtype=float))
                kappa = np.asarray(handle["kappa"], dtype=float)
            if kappa.ndim == 1:
                kappa = kappa.reshape(1, -1)
            for temperature, tensor in zip(temperatures, kappa):
                if tensor.size < 3:
                    continue
                diagonal = np.asarray(tensor[:3], dtype=float)
                records.append(
                    {
                        "file": str(path.relative_to(root)),
                        "method": "3ph-phono3py",
                        "temperature_K": float(temperature),
                        "kxx_W_mK": float(diagonal[0]),
                        "kyy_W_mK": float(diagonal[1]),
                        "kzz_W_mK": float(diagonal[2]),
                        "kavg_W_mK": float(np.mean(diagonal)),
                    }
                )
        except Exception as exc:
            warnings.append(f"Could not summarize {path.relative_to(root)}: {exc}")
    records.extend(_read_fourphonon_transport(root, warnings))
    return records


def _read_fourphonon_transport(root, warnings):
    """Read normalized 3ph+4ph tensors from FourPhonon summaries."""
    records = []
    for path in sorted(root.rglob("fourphonon-summary.yaml")):
        data = _read_yaml_if_present(path, warnings)
        if not isinstance(data, dict):
            continue
        primary = data.get("primary")
        solution = (data.get("solutions") or {}).get(primary)
        if not isinstance(solution, dict):
            continue
        temperatures = solution.get("temperatures") or []
        diagonals = solution.get("diagonal") or []
        if len(temperatures) != len(diagonals):
            warnings.append(
                f"Could not summarize {path.relative_to(root)}: "
                "temperature and tensor lengths differ"
            )
            continue
        for temperature, diagonal in zip(temperatures, diagonals):
            try:
                tensor = np.asarray(diagonal, dtype=float)
                if tensor.size < 3 or not np.all(np.isfinite(tensor[:3])):
                    raise ValueError("invalid diagonal tensor")
                temperature = float(temperature)
                if not np.isfinite(temperature):
                    raise ValueError("invalid temperature")
            except (TypeError, ValueError) as exc:
                warnings.append(
                    f"Could not summarize {path.relative_to(root)}: {exc}"
                )
                continue
            records.append(
                {
                    "file": str(path.relative_to(root)),
                    "method": f"3+4ph-{primary}",
                    "temperature_K": temperature,
                    "kxx_W_mK": float(tensor[0]),
                    "kyy_W_mK": float(tensor[1]),
                    "kzz_W_mK": float(tensor[2]),
                    "kavg_W_mK": float(np.mean(tensor[:3])),
                }
            )
    return records


def _find_artifacts(root):
    paths = []
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        if path.name in ARTIFACT_NAMES or path.name.startswith("kappa-m"):
            paths.append(str(path.relative_to(root)))
    return sorted(paths)


def _find_figures(root):
    return sorted(
        str(path.relative_to(root))
        for path in root.rglob("*")
        if path.is_file() and path.suffix.lower() in {".png", ".pdf", ".svg"}
    )


def _finite_or_none(value):
    if value is None:
        return None
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if np.isfinite(value) else None


def _render_markdown(report):
    lines = [
        "# NEP-kappa report",
        "",
        f"- Result directory: `{report['result_directory']}`",
        f"- Generated: {report['generated_at']}",
    ]
    provenance = report.get("provenance") or {}
    if provenance:
        lines.extend(
            [
                f"- Run status: {provenance.get('status', 'unknown')}",
                f"- Command: `{provenance.get('command', 'unknown')}`",
            ]
        )
    qha = report.get("qha")
    if qha:
        lines.extend(["", "## QHA equilibrium volume", "", "| T (K) | V (Å³) |", "|---:|---:|"])
        lines.extend(
            f"| {row['temperature_K']:.6g} | {row['volume_A3']:.8g} |"
            for row in qha["points"]
        )
    if report["sscha"]:
        lines.extend(
            [
                "",
                "## SSCHA",
                "",
                "| T (K) | status | max drift (eV/Å²) | result |",
                "|---:|:---|---:|:---|",
            ]
        )
        for row in report["sscha"]:
            drift = row["max_translational_drift_eV_per_A2"]
            drift_text = "—" if drift is None else f"{drift:.3e}"
            lines.append(
                f"| {row['temperature_K']:.6g} | {row['status'] or 'unknown'} | "
                f"{drift_text} | `{row['file']}` |"
            )
    if report["transport"]:
        lines.extend(
            [
                "",
                "## Thermal conductivity",
                "",
                "| method | T (K) | κxx | κyy | κzz | κavg (W m⁻¹ K⁻¹) | result |",
                "|:---|---:|---:|---:|---:|---:|:---|",
            ]
        )
        for row in report["transport"]:
            lines.append(
                f"| {row.get('method', 'unknown')} | {row['temperature_K']:.6g} | "
                f"{row['kxx_W_mK']:.6g} | "
                f"{row['kyy_W_mK']:.6g} | {row['kzz_W_mK']:.6g} | "
                f"{row['kavg_W_mK']:.6g} | `{row['file']}` |"
            )
    lines.extend(
        [
            "",
            "## Files",
            "",
            f"- Key artifacts: {len(report['artifacts'])}",
            f"- Figures: {len(report['figures'])}",
        ]
    )
    if report["warnings"]:
        lines.extend(["", "## Warnings", ""])
        lines.extend(f"- {warning}" for warning in report["warnings"])
    return "\n".join(lines) + "\n"
