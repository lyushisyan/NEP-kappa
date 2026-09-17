"""Reproducibility metadata and stable input hashing for NEP-kappa."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
from importlib import metadata
import json
import os
from pathlib import Path
import platform
import socket
import subprocess
import sys

import yaml

from nepkappa import __version__


PROVENANCE_SCHEMA_VERSION = 1
DEPENDENCIES = (
    "ase",
    "calorine",
    "h5py",
    "hiphive",
    "matplotlib",
    "mace-torch",
    "numpy",
    "phonopy",
    "phono3py",
    "PyYAML",
    "seekpath",
    "trainstation",
)
ARTIFACT_NAMES = (
    "POSCAR_relaxed",
    "phono3py.yaml",
    "phono3py_disp.yaml",
    "fc2.hdf5",
    "fc3.hdf5",
    "FORCE_CONSTANTS_2ND",
    "FORCE_CONSTANTS_3RD",
    "FORCE_CONSTANTS_4TH",
    "hiphive_model.fcp",
    "e-v.dat",
    "volume-temperature.dat",
    "thermal_expansion.dat",
    "bulk_modulus-temperature.dat",
    "gibbs-temperature.dat",
    "Cp-temperature_polyfit.dat",
    "entropy-volume.dat",
    "Cv-volume.dat",
    "dsdv-temperature.dat",
    "gruneisen-temperature.dat",
    "helmholtz-volume.dat",
    "volume-temperature.pdf",
    "thermal_expansion.pdf",
    "bulk_modulus-temperature.pdf",
    "Cp-temperature_polyfit.pdf",
    "qha-summary.yaml",
)


def utc_now():
    """Return an ISO-8601 UTC timestamp."""
    return datetime.now(timezone.utc).isoformat()


def file_sha256(path, chunk_size=1024 * 1024):
    """Return the SHA256 digest of a file without loading it all into memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def file_identity(path):
    """Return stable provenance fields for a file, including its content hash."""
    path = Path(path).expanduser()
    identity = {"path": str(path.resolve(strict=False)), "exists": path.is_file()}
    if path.is_file():
        stat = path.stat()
        identity.update({"size": stat.st_size, "sha256": file_sha256(path)})
    return identity


def canonical_data(value):
    """Convert common configuration values into deterministic JSON data."""
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {
            str(key): canonical_data(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (list, tuple)):
        return [canonical_data(item) for item in value]
    if hasattr(value, "tolist"):
        return canonical_data(value.tolist())
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def data_sha256(value):
    """Return a deterministic SHA256 digest for JSON-compatible data."""
    encoded = json.dumps(
        canonical_data(value), sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def installed_versions():
    """Return versions of the scientific packages that define the result."""
    versions = {"nepkappa": __version__}
    for distribution in DEPENDENCIES:
        try:
            versions[distribution] = metadata.version(distribution)
        except metadata.PackageNotFoundError:
            versions[distribution] = None
    return versions


def git_state(cwd=None):
    """Return commit and dirty-state metadata when running from a Git checkout."""
    cwd = Path(cwd or os.getcwd())
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=cwd,
            text=True,
            capture_output=True,
            check=True,
            timeout=2,
        ).stdout.strip()
        dirty = bool(
            subprocess.run(
                ["git", "status", "--porcelain"],
                cwd=cwd,
                text=True,
                capture_output=True,
                check=True,
                timeout=2,
            ).stdout.strip()
        )
    except (FileNotFoundError, subprocess.SubprocessError):
        return None
    return {"commit": commit, "dirty": dirty}


class ProvenanceRecorder:
    """Maintain an atomic ``provenance.yaml`` throughout one CLI command."""

    def __init__(self, output_dir, command, config_path, config):
        self.output_dir = Path(output_dir).resolve()
        self.path = self.output_dir / "provenance.yaml"
        self.command = str(command)
        self.config_path = Path(config_path).resolve()
        self.config = canonical_data(vars(config))
        self.document = {}
        self.history_path = None

    def start(self):
        """Write initial metadata before the expensive workflow begins."""
        self.output_dir.mkdir(parents=True, exist_ok=True)
        inputs = {"config": file_identity(self.config_path)}
        for name in (
            "poscar",
            "nep_model",
            "potcar_path",
            "scph_initial_fc2",
            "scph_born",
            "scph_transport_fc3",
            "scph_transport_metadata",
        ):
            value = self.config.get(name)
            if value and Path(str(value)).is_file():
                inputs[name] = file_identity(value)
            elif value:
                inputs[name] = {
                    "path": str(Path(str(value)).resolve(strict=False)),
                    "exists": Path(str(value)).exists(),
                    "kind": "directory" if Path(str(value)).is_dir() else "missing",
                }
        calculator_name = str(self.config.get("calculator", "nep")).lower()
        force_commands = {
            "run",
            "relax",
            "fc",
            "fc2",
            "fc2fc3",
            "fc4",
            "qha",
            "scph",
        }
        if (
            calculator_name not in {"nep", "vasp"}
            and self.command in force_commands
        ):
            # Imported lazily to avoid a provenance/calculator module cycle.
            from nepkappa.calculators import calculator_backend_from_config

            backend = calculator_backend_from_config(self.config)
            inputs["external_calculator"] = backend.cache_inputs()
        preexisting = [self.output_dir / name for name in ARTIFACT_NAMES]
        preexisting.extend(sorted(self.output_dir.glob("kappa-m*.hdf5")))
        existing_inputs = {
            path.name: file_identity(path) for path in preexisting if path.is_file()
        }
        if existing_inputs:
            inputs["preexisting_results"] = existing_inputs

        started_at = utc_now()
        safe_command = "".join(
            character if character.isalnum() or character in "-_" else "-"
            for character in self.command
        )
        run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
        history_dir = self.output_dir / "provenance-history"
        history_dir.mkdir(parents=True, exist_ok=True)
        self.history_path = history_dir / f"{run_id}-{safe_command}.yaml"

        self.document = {
            "schema_version": PROVENANCE_SCHEMA_VERSION,
            "status": "running",
            "command": self.command,
            "invocation": ["nepkappa", self.command, str(self.config_path)],
            "started_at": started_at,
            "working_directory": str(Path.cwd()),
            "configuration": self.config,
            "inputs": inputs,
            "software": installed_versions(),
            "environment": {
                "python": sys.version.splitlines()[0],
                "platform": platform.platform(),
                "hostname": socket.gethostname(),
            },
            "git": git_state(),
            "artifacts": {},
        }
        self._write()
        return self.path

    def finish(self, status, *, error=None, stage_timings=None):
        """Finalize status, timings, and hashes of scientific result artifacts."""
        self.document["status"] = str(status)
        self.document["finished_at"] = utc_now()
        if error:
            self.document["error"] = str(error)
        if stage_timings:
            self.document["stage_timings_seconds"] = {
                label: float(seconds) for label, seconds in stage_timings
            }
        self.document["artifacts"] = self._artifact_identities()
        self._write()
        return self.path

    def _artifact_identities(self):
        paths = [self.output_dir / name for name in ARTIFACT_NAMES]
        paths.extend(sorted(self.output_dir.glob("kappa-m*.hdf5")))
        paths.extend(sorted((self.output_dir / "plots").glob("*.png")))
        for pattern in (
            "point.yaml",
            "POSCAR",
            "force_constants.hdf5",
            "phonopy_params.yaml",
            "thermal_properties.yaml",
        ):
            paths.extend(sorted((self.output_dir / "qha-points").glob(f"*/{pattern}")))
        scph_workdir = self.config.get("scph_workdir", "phonopy-sscha")
        paths.extend(sorted((self.output_dir / scph_workdir).rglob("*")))
        return {
            str(path.relative_to(self.output_dir)): file_identity(path)
            for path in paths
            if path.is_file()
        }

    def _write(self):
        content = yaml.safe_dump(self.document, sort_keys=False, allow_unicode=True)
        for path in (self.path, self.history_path):
            if path is None:
                continue
            temporary = path.with_suffix(".yaml.tmp")
            temporary.write_text(content, encoding="utf-8")
            temporary.replace(path)
