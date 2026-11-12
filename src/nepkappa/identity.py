"""Deterministic data and file identities shared by caches and provenance."""

import hashlib
import json
from pathlib import Path

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
