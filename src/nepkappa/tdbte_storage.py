"""Versioned, checksummed energy-shell chunks with bounded event memory."""
from contextlib import contextmanager
import hashlib
import json
from pathlib import Path

import numpy as np

from .tdbte import ShellKernel

FORMAT = 'nepkappa-tdbte-shell-v1'
EVENT_FIELDS = ('parent', 'daughter1', 'daughter2', 'weights', 'coefficients')


def sha256(path):
    """Hash large artifacts without loading the file into memory."""
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path, data):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.partial')
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def _local_path(root, name):
    if not isinstance(name, str) or Path(name).is_absolute():
        raise ValueError('Chunk paths must be relative to their manifest')
    path = (root / name).resolve()
    if root.resolve() not in path.parents:
        raise ValueError('Chunk path escapes its manifest directory')
    return path


class ChunkWriter:
    """Write validated chunks; publish a manifest only after completion."""

    def __init__(self, directory, frequency, zero, metadata):
        self.root = Path(directory)
        self.root.mkdir(parents=True, exist_ok=False)
        self.frequency = np.asarray(frequency)
        self.zero = np.asarray(zero)
        self.metadata = dict(metadata)
        self.chunks = []
        np.savez(self.root / 'modes.npz', frequency_THz=frequency, zero_modes=zero)

    def append(self, kernel):
        kernel.validate()
        if not np.array_equal(kernel.frequency, self.frequency) or not np.array_equal(kernel.zero, self.zero):
            raise ValueError('Chunk modes differ from the manifest modes')
        directory = self.root / f'events-{len(self.chunks):06d}'
        directory.mkdir()
        hashes = {}
        for name in EVENT_FIELDS:
            path = directory / f'{name}.npy'
            np.save(path, getattr(kernel, name), allow_pickle=False)
            hashes[name] = sha256(path)
        self.chunks.append({'directory': directory.name, 'events': len(kernel.parent),
                            'sha256': hashes})

    def finish(self):
        if not self.chunks:
            raise ValueError('No nonzero scattering events; refusing an empty kernel')
        manifest = dict(self.metadata)
        manifest.update(format=FORMAT, complete=True, chunks=self.chunks,
                        surface_quadrature_points=sum(c['events'] for c in self.chunks),
                        modes={'file': 'modes.npz', 'sha256': sha256(self.root / 'modes.npz')})
        path = self.root / 'manifest.json'
        atomic_json(path, manifest)
        return path


class ChunkedKernel:
    """Read one memory-mapped event block at a time, never an Nmode² matrix.

    Checksums and shell identities are checked once, before integration. Kernel
    files must remain immutable during a run. Mappings are closed after each
    block so a large manifest does not exhaust file descriptors or pin all data.
    """

    def __init__(self, path):
        from .tdbte_accumulate import accumulate
        self.accumulate = accumulate
        self.backend = 'numba' if accumulate is not None else 'numpy'
        self.path = Path(path).resolve()
        self.root = self.path.parent
        self.metadata = json.loads(self.path.read_text())
        meta = self.metadata
        if meta.get('format') != FORMAT or meta.get('complete') is not True:
            raise ValueError('Unsupported or incomplete TD-BTE chunk manifest')
        self.chunks = meta.get('chunks')
        if not isinstance(self.chunks, list) or not self.chunks:
            raise ValueError('No completed chunks in manifest')
        modes = _local_path(self.root, meta['modes']['file'])
        if sha256(modes) != meta['modes']['sha256']:
            raise ValueError('Mode checksum mismatch')
        with np.load(modes, allow_pickle=False) as data:
            self.frequency = data['frequency_THz'].copy()
            self.zero = data['zero_modes'].copy()
        count, seen = 0, set()
        for entry in self.chunks:
            directory = _local_path(self.root, entry['directory'])
            if directory in seen:
                raise ValueError('Duplicate chunk in manifest')
            seen.add(directory)
            for name in EVENT_FIELDS:
                path = _local_path(self.root, str(Path(entry['directory']) / f'{name}.npy'))
                if sha256(path) != entry['sha256'][name]:
                    raise ValueError(f'Chunk checksum mismatch: {path}')
            with self._read(entry) as kernel:
                kernel.validate()
                if len(kernel.parent) != entry['events']:
                    raise ValueError('Chunk event count mismatch')
                if self.accumulate is not None:
                    self._check_compiled(kernel)
                count += len(kernel.parent)
        if count != meta.get('surface_quadrature_points'):
            raise ValueError('Incomplete event collection')
        print(f'TD-BTE kernel validated: {count} events in {len(self.chunks)} chunks; '
              f'{self.backend} reductions', flush=True)

    @contextmanager
    def _read(self, entry):
        arrays = []
        try:
            for name in EVENT_FIELDS:
                path = _local_path(self.root, str(Path(entry['directory']) / f'{name}.npy'))
                arrays.append(np.load(path, mmap_mode='r', allow_pickle=False))
            yield ShellKernel(self.frequency, *arrays, self.zero)
        finally:
            for array in arrays:
                array._mmap.close()

    def bose(self, temperature):
        return ShellKernel.bose(self, temperature)

    def _check_compiled(self, kernel):
        # Verify the installed compiler against NumPy on each block before use.
        rng = np.random.default_rng(20260927)
        n = kernel.bose(300) * (1 + .01 * rng.uniform(-1, 1, len(self.frequency)))
        eta = np.zeros_like(n)
        eta[~self.zero] = np.log1p(n[~self.zero]) - np.log(n[~self.zero])
        rhs = np.zeros_like(n)
        self.accumulate(n, eta, kernel.parent, kernel.daughter1, kernel.daughter2,
                        kernel.weights, kernel.coefficients, rhs)
        rhs[self.zero] = 0
        reference = kernel.collision(n)[0]
        error = np.linalg.norm(rhs - reference) / max(np.linalg.norm(reference), 1e-30)
        if not np.isfinite(error) or error > 1e-10:
            raise ValueError(f'Compiled/NumPy collision disagreement: {error}')

    def collision(self, occupation):
        n = np.asarray(occupation, dtype=float)
        if n.shape != self.frequency.shape or not np.isfinite(n).all() or np.any(n[~self.zero] <= 0):
            raise ValueError('Active occupations must be finite and positive; no clipping is applied')
        eta = np.zeros_like(n)
        eta[~self.zero] = np.log1p(n[~self.zero]) - np.log(n[~self.zero])
        rhs = np.zeros_like(self.frequency, dtype=float)
        for entry in self.chunks:
            with self._read(entry) as kernel:
                if self.accumulate is None:
                    rhs += kernel.collision(n)[0]
                else:
                    self.accumulate(n, eta, kernel.parent, kernel.daughter1, kernel.daughter2,
                                    kernel.weights, kernel.coefficients, rhs)
        rhs[self.zero] = 0
        if not np.isfinite(rhs).all():
            raise ValueError('Nonfinite accumulated collision derivative')
        return rhs, None
