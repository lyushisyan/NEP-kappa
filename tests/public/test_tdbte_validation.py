"""Small independent checks of the opt-in physical-audit mathematics.

These synthetic checks do not establish physical scattering-rate accuracy.
"""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

from nepkappa.tdbte import ShellKernel
from nepkappa.tdbte_storage import ChunkWriter, ChunkedKernel

_path = Path(__file__).resolve().parents[2] / 'tools/validate_tdbte_kernel.py'
if not _path.is_file():
    pytest.skip('Optional local TD-BTE audit tool is unavailable', allow_module_level=True)
_spec = importlib.util.spec_from_file_location('tdbte_validation', _path)
audit = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(audit)


def test_repeated_daughters_are_combined_in_occupation_jacobian():
    k = ShellKernel(np.array([2., 1.]), np.array([0]), np.array([[1]*4]),
                    np.array([[1]*4]), np.full((1, 4), .25), np.array([.2]),
                    np.zeros(2, dtype=bool))
    n, b, c = audit.equilibrium_event_data(k, 300.)
    reverse = k.coefficients * (1 + n[k.parent]) * b * c
    B = audit.incidence(k)
    np.testing.assert_array_equal(B.toarray(), [[-1], [2]])
    diagonal = (B.power(2) @ reverse) / (n * (1 + n))
    eps = n * 1e-5
    measured = []
    for i in range(2):
        delta = np.zeros(2)
        delta[i] = eps[i]
        measured.append(-(k.collision(n+delta)[0][i] - k.collision(n-delta)[0][i]) / (2*eps[i]))
    np.testing.assert_allclose(diagonal, measured, rtol=1e-8)
    # A tagged daughter sees one leg; the diagonal includes both cross terms.
    tagged_daughter = 2 * .2 * (n[1] - n[0])
    np.testing.assert_allclose(diagonal[1], 2*tagged_daughter)


@pytest.mark.parametrize('level', [0, 1, 2])
def test_surface_refinement_preserves_polynomial_moments(level):
    base = np.full((3, 3), 1/6) + np.eye(3)/2
    nodes = audit.refinement_transform(level) @ base
    assert np.min(nodes) >= -1e-14
    np.testing.assert_allclose(nodes.sum(axis=1), 1)
    # Mean x^2 and xy on a unit triangle, independently integrated analytically.
    np.testing.assert_allclose(np.mean(nodes[:, 1]**2), 1/6)
    np.testing.assert_allclose(np.mean(nodes[:, 1]*nodes[:, 2]), 1/12)


def test_surface_triangles_split_across_event_chunks(tmp_path):
    k = ShellKernel(np.array([5., 1., 2., 3., 4.]), np.zeros(3, dtype=int),
        np.tile([1, 2, 1, 2], (3, 1)), np.tile([3, 4, 3, 4], (3, 1)),
        np.full((3, 4), .25), np.full(3, .2/3), np.zeros(5, dtype=bool))
    writer = ChunkWriter(tmp_path/'kernel', k.frequency, k.zero, {})
    for start, stop in ((0, 2), (2, 3)):
        writer.append(ShellKernel(k.frequency, k.parent[start:stop], k.daughter1[start:stop],
            k.daughter2[start:stop], k.weights[start:stop], k.coefficients[start:stop], k.zero))
    kernel = ChunkedKernel(writer.finish())
    _, b, c = audit.equilibrium_event_data(k, 300.)
    expected = np.zeros(5)
    expected[0] = np.sum(k.coefficients * (1 + b + c))
    for level in (0, 1, 2):
        np.testing.assert_allclose(audit.shell_rates(kernel, 300., level), expected, rtol=1e-14)


def test_incomplete_surface_triangle_fails(tmp_path):
    k = ShellKernel(np.array([3., 1., 2.]), np.array([0]), np.array([[1]*4]),
        np.array([[2]*4]), np.full((1, 4), .25), np.array([.2]), np.zeros(3, dtype=bool))
    writer = ChunkWriter(tmp_path/'kernel', k.frequency, k.zero, {})
    writer.append(k)
    with pytest.raises(ValueError, match='Incomplete'):
        list(audit.triangle_groups(ChunkedKernel(writer.finish())))
