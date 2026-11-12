"""Analytic surface moments for the packaged energy-shell builder."""
import itertools

import numpy as np
import pytest

from nepkappa.tdbte_quadrature import section_quadrature
VERTICES = np.array([[0., 0, 0], [1., 0, 0], [0., 1, 0], [0., 0, 1]])


@pytest.mark.parametrize('height', [.2, .5, .8])
def test_exact_plane_section_polynomial_moments(height):
    # delta(z-h) in the unit simplex: right triangle of side a=1-h.
    a = 1-height
    points = section_quadrature(VERTICES, VERTICES[:, 2]-height, 'symmetric-degree2')
    xyz = np.array([w @ VERTICES for w, _ in points])
    weights = np.array([v for _, v in points])
    measured = weights @ np.array([np.ones(len(xyz)), xyz[:, 0],
                                  xyz[:, 0]**2, xyz[:, 0]*xyz[:, 1]]).T
    np.testing.assert_allclose(measured, [a*a/2, a**3/6, a**4/12, a**4/24], rtol=1e-13)
    assert np.all(weights > 0)


def test_nonlinear_integral_independent_of_vertex_order():
    d = np.array([-1., -2., 3., 4.])  # quadrilateral section
    values = []
    for permutation in itertools.permutations(range(4)):
        p = np.array(permutation)
        points = section_quadrature(VERTICES[p], d[p], 'symmetric-degree2')
        values.append(sum(v*np.exp((w @ VERTICES[p]) @ [1., 2., 3.]) for w, v in points))
    np.testing.assert_allclose(values, values[0], rtol=1e-14)


def test_boundary_half_weight_and_flat_rejection():
    points = section_quadrature(VERTICES, VERTICES[:, 2], 'symmetric-degree2')
    np.testing.assert_allclose(sum(v for _, v in points), .25)
    with pytest.raises(ValueError, match='Flat resonant'):
        section_quadrature(VERTICES, np.zeros(4), 'symmetric-degree2')
