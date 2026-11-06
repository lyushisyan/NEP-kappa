"""Positive surface quadrature for linear three-phonon energy shells.

This is a weak, entropy-interpolated discretization, not phono3py's native
tetrahedron collision operator. Coordinates are in mesh-cell units.
"""
import numpy as np


def interval_quadrature(vertices, detuning):
    """Integrate delta(linear detuning) on a 1D cell in mesh-cell units.

    A root at a shared endpoint receives half weight from either cell. Flat
    resonances are singular, not zero scattering, and must not be discarded.
    The native vertex already contains 1/Nq, so no additional mesh factor is
    included here. Near-tangent roots require independent mesh convergence.
    """
    vertices, d = np.asarray(vertices, float), np.asarray(detuning, float)
    if (vertices.shape != (2, 3) or d.shape != (2,)
            or not np.isfinite(vertices).all() or not np.isfinite(d).all()):
        raise ValueError('Invalid interval vertices/detuning')
    length = np.linalg.norm(vertices[1] - vertices[0])
    if length <= 0:
        raise ValueError('Degenerate interval')
    slope = d[1] - d[0]
    if abs(slope) < 1e-13:
        if np.max(abs(d)) < 1e-12:
            raise ValueError('Flat resonant interval: undefined delta integral')
        return []
    if min(d) > 0 or max(d) < 0:
        return []
    t = -d[0] / slope
    boundary = .5 if d[0] == 0 or d[1] == 0 else 1.
    return [(np.array([1 - t, t]), boundary * length / abs(slope))]


def section_quadrature(vertices, detuning, rule='symmetric-degree2'):
    """Integrate delta(linear detuning) over a tetrahedron.

    Return (four barycentric weights, surface weight) at each quadrature node.
    A resonant shared face has half weight; a flat resonant volume is singular.
    """
    if rule not in ('centroid', 'symmetric-centroid', 'symmetric-degree2'):
        raise ValueError('Unknown section quadrature rule')
    vertices = np.asarray(vertices, float)
    d = np.asarray(detuning, float)
    gradient = np.linalg.solve(vertices[1:] - vertices[0], d[1:] - d[0])
    norm = np.linalg.norm(gradient)
    if norm < 1e-13:
        if np.max(abs(d)) < 1e-12:
            raise ValueError('Flat resonant tetrahedron: undefined delta integral')
        return []
    points = []
    eye = np.eye(4)
    for i in range(4):
        if d[i] == 0:
            points.append(eye[i])
        for j in range(i + 1, 4):
            if d[i] * d[j] < 0:
                t = d[i] / (d[i] - d[j])
                points.append((1 - t) * eye[i] + t * eye[j])
    if len(points) < 3:
        return []
    bary = np.asarray(points)
    xyz = bary @ vertices
    center = xyz.mean(axis=0)
    first = xyz[0] - center
    first /= np.linalg.norm(first)
    second = np.cross(gradient / norm, first)
    order = np.argsort(np.arctan2((xyz - center) @ second, (xyz - center) @ first))
    bary = bary[order]
    boundary_weight = .5 if np.count_nonzero(d == 0) == 3 else 1.
    if rule == 'centroid':
        triangles = [bary[[0, i, i + 1]] for i in range(1, len(bary) - 1)]
    else:
        triangles = [np.array([bary.mean(axis=0), bary[i], bary[(i + 1) % len(bary)]])
                     for i in range(len(bary))]
    nodes = (np.full((3, 3), 1 / 6) + np.eye(3) / 2
             if rule == 'symmetric-degree2' else np.array([[1 / 3] * 3]))
    result = []
    for triangle in triangles:
        v = triangle @ vertices
        area = np.linalg.norm(np.cross(v[1] - v[0], v[2] - v[0])) / 2
        if area > 0:
            for node in nodes:
                result.append((node @ triangle, boundary_weight * area / norm / len(nodes)))
    return result
