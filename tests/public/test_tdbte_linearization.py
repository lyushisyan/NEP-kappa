"""Independent event-matrix proof checks, including nontrivial interpolation."""
from dataclasses import replace

import numpy as np
from nepkappa.tdbte import ShellKernel


def interpolated_kernel():
    return ShellKernel(
        np.array([5., 1., 2., 3., 4.]), np.array([0]),
        np.array([[1, 2, 1, 2]]), np.array([[3, 4, 3, 4]]),
        np.full((1, 4), .25), np.array([.2]), np.zeros(5, dtype=bool))


def test_interpolated_jacobian_from_stoichiometry():
    k = interpolated_kernel()
    k.validate()
    n = k.bose(300)
    eta = np.log1p(1/n)
    b = 1 / np.expm1((eta[1]+eta[2])/2)
    c = 1 / np.expm1((eta[3]+eta[4])/2)
    B = np.array([-1., .5, .5, .5, .5])
    R = .2*(1+n[0])*b*c
    D = n*(1+n)
    J = -R*np.outer(B, B)/D[None, :]
    numerical = np.empty_like(J)
    for i in range(len(n)):
        delta = np.zeros_like(n)
        delta[i] = n[i]*1e-5
        numerical[:, i] = (k.collision(n+delta)[0]-k.collision(n-delta)[0])/(2*delta[i])
    np.testing.assert_allclose(numerical, J, rtol=1e-8, atol=1e-11)
    np.testing.assert_allclose(k.frequency @ J, 0, atol=1e-13)
    np.testing.assert_allclose(J @ (D*k.frequency), 0, atol=1e-13)
    L = -J*np.sqrt(D[None, :]/D[:, None])
    np.testing.assert_allclose(L, L.T, atol=1e-13)
    assert np.linalg.eigvalsh(L).min() > -1e-13


def test_interpolated_entropy_identity_and_rate_scaling():
    k = interpolated_kernel()
    B = np.array([-1., .5, .5, .5, .5])
    for n in np.random.default_rng(55).uniform(.1, 3, (20, 5)):
        rhs, flux = k.collision(n)
        eta = np.log1p(1/n)
        np.testing.assert_allclose(eta @ rhs, flux[0]*(eta @ B), atol=1e-13)
        assert eta @ rhs >= -1e-13
        np.testing.assert_allclose(k.frequency @ rhs, 0, atol=1e-13)
        scaled = replace(k, coefficients=10*k.coefficients)
        np.testing.assert_allclose(scaled.collision(n)[0], 10*rhs, rtol=1e-13)
