"""A transport-equivalent odd kernel is not a full population generator."""
import numpy as np


def test_transport_rearrangement_preserves_odd_but_not_even_dynamics():
    # Exact resonant event and its time-reversed q partner; E=(3,1,2).
    energy = np.tile([3., 1., 2.], 2)
    n = 1 / np.expm1(energy / 5.)
    root_d = np.sqrt(n * (1+n))
    B = np.zeros((6, 2))
    B[:3, 0] = [-1, 1, 1]
    B[3:, 1] = [-1, 1, 1]
    weighted = B / root_d[:, None]
    full = weighted @ weighted.T
    inverse = np.eye(6)[[3, 4, 5, 0, 1, 2]]
    # Replace negative cross terms with positive terms at inverted columns.
    negative = np.minimum(full, 0)
    transport = full - negative - negative @ inverse
    odd = (np.eye(6)-inverse)/2
    even = (np.eye(6)+inverse)/2
    np.testing.assert_allclose((transport-full) @ odd, 0, atol=1e-14)
    assert np.linalg.norm((transport-full) @ even) > .1
    thermal = root_d*energy
    np.testing.assert_allclose(full @ thermal, 0, atol=1e-14)
    assert np.linalg.norm(transport @ thermal) > .1
    assert np.linalg.eigvalsh(full).min() > -1e-13
