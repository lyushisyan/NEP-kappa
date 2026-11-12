"""Exercise the reference and compiled reductions in separate CI jobs."""
import os

import pytest


@pytest.fixture(scope='session', autouse=True)
def collision_backend():
    from nepkappa import tdbte_accumulate as reductions

    requested = os.environ.get('NEPKAPPA_TEST_COLLISION_BACKEND')
    if requested is None:
        yield
        return
    if requested not in {'numpy', 'numba'}:
        pytest.fail(f'Unknown test collision backend: {requested}')
    if requested == 'numba':
        assert reductions.accumulate is not None, 'CI requested Numba, but it is unavailable'
        yield
        return
    # HiPhive also depends on Numba. Keep that dependency available while forcing
    # TD-BTE through its NumPy reference path, without a production-only switch.
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(reductions, 'accumulate', None)
        patch.setattr(reductions, 'accumulate_batch', None)
        yield
