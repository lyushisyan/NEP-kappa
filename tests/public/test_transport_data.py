import numpy as np
import pytest
from nepkappa.transport_data import read_conductivity_components


@pytest.mark.parametrize('change', ['missing', 'mismatch', 'shape', 'nonfinite'])
def test_rejects_incomplete_or_inconsistent_wigner_components(change):
    data = {'kappa': np.ones((1, 6))*3, 'kappa_intra_exact': np.ones((1, 6))*2,
            'kappa_inter': np.ones((1, 6))}
    if change == 'missing':
        del data['kappa_inter']
    elif change == 'mismatch':
        data['kappa'] *= 10
    elif change == 'shape':
        data['kappa_inter'] = np.ones((2, 6))
    else:
        data['kappa_inter'][0, 0] = np.nan
    with pytest.raises(ValueError):
        read_conductivity_components(data)
