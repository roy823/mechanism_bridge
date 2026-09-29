import numpy as np
import pytest
from mechbridge.kinetics import tst_unimolecular,rate_matrix

def test_free_energy_sensitivity_and_detailed_balance():
    ratio=tst_unimolecular(10,298)/tst_unimolecular(11,298)
    assert 5.3<ratio<5.5
    base={'barrier_kind':'activation_gibbs_free_energy','molecularity':1,'connectivity_verified':True}
    q=rate_matrix(['a','b'],[{**base,'from':'a','to':'b','delta_g_kcal_mol':10},
                            {**base,'from':'b','to':'a','delta_g_kcal_mol':11}],298)
    assert np.allclose(q.sum(axis=1),0)
    assert q[0,1]/q[1,0]==pytest.approx(ratio)

def test_electronic_barrier_cannot_be_silently_used_for_rates():
    with pytest.raises(ValueError):
        rate_matrix(['a','b'],[{'from':'a','to':'b','barrier_kind':'electronic_energy','molecularity':1,
                                'delta_g_kcal_mol':10,'connectivity_verified':True}],298)
