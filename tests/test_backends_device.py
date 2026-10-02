import pytest

from mechbridge.backends import PySCFCalculator


def test_pyscf_calculator_device_is_explicit():
    assert PySCFCalculator(0, 1).device == 'cpu'
    with pytest.raises(ValueError, match='Unknown device'):
        PySCFCalculator(0, 1, device='tpu')
