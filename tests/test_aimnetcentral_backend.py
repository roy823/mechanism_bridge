import numpy as np
import pytest
from mechbridge.aimnetcentral_backend import AIMNetCentralPotential


def test_model_family_table_is_explicit():
    from mechbridge.aimnetcentral_backend import FAMILIES
    assert set(FAMILIES)=={'aimnet2','aimnet2-2025','aimnet2-nse','aimnet2-rxn'}
    assert FAMILIES['aimnet2-nse']['open_shell']
    assert not FAMILIES['aimnet2-rxn']['open_shell']


def test_system_validation_without_loading_weights():
    model=object.__new__(AIMNetCentralPotential)
    model.family='aimnet2';model.supported_species={1,6,7,8,9,15,16,17,35,53};model.open_shell=False
    model.validate_system([6,17,1],0,1)
    with pytest.raises(ValueError,match='outside'):model.validate_system([6,26],0,1)
    with pytest.raises(ValueError,match='singlets'):model.validate_system([6],0,2)
    model.family='aimnet2-rxn';model.supported_species={1,6,7,8}
    with pytest.raises(ValueError,match='net-neutral'):model.validate_system([6],1,1)
    model.family='aimnet2-nse';model.open_shell=True
    model.validate_system([6,1,1,1],0,2)
