import pytest
from src.rotor.parameters import Params
from src.rotor.hover import hover_model

def test_proposal_hover_screening_case():
    p = Params()
    h = hover_model(0.25,0.030,300.0,3,p)
    assert h["T_req"] == pytest.approx(12.876, abs=0.01)
    assert h["sigma"] == pytest.approx(0.1146, abs=1e-4)
    assert h["Pi"] == pytest.approx(66.7, abs=0.2)
    assert h["P_sh"] == pytest.approx(94.2, abs=0.3)
    assert h["FM"] == pytest.approx(0.708, abs=0.003)
    assert p.Nr*h["P_el"] == pytest.approx(443.0, abs=2.0)
