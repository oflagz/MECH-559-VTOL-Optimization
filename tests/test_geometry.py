import pytest
from src.rotor.geometry import rotor_geometry

def test_screening_geometry():
    area, solidity = rotor_geometry(0.25,0.030,3)
    assert area == pytest.approx(0.19635, abs=1e-5)
    assert solidity == pytest.approx(0.11459, abs=1e-5)
