from src.rotor.validation import summary

def test_validation_summary():
    data = [{"converged":True,"thrust_error_n":-2.0,"power_error_w":4.0},
            {"converged":False,"thrust_error_n":1.0,"power_error_w":-2.0}]
    s = summary(data)
    assert s["n"] == 2
    assert s["converged"] == 1
    assert s["thrust_mae_n"] == 1.5
    assert s["power_mae_w"] == 3.0
