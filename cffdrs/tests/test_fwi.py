import pytest
from cffdrs.fwi import ffmc


def test_ffmc(fwi_test_data):
    """Test ffmc function using sequential test data from fwi_test_data.csv"""
    initial_ffmc = 85.0
    current_ffmc = initial_ffmc

    for row in fwi_test_data:
        temp = row["TEMP"]
        rh = row["RH"]
        ws = row["WS"]
        prec = row["PREC"]
        expected_ffmc = row["FFMC"]

        result = ffmc(current_ffmc, temp, rh, ws, prec)
        # Use a small tolerance for floating point comparison
        assert pytest.approx(result, abs=0.1) == expected_ffmc, (
            f"Failed for row: {row}, got {result}, expected {expected_ffmc}"
        )

        current_ffmc = result
