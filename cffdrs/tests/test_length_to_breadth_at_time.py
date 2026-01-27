import pytest
from cffdrs.length_to_breadth_at_time import length_to_breadth_at_time
from cffdrs.tests.conftest import float_or_nan


csv_schema = {
    "FUELTYPE": str,
    "LB": float,
    "HR": int,
    "CFB": float,
    "LengthToBreadthRatioAtTime": float_or_nan,
}


def test_length_to_breadth_at_time(load_csv):
    test_data = load_csv("cffdrs/tests/data/LengthToBreadthRatioAtTime.csv", csv_schema)

    for row in test_data:
        fuel_type = row["FUELTYPE"]
        lb = row["LB"]
        hr = row["HR"]
        cfb = row["CFB"]

        expected_ratio = row["LengthToBreadthRatioAtTime"]

        calculated_ratio = length_to_breadth_at_time(fuel_type, lb, hr, cfb)

        assert pytest.approx(expected_ratio, abs=0.01, nan_ok=True) == calculated_ratio, (
            f"Failed for row: {row} - LengthToBreadthRatioAtTime: {calculated_ratio}"
        )
