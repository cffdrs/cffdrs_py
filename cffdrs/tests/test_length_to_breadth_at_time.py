import pytest
from cffdrs.length_to_breadth_at_time import (
    length_to_breadth_at_time,
    _length_to_breadth_at_time,
)
from cffdrs.constants import FUEL_TYPE_CODES
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


def test_length_to_breadth_at_time_equivalence(load_csv):
    """_length_to_breadth_at_time (int fuel_type_code) must match the string version."""
    test_data = load_csv("cffdrs/tests/data/LengthToBreadthRatioAtTime.csv", csv_schema)

    for row in test_data:
        fuel_type = row["FUELTYPE"]
        lb = row["LB"]
        hr = row["HR"]
        cfb = row["CFB"]

        expected = length_to_breadth_at_time(fuel_type, lb, hr, cfb)
        calculated = _length_to_breadth_at_time(FUEL_TYPE_CODES[fuel_type], lb, hr, cfb)

        assert pytest.approx(expected, abs=1e-9, nan_ok=True) == calculated, (
            f"Failed for row: {row} - LBt core: {calculated} vs scalar: {expected}"
        )
