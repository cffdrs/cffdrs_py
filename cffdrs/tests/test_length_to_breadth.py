import pytest
from cffdrs.tests.conftest import float_or_nan
from cffdrs.length_to_breadth import length_to_breadth



csv_schema = {
    "FUELTYPE": str,
    "WSV": float,
    "LengthToBreadthRatio": float_or_nan
}

def test_length_to_breadth(load_csv):

    test_data = load_csv("cffdrs/tests/data/LengthToBreadthRatio.csv", csv_schema)

    for row in test_data:
        fuel_type = row["FUELTYPE"]
        wsv = row["WSV"]

        expected_ratio = row["LengthToBreadthRatio"]

        calculated_ratio = length_to_breadth(fuel_type, wsv)

        assert pytest.approx(expected_ratio, abs=0.01, nan_ok=True) == calculated_ratio, f"Failed for row: {row} - Length to Breadth Ratio: {calculated_ratio}"