import pytest
from cffdrs.distance_at_time import distance_at_time, _distance_at_time
from cffdrs.constants import FUEL_TYPE_CODES
from cffdrs.tests.conftest import float_or_nan


csv_schema = {
    "FUELTYPE": str,
    "ROSeq": float,
    "HR": int,
    "CFB": float,
    "DistanceAtTime": float_or_nan,
}


def test_distance_at_time(load_csv):
    test_data = load_csv("cffdrs/tests/data/DistanceAtTime.csv", csv_schema)

    for row in test_data:
        fuel_type = row["FUELTYPE"]
        roseq = row["ROSeq"]
        hr = row["HR"]
        cfb = row["CFB"]

        expected_distance = row["DistanceAtTime"]

        calculated_distance = distance_at_time(fuel_type, roseq, hr, cfb)

        assert pytest.approx(expected_distance, abs=0.1, nan_ok=True) == calculated_distance, (
            f"Failed for row: {row} - DistanceAtTime: {calculated_distance}"
        )


def test_distance_at_time_equivalence(load_csv):
    """_distance_at_time (int fuel_type_code) must match distance_at_time (string)."""
    test_data = load_csv("cffdrs/tests/data/DistanceAtTime.csv", csv_schema)

    for row in test_data:
        fuel_type = row["FUELTYPE"]
        roseq = row["ROSeq"]
        hr = row["HR"]
        cfb = row["CFB"]

        expected = distance_at_time(fuel_type, roseq, hr, cfb)
        calculated = _distance_at_time(FUEL_TYPE_CODES[fuel_type], roseq, hr, cfb)

        assert pytest.approx(expected, abs=1e-9, nan_ok=True) == calculated, (
            f"Failed for row: {row} - DISTt core: {calculated} vs scalar: {expected}"
        )
