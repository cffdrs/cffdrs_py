import pytest
from cffdrs.distance_at_time import distance_at_time
from cffdrs.tests.conftest import float_or_nan


csv_schema = {
    "FUELTYPE": str,
    "ROSeq": float,
    "HR": int,
    "CFB": float,
    "expected_DistanceAtTime": float_or_nan
}

def test_distance_at_time(load_csv):
    test_data = load_csv("cffdrs/tests/data/distance_at_time_data.csv", csv_schema)

    for row in test_data:
        fuel_type = row["FUELTYPE"]
        roseq = row["ROSeq"]
        hr = row["HR"]
        cfb = row["CFB"]

        expected_distance = row["expected_DistanceAtTime"]

        calculated_distance = distance_at_time(fuel_type, roseq, hr, cfb)

        assert pytest.approx(expected_distance, 0.01, nan_ok=True) == calculated_distance, f"Failed for row: {row} - DistanceAtTime: {calculated_distance}"