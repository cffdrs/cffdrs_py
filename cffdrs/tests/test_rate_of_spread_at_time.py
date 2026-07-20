import pytest
from cffdrs.tests.conftest import float_or_nan
from cffdrs.constants import FUEL_TYPE_CODES
from cffdrs.rate_of_spread_at_time import rate_of_spread_at_time, _rate_of_spread_at_time


csv_schema = {
    "FUELTYPE": str,
    "ROSeq": float,
    "HR": int,
    "CFB": float,
    "RateOfSpreadAtTime": float_or_nan,
}


def test_rate_of_spread_at_time(load_csv):
    test_data = load_csv("cffdrs/tests/data/RateOfSpreadAtTime.csv", csv_schema)

    for row in test_data:
        fuel_type = row["FUELTYPE"]
        roseq = row["ROSeq"]
        hr = row["HR"]
        cfb = row["CFB"]

        expected_ros_at_time = row["RateOfSpreadAtTime"]

        calculated_ros_at_time = rate_of_spread_at_time(fuel_type, roseq, hr, cfb)

        assert (
            pytest.approx(expected_ros_at_time, abs=0.01, nan_ok=True) == calculated_ros_at_time
        ), f"Failed for row: {row} - ROS at Time: {calculated_ros_at_time}"


def test_rate_of_spread_at_time_equivalence(load_csv):
    """_rate_of_spread_at_time (int fuel_type_code) must match the string version."""
    test_data = load_csv("cffdrs/tests/data/RateOfSpreadAtTime.csv", csv_schema)

    for row in test_data:
        fuel_type = row["FUELTYPE"]
        roseq = row["ROSeq"]
        hr = row["HR"]
        cfb = row["CFB"]

        expected = rate_of_spread_at_time(fuel_type, roseq, hr, cfb)
        calculated = _rate_of_spread_at_time(FUEL_TYPE_CODES[fuel_type], roseq, hr, cfb)

        assert pytest.approx(expected, abs=1e-9, nan_ok=True) == calculated, (
            f"Failed for row: {row} - ROSt core: {calculated} vs scalar: {expected}"
        )
