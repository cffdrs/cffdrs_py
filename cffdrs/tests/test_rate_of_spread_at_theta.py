import pytest
from cffdrs.tests.conftest import float_or_nan
from cffdrs.rate_of_spread_at_theta import rate_of_spread_at_theta


csv_schema = {
    "ROS": float,
    "FROS": float,
    "BROS": float,
    "THETA": float,
    "RateOfSpreadAtTheta": float_or_nan,
}


def test_rate_of_spread_at_theta(load_csv):
    test_data = load_csv("cffdrs/tests/data/RateOfSpreadAtTheta.csv", csv_schema)

    for row in test_data:
        ROS = row["ROS"]
        FROS = row["FROS"]
        BROS = row["BROS"]
        THETA = row["THETA"]

        expected_ros_at_theta = row["RateOfSpreadAtTheta"]

        calculated_ros_at_theta = rate_of_spread_at_theta(ROS, FROS, BROS, THETA)

        assert (
            pytest.approx(expected_ros_at_theta, abs=0.01, nan_ok=True) == calculated_ros_at_theta
        ), f"Failed for row: {row} - ROS at Theta: {calculated_ros_at_theta}"
