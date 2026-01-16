import pytest

from cffdrs.flank_rate_of_spread import flank_rate_of_spread


csv_schema = {
    "ROS": float,
    "BROS": float,
    "LB": float,
    "expected_FlankRateOfSpread": float
}


def test_flank_rate_of_spread(load_csv):
    test_data = load_csv("cffdrs/tests/data/flank_rate_of_spread_data.csv", csv_schema)

    for row in test_data:
        ros = row["ROS"]
        bros = row["BROS"]
        lb = row["LB"]

        expected_fros = row["expected_FlankRateOfSpread"]

        calculated_fros = flank_rate_of_spread(ros, bros, lb)

        assert pytest.approx(expected_fros, 0.01) == calculated_fros, f"Failed for row: {row} - FROS: {calculated_fros}"