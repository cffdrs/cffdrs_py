import pytest

from cffdrs.flank_rate_of_spread import flank_rate_of_spread


csv_schema = {
    "ROS": float,
    "BROS": float,
    "LB": float,
    "FlankRateOfSpread": float
}


def test_flank_rate_of_spread(load_csv):
    test_data = load_csv("cffdrs/tests/data/FlankRateOfSpread.csv", csv_schema)

    for row in test_data:
        ros = row["ROS"]
        bros = row["BROS"]
        lb = row["LB"]

        expected_fros = row["FlankRateOfSpread"]

        calculated_fros = flank_rate_of_spread(ros, bros, lb)

        assert pytest.approx(expected_fros, abs=0.01) == calculated_fros, f"Failed for row: {row} - FROS: {calculated_fros}"