import pytest
from cffdrs.lros import lros


csv_schema = {
    "Ros": float,
    "Direction": float,
    "T1": float,
    "LengthT1T2": float,
    "T2": float,
    "LengthT1T3": float,
    "T3": float,
    "LengthT2T3": float,
    "BearingT1T2": float,
    "BearingT1T3": float,
}


def test_lros(load_csv):
    test_data = load_csv("cffdrs/tests/data/SimardRateOfSpreadLine.csv", schema=csv_schema)

    for idx, row in enumerate(test_data, start=2):
        expected_ros = row["Ros"]
        expected_direction = row["Direction"]
        result = lros(
            row["T1"],
            row["LengthT1T2"],
            row["T2"],
            row["LengthT1T3"],
            row["T3"],
            row["LengthT2T3"],
            row["BearingT1T2"],
            row["BearingT1T3"],
        )
        ros_calculated = result.ros
        direction_calculated = result.direction

        assert pytest.approx(expected_ros, abs=0.01) == ros_calculated, (
            f"Failed ROS for row {idx}: {row}. Result: {ros_calculated}, Expected: {expected_ros}"
        )
        assert pytest.approx(expected_direction, abs=0.1) == direction_calculated, (
            f"Failed Direction for row {idx}: {row}. Result: {direction_calculated}, Expected: {expected_direction}"
        )
