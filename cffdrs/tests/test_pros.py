import pytest
from cffdrs.pros import pros

csv_schema = {
    "Ros": float,
    "Direction": float,
    "T1": float,
    "Long1": float,
    "Lat1": float,
    "T2": float,
    "Long2": float,
    "Lat2": float,
    "T3": float,
    "Long3": float,
    "Lat3": float,
}


def test_pros(load_csv):
    test_data = load_csv("cffdrs/tests/data/SimardRateOfSpreadPoint.csv", schema=csv_schema)

    for idx, row in enumerate(test_data, start=2):
        expected_ros = row["Ros"]
        expected_direction = row["Direction"]
        result = pros(
            row["T1"],
            row["Long1"],
            row["Lat1"],
            row["T2"],
            row["Long2"],
            row["Lat2"],
            row["T3"],
            row["Long3"],
            row["Lat3"],
        )
        ros_calculated = result["ros"]
        direction_calculated = result["direction"]

        assert pytest.approx(expected_ros, abs=0.01) == ros_calculated, (
            f"Failed ROS for row {idx}: {row}. Result: {ros_calculated}, Expected: {expected_ros}"
        )
        assert pytest.approx(expected_direction, rel=0.1) == direction_calculated, (
            f"Failed Direction for row {idx}: {row}. Result: {direction_calculated}, Expected: {expected_direction}"
        )
