import pytest
from cffdrs.fire_intensity import fire_intensity


csv_schema = {
    "FC": int,
    "ROS": float,
    "expected_FireIntensity": int
}

def test_fire_intensity(load_csv):
    test_data = load_csv("cffdrs/tests/data/fire_intensity_data.csv", csv_schema)

    for row in test_data:
        FC = row["FC"]
        ROS = row["ROS"]

        expected_FireIntensity = row["expected_FireIntensity"]

        calculated_FireIntensity = fire_intensity(FC, ROS)

        assert pytest.approx(expected_FireIntensity) == calculated_FireIntensity, f"Failed for row: {row} - Fire Intensity: {calculated_FireIntensity}"