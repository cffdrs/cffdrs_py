import math

import pytest
from cffdrs.crown_base_height import crown_base_height
from cffdrs.tests.conftest import float_or_nan


csv_schema = {
    "fuel_type": str,
    "cbh": float,
    "sd": float,
    "sh": float,
    "expected_crown_base_height": float_or_nan
}


def test_crown_base_height(load_csv):
    test_data = load_csv("cffdrs/tests/data/crown_base_height_data.csv", csv_schema)

    for row in test_data:
        fuel_type = row["fuel_type"]
        cbh = row["cbh"]
        sd = row["sd"]
        sh = row["sh"]

        expected_cbh = row["expected_crown_base_height"]

        calculated_cbh = crown_base_height(fuel_type, cbh, sd, sh)

        assert pytest.approx(expected_cbh, 0.01, nan_ok=True) == calculated_cbh, f"Failed for row: {row} - CBH: {calculated_cbh}"