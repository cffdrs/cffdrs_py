import pytest
from cffdrs.crown_base_height import crown_base_height
from cffdrs.tests.conftest import float_or_nan


csv_schema = {
    "FUELTYPE": str,
    "CBH": float,
    "SD": float,
    "SH": float,
    "CrownBaseHeight": float_or_nan,
}


def test_crown_base_height(load_csv):
    test_data = load_csv("cffdrs/tests/data/CrownBaseHeight.csv", csv_schema)

    for row in test_data:
        fuel_type = row["FUELTYPE"]
        cbh = row["CBH"]
        sd = row["SD"]
        sh = row["SH"]

        expected_cbh = row["CrownBaseHeight"]

        calculated_cbh = crown_base_height(fuel_type, cbh, sd, sh)

        assert pytest.approx(expected_cbh, abs=0.01, nan_ok=True) == calculated_cbh, (
            f"Failed for row: {row} - CBH: {calculated_cbh}"
        )
