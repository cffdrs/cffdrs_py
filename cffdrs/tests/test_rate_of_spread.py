
import pytest
from cffdrs.rate_of_spread import rate_of_spread
from cffdrs.tests.conftest import float_or_nan


csv_schema = {
    "FUELTYPE": str,
    "ISI": float,
    "BUI": float,
    "FMC": float,
    "SFC": int,
    "PC": int,
    "PDF": int,
    "CC": int,
    "CBH": float,
    "RateOfSpread": float_or_nan
}


def test_rate_of_spread(load_csv):
    test_data = load_csv("cffdrs/tests/data/RateOfSpread.csv", csv_schema)

    for row in test_data:
        fuel_type = row["FUELTYPE"]
        isi = row["ISI"]
        bui = row["BUI"]
        fmc = row["FMC"]
        sfc = row["SFC"]
        pc = row["PC"]
        pdf = row["PDF"]
        cc = row["CC"]
        cbh = row["CBH"]

        expected_ros = row["RateOfSpread"]

        calculated_ros = rate_of_spread(
            fuel_type,
            isi,
            bui,
            fmc,
            sfc,
            pc,
            pdf,
            cc,
            cbh
        )

        assert pytest.approx(expected_ros, 0.01, nan_ok=True) == calculated_ros, f"Failed for row: {row} - ROS: {calculated_ros}"