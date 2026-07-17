import pytest
from cffdrs.rate_of_spread import rate_of_spread, rate_of_spread_core
from cffdrs.constants import FUEL_TYPE_CODES
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
    "RateOfSpread": float_or_nan,
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

        calculated_ros = rate_of_spread(fuel_type, isi, bui, fmc, sfc, pc, pdf, cc, cbh)

        assert pytest.approx(expected_ros, abs=0.01, nan_ok=True) == calculated_ros, (
            f"Failed for row: {row} - ROS: {calculated_ros}"
        )


def test_rate_of_spread_core(load_csv):
    """rate_of_spread_core (int fuel_type_code, recursion unrolled) must match rate_of_spread."""
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

        expected = rate_of_spread(fuel_type, isi, bui, fmc, sfc, pc, pdf, cc, cbh)
        calculated = rate_of_spread_core(
            FUEL_TYPE_CODES[fuel_type], isi, bui, fmc, sfc, pc, pdf, cc, cbh
        )

        assert pytest.approx(expected, abs=1e-9, nan_ok=True) == calculated, (
            f"Failed for row: {row} - ROS core: {calculated} vs scalar: {expected}"
        )
