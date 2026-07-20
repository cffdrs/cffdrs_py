import pytest

from cffdrs.back_rate_of_spread import back_rate_of_spread, _back_rate_of_spread
from cffdrs.constants import FUEL_TYPE_CODES
from cffdrs.tests.conftest import float_or_nan

csv_schema = {
    "FUELTYPE": str,
    "FFMC": float,
    "BUI": float,
    "WSV": float,
    "FMC": float,
    "SFC": float,
    "PC": float,
    "PDF": float,
    "CC": float,
    "CBH": float,
    "BackRateOfSpread": float_or_nan,
}


def test_back_rate_of_spread(load_csv):
    """Test back_rate_of_spread function with data from BackRateOfSpread.csv"""
    test_data = load_csv("cffdrs/tests/data/BackRateOfSpread.csv", csv_schema)

    for row in test_data:
        result = back_rate_of_spread(
            fuel_type=row["FUELTYPE"],
            ffmc=row["FFMC"],
            bui=row["BUI"],
            wsv=row["WSV"],
            fmc=row["FMC"],
            sfc=row["SFC"],
            pc=row["PC"],
            pdf=row["PDF"],
            cc=row["CC"],
            cbh=row["CBH"],
        )
        assert pytest.approx(result, abs=0.01, nan_ok=True) == row["BackRateOfSpread"], (
            f"Failed for fuel_type={row['FUELTYPE']}, "
            f"ffmc={row['FFMC']}, bui={row['BUI']}, wsv={row['WSV']}, "
            f"fmc={row['FMC']}, sfc={row['SFC']}, pc={row['PC']}, "
            f"pdf={row['PDF']}, cc={row['CC']}, cbh={row['CBH']}. "
            f"Expected {row['BackRateOfSpread']}, got {result}."
        )


def test_back_rate_of_spread_equivalence(load_csv):
    """_back_rate_of_spread (int fuel_type_code) must match back_rate_of_spread (string)."""
    test_data = load_csv("cffdrs/tests/data/BackRateOfSpread.csv", csv_schema)

    for row in test_data:
        kwargs = dict(
            ffmc=row["FFMC"],
            bui=row["BUI"],
            wsv=row["WSV"],
            fmc=row["FMC"],
            sfc=row["SFC"],
            pc=row["PC"],
            pdf=row["PDF"],
            cc=row["CC"],
            cbh=row["CBH"],
        )
        expected = back_rate_of_spread(fuel_type=row["FUELTYPE"], **kwargs)
        calculated = _back_rate_of_spread(FUEL_TYPE_CODES[row["FUELTYPE"]], **kwargs)

        assert pytest.approx(expected, abs=1e-9, nan_ok=True) == calculated, (
            f"Failed for row: {row} - BROS core: {calculated} vs scalar: {expected}"
        )
