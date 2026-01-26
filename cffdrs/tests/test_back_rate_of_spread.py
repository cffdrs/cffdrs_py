import pytest

from cffdrs.back_rate_of_spread import back_rate_of_spread
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
    test_data = load_csv(
        "cffdrs/tests/data/BackRateOfSpread.csv",
        csv_schema
    )

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
            cbh=row["CBH"]
        )
        assert pytest.approx(result, abs=0.01, nan_ok=True) == row["BackRateOfSpread"], (
            f"Failed for fuel_type={row['FUELTYPE']}, "
            f"ffmc={row['FFMC']}, bui={row['BUI']}, wsv={row['WSV']}, "
            f"fmc={row['FMC']}, sfc={row['SFC']}, pc={row['PC']}, "
            f"pdf={row['PDF']}, cc={row['CC']}, cbh={row['CBH']}. "
            f"Expected {row['BackRateOfSpread']}, got {result}."
        )