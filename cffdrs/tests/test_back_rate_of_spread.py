import pytest

from cffdrs.back_rate_of_spread import back_rate_of_spread
from cffdrs.tests.conftest import float_or_nan

csv_schema = {
    "fuel_type": str,
    "ffmc": float,
    "bui": float,
    "wsv": float,
    "fmc": float,
    "sfc": float,
    "pc": float,
    "pdf": float,
    "cc": float,    
    "cbh": float,
    "expected_back_rate_of_spread": float_or_nan,
}

def test_back_rate_of_spread(load_csv):
    """Test back_rate_of_spread function with data from back_rate_of_spread.csv"""
    test_data = load_csv(
        "cffdrs/tests/data/back_rate_of_spread.csv",
        csv_schema
    )

    for row in test_data:
        result = back_rate_of_spread(
            fuel_type=row["fuel_type"],
            ffmc=row["ffmc"],
            bui=row["bui"],
            wsv=row["wsv"],
            fmc=row["fmc"],
            sfc=row["sfc"],
            pc=row["pc"],
            pdf=row["pdf"],
            cc=row["cc"],
            cbh=row["cbh"]
        )
        assert pytest.approx(result, abs=0.01, nan_ok=True) == row["expected_back_rate_of_spread"], (
            f"Failed for fuel_type={row['fuel_type']}, "
            f"ffmc={row['ffmc']}, bui={row['bui']}, wsv={row['wsv']}, "
            f"fmc={row['fmc']}, sfc={row['sfc']}, pc={row['pc']}, "
            f"pdf={row['pdf']}, cc={row['cc']}, cbh={row['cbh']}. "
            f"Expected {row['expected_back_rate_of_spread']}, got {result}."
        )