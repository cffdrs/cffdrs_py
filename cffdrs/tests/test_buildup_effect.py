
import math
import pytest

from cffdrs.buildup_effect import buildup_effect
from cffdrs.tests.conftest import float_or_nan


csv_schema = {
    "fuel_type": str,
    "bui": float,
    "expected_buildup_effect": float_or_nan,
}

def test_buildup_effect(load_csv):
    """Test the buildup_effect function with test data from CSV."""
    test_data = load_csv(
        "cffdrs/tests/data/buildup_effect.csv",
        csv_schema,
    )

    for row in test_data:
        bui = row["bui"]
        fuel_type = row["fuel_type"]
        expected = row["expected_buildup_effect"]

        result = buildup_effect(fuel_type, bui)

        assert pytest.approx(result, abs=0.01, nan_ok=True) == expected, f"Failed for row: {row} - BE: {result}"