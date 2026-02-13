import pytest

from cffdrs.buildup_effect import buildup_effect
from cffdrs.tests.conftest import float_or_nan


csv_schema = {
    "FUELTYPE": str,
    "BUI": float,
    "BuildupEffect": float_or_nan,
}


def test_buildup_effect(load_csv):
    """Test the buildup_effect function with test data from CSV."""
    test_data = load_csv(
        "cffdrs/tests/data/BuildupEffect.csv",
        csv_schema,
    )

    for row in test_data:
        bui = row["BUI"]
        fuel_type = row["FUELTYPE"]
        expected = row["BuildupEffect"]

        result = buildup_effect(fuel_type, bui)

        assert pytest.approx(result, abs=0.01, nan_ok=True) == expected, (
            f"Failed for row: {row} - BE: {result}"
        )
