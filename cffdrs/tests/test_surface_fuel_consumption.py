

import pytest
from cffdrs.surface_fuel_consumption import surface_fuel_consumption


csv_schema = {
    "FUELTYPE": str,
    "FFMC": float,
    "BUI": float,
    "PC": float,
    "GFL": float,
    "SurfaceFuelConsumption": float,
}

def test_surface_fuel_consumption(load_csv):
    data = load_csv("cffdrs/tests/data/SurfaceFuelConsumption.csv", csv_schema)

    for row in data:
        fuel_type = row["FUELTYPE"]
        ffmc = row["FFMC"]
        bui = row["BUI"]
        pc = row["PC"]
        gfl = row["GFL"]

        expected_consumption = row["SurfaceFuelConsumption"]

        calculated_consumption = surface_fuel_consumption(fuel_type, ffmc, bui, pc, gfl)

        assert pytest.approx(expected_consumption, abs=0.01) == calculated_consumption, f"Failed for row: {row} - Consumption: {calculated_consumption}"