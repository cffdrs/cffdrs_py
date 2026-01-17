import pytest

from cffdrs.grass_fuel_moisture_code import grass_fuel_moisture_code
from cffdrs.grass_fuel_moisture import grass_fuel_moisture

csv_schema = {
    "temp": float,
    "rh": float,
    "ws": float,
    "prec": float,
    "isol": float,
    "mon": int,
    "expected_grass_fuel_moisture": float,
    "expected_grass_fuel_moisture_code": float
}

def test_grass_fuel_moisture(load_csv):

    test_data = load_csv("cffdrs/tests/data/grass_fuel_moisture_data.csv", csv_schema)

    for row in test_data:
        temp = row["temp"]
        rh = row["rh"]
        ws = row["ws"]
        prec = row["prec"]
        isol = row["isol"]

        expected_gfmc = row["expected_grass_fuel_moisture_code"]
        expected_gfm = row["expected_grass_fuel_moisture"]

        mc = grass_fuel_moisture(
            temp=temp, rh=rh, ws=ws, prec=prec, isol=isol,
            gfmc_old=85
        )

        calculated_gfmc = grass_fuel_moisture_code(mc)

        assert pytest.approx(expected_gfm, abs=0.01) == mc, f"Failed GFM for row: {row} - GFM: {mc}"
        assert pytest.approx(expected_gfmc, abs=0.01) == calculated_gfmc, f"Failed GFMC for row: {row} - GFMC: {calculated_gfmc}"