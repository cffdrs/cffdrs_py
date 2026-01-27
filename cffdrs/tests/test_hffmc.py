import pytest
from cffdrs.hourly_fine_fuel_moisture_code import hourly_fine_fuel_moisture_code
from cffdrs.tests.conftest import float_or_nan

csv_schema = {
    "temp": float,
    "rh": float,
    "ws": float,
    "prec": float,
    "ffmc_old": float,
    "time.step": float,
    "HourlyFineFuelMoistureCode": float_or_nan,
}


def test_hourly_fine_fuel_moisture_code_cases(load_csv):
    hffmc_test_data = load_csv("cffdrs/tests/data/HourlyFineFuelMoistureCode.csv", csv_schema)
    for row in hffmc_test_data:
        temp = row["temp"]
        rh = row["rh"]
        ws = row["ws"]
        prec = row["prec"]
        ffmc_old = row["ffmc_old"]
        time_step = row["time.step"]
        expected_hffmc = row["HourlyFineFuelMoistureCode"]

        if rh > 100 or rh < 0 or prec < 0:
            with pytest.raises(ValueError):
                hourly_fine_fuel_moisture_code(temp, rh, ws, prec, ffmc_old, time_step)
        else:
            result = hourly_fine_fuel_moisture_code(temp, rh, ws, prec, ffmc_old, time_step)
            assert pytest.approx(expected_hffmc, abs=0.01) == result
