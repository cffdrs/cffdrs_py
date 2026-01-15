import pytest
from cffdrs.hourly_fine_fuel_moisture_code import hourly_fine_fuel_moisture_code


@pytest.mark.parametrize(
    "temp, rh, ws, prec, ffmc_old, time_step",
    [
        (-30, -10, 0, -10, 0, 0),
        (42.9, -10, 0, -10, 0, 0),
        (25.7, -2.71, 0, -10, 0, 0),
    ],
)
def test_invalid_precip_raises_error(temp, rh, ws, prec, ffmc_old, time_step):
    with pytest.raises(ValueError):
        hourly_fine_fuel_moisture_code(temp, rh, ws, prec, ffmc_old, time_step)


@pytest.mark.parametrize(
    "temp, rh, ws, prec, ffmc_old, time_step",
    [
        (37.5, -10, 0, 55.61, 0, 0),
        (20.3, -2.71, 0, 55.61, 0, 0),
        (-19, 106.64, 0, 55.61, 0, 0),
        (53.9, 106.64, 0, 55.61, 0, 0),
    ],
)
def test_invalid_rh_raises_error(temp, rh, ws, prec, ffmc_old, time_step):
    with pytest.raises(ValueError):
        hourly_fine_fuel_moisture_code(temp, rh, ws, prec, ffmc_old, time_step)


def test_hourly_fine_fuel_moisture_code_cases(
    hffmc_test_data
):
    for row in hffmc_test_data:
        temp = row["temp"]
        rh = row["rh"]
        ws = row["ws"]
        prec = row["prec"]
        ffmc_old = row["ffmc_old"]
        time_step = row["time_step"]
        expected_hffmc = row["expected_hffmc"]
    result = hourly_fine_fuel_moisture_code(temp, rh, ws, prec, ffmc_old, time_step)
    assert pytest.approx(expected_hffmc, abs=0.01) == result
