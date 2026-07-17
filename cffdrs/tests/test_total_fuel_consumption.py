import pytest
from cffdrs.total_fuel_consumption import (
    total_fuel_consumption,
    total_fuel_consumption_core,
    crown_fuel_consumption_core,
)
from cffdrs.constants import FUEL_TYPE_CODES

csv_schema = {
    "FUELTYPE": str,
    "CFL": float,
    "CFB": float,
    "SFC": int,
    "PC": int,
    "PDF": int,
    "option": str,
    "TotalFuelConsumption": float,
}


def test_total_fuel_consumption(load_csv):
    data = load_csv("cffdrs/tests/data/TotalFuelConsumption.csv", csv_schema)

    for row in data:
        fuel_type = row["FUELTYPE"]
        cfl = row["CFL"]
        cfb = row["CFB"]
        sfc = row["SFC"]
        pc = row["PC"]
        pdf = row["PDF"]
        option = row["option"]

        expected_consumption = row["TotalFuelConsumption"]

        calculated_consumption = total_fuel_consumption(fuel_type, cfl, cfb, sfc, pc, pdf, option)

        assert pytest.approx(expected_consumption, abs=0.01) == calculated_consumption, (
            f"Failed for row: {row} - Consumption: {calculated_consumption}"
        )


def test_total_fuel_consumption_core(load_csv):
    """
    crown_fuel_consumption_core/total_fuel_consumption_core (int fuel_type_code,
    no "option" string) must match total_fuel_consumption(..., option=...).
    """
    data = load_csv("cffdrs/tests/data/TotalFuelConsumption.csv", csv_schema)

    for row in data:
        fuel_type = row["FUELTYPE"]
        cfl = row["CFL"]
        cfb = row["CFB"]
        sfc = row["SFC"]
        pc = row["PC"]
        pdf = row["PDF"]
        option = row["option"]
        fuel_type_code = FUEL_TYPE_CODES[fuel_type]

        expected = total_fuel_consumption(fuel_type, cfl, cfb, sfc, pc, pdf, option)
        if option == "CFC":
            calculated = crown_fuel_consumption_core(fuel_type_code, cfl, cfb, pc, pdf)
        else:
            calculated = total_fuel_consumption_core(fuel_type_code, cfl, cfb, sfc, pc, pdf)

        assert pytest.approx(expected, abs=1e-9) == calculated, (
            f"Failed for row: {row} - core: {calculated} vs scalar: {expected}"
        )
