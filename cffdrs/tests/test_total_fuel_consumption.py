
import pytest
from cffdrs.total_fuel_consumption import total_fuel_consumption

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


        calculated_consumption = total_fuel_consumption(
            fuel_type, cfl, cfb, sfc, pc, pdf, option
        )

        assert pytest.approx(expected_consumption, 0.01) == calculated_consumption, f"Failed for row: {row} - Consumption: {calculated_consumption}"