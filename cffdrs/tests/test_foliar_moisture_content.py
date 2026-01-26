from cffdrs.foliar_moisture_content import foliar_moisture_content

import pytest


csv_schema = {
    "LAT": float,
    "LONG": float,
    "ELV": float,
    "DJ": int,
    "D0": int,
    "FoliarMoistureContent": float
}

def test_foliar_moisture_content(load_csv):
    test_data = load_csv("cffdrs/tests/data/FoliarMoistureContent.csv", csv_schema)

    for row in test_data:
        lat = row["LAT"]
        long = row["LONG"]
        elv = row["ELV"]
        dj = row["DJ"]
        d0 = row["D0"]

        expected_fmc = row["FoliarMoistureContent"]

        calculated_fmc = foliar_moisture_content(lat, long, elv, dj, d0)

        assert pytest.approx(expected_fmc, abs=0.01) == calculated_fmc, f"Failed for row: {row} - FMC: {calculated_fmc}"