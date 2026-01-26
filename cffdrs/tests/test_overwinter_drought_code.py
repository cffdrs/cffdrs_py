from datetime import datetime
import pytest
from cffdrs.overwinter_drought_code import overwinter_drought_code


csv_schema = {
    "DCf": float,
    "rw": float,
    "a": float,
    "b": float,
    "OverwinterDroughtCode": float,
}


def test_overwinter_drought_code(load_csv):
    test_data = load_csv("cffdrs/tests/data/OverwinterDroughtCode.csv", csv_schema)
    for row in test_data:
        dcf = row["DCf"]
        rw = row["rw"]
        a = row["a"]
        b = row["b"]

        expected_dc_start = row["OverwinterDroughtCode"]

        calculated_dc_start = overwinter_drought_code(dcf=dcf, rw=rw, a=a, b=b)

        assert pytest.approx(expected_dc_start, abs=0.1, nan_ok=True) == calculated_dc_start, (
            f"Failed for row: {row} - OverwinterDroughtCode: {calculated_dc_start}"
        )


wDC_schema = {
    "id": int,
    "lat": float,
    "long": float,
    "yr": int,
    "mon": int,
    "day": int,
    "temp": float,
    "rh": float,
    "ws": float,
    "prec": float,
    "tmax": float,
}

wDC_fs_schema = {
    "id": int,
    "lat": float,
    "long": float,
    "yr": int,
    "mon": int,
    "day": int,
    "fsdatetype": str,
}


@pytest.mark.parametrize("dcf, expected_dc", [(500, 385.07), (250, 225.94)])
def test_overwinter_drought_code_wDC_4_and_5(dcf, expected_dc, load_csv):
    """
    Specific tests to mimin the behaviour of the R code for test_overwinter_drought_code_wDC_4_and_5


    """
    test_wDC = load_csv("cffdrs/tests/data/test_wDC.csv", wDC_schema)
    test_wDC_fs = load_csv("cffdrs/tests/data/test_wDC_fs.csv", wDC_fs_schema)

    # test_wDC contains precip data for dates. Add date field for easy filtering
    input_2 = [row for row in test_wDC if row["id"] == 2]
    for row in input_2:
        row["date"] = datetime(int(row["yr"]), int(row["mon"]), int(row["day"]))

    # get fire season date from test_wDC_fs
    input_2_fs = [row for row in test_wDC_fs if row["id"] == 2]
    for row in input_2_fs:
        row["date"] = datetime(int(row["yr"]), int(row["mon"]), int(row["day"]))
    winterStartDate = input_2_fs[1]["date"]  # second row, index 1
    winterEndDate = input_2_fs[2]["date"]  # third row, index 2

    # Filter and sum precip for dates between winterStartDate and winterEndDate
    filtered = [
        row for row in input_2 if row["date"] > winterStartDate and row["date"] < winterEndDate
    ]
    curYr_prec = sum((row["prec"]) for row in filtered)

    result = overwinter_drought_code(dcf=dcf, rw=curYr_prec)

    assert pytest.approx(result, abs=0.01) == expected_dc
