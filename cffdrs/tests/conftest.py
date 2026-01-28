import csv
import math
import pytest
from cffdrs.models import FBPInput


def float_or_nan(value: str) -> float:
    value = value.strip()

    if not value or value.upper() == "NA":
        return math.nan

    if value.upper() == "INF":
        return math.inf

    if value.upper() == "-INF":
        return -math.inf

    return float(value)


def string_or_none(value: str) -> str:
    value = value.strip()
    if not value or value.upper() == "NA":
        return None
    return str(value)


@pytest.fixture
def load_csv():
    def _load_csv(path, schema):
        rows = []
        with open(path, newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                rows.append({key: schema[key](value) for key, value in row.items()})
        return rows

    return _load_csv


def run_csv_test(results, expected, tol=0.01, nan_ok=True):
    """
    Reusable test runner for CSV-based assertions.

    :param results: List of dicts with calculated results (keys should match expected keys case-insensitively)
    :param expected: List of dicts with expected values from CSV
    :param tol: Absolute tolerance for approx comparison
    :param nan_ok: Whether NaN values are acceptable
    """
    failures = []

    for idx, (result, exp) in enumerate(zip(results, expected), start=1):
        for key in exp.keys():
            calc = result.get(key.lower(), result.get(key))  # Handle case mismatch
            exp_val = exp[key]

            if not pytest.approx(exp_val, abs=tol, nan_ok=nan_ok) == calc:
                failures.append({"row": idx, "key": key, "expected": exp_val, "calculated": calc})
                print(f"Row {idx} failed for {key}: {calc} vs {exp_val}")

    if failures:
        print("\nSummary of all failures:")
        for f in failures:
            print(f"Row {f['row']}: {f['key']} - expected {f['expected']}, got {f['calculated']}")
        pytest.fail(f"{len(failures)} failures in CSV test. See printed details above.")


@pytest.fixture
def fwi_test_data():
    """Fixture to load test data from fwi_test_data.csv"""
    with open("cffdrs/tests/data/fwi_test_data.csv", "r") as f:
        reader = csv.DictReader(f)
        data = []
        for row in reader:
            # convert date components to int
            row["YR"] = int(row["YR"])
            row["MON"] = int(row["MON"])
            row["DAY"] = int(row["DAY"])
            # convert all other numeric columns to float
            row["LONG"] = float(row["LONG"])
            row["LAT"] = float(row["LAT"])
            row["TEMP"] = float(row["TEMP"])
            row["RH"] = float(row["RH"])
            row["WS"] = float(row["WS"])
            row["PREC"] = float(row["PREC"])
            row["FFMC"] = float(row["FFMC"])
            row["DMC"] = float(row["DMC"])
            row["DC"] = float(row["DC"])
            row["ISI"] = float(row["ISI"])
            row["BUI"] = float(row["BUI"])
            row["FWI"] = float(row["FWI"])
            row["DSR"] = float(row["DSR"])
            data.append(row)
        return data


#############
# cfb_calc test data fixture

cfb_calc_csv_schema = {
    "fuel_type": str,
    "fmc": float,
    "sfc": float,
    "ros": float,
    "cbh": float,
    "expected_cfb": float_or_nan,
    "expected_critical_surface_intensity": float_or_nan,
    "expected_critical_surface_rate_of_spread": float_or_nan,
}


@pytest.fixture
def cfb_calc_test_data(load_csv):
    """Fixture to load test data from crown_fraction_burned.csv"""
    return load_csv(
        "cffdrs/tests/data/cfb_calc_data.csv",
        cfb_calc_csv_schema,
    )


###########
# c6_calc test data fixture

c6_calc_csv_schema = {
    "fuel_type": str,
    "isi": float,
    "bui": float,
    "fmc": float,
    "sfc": float,
    "cbh": float,
    "ros": float,
    "cfb": float_or_nan,
    "rsc": float,
    "option": str,
    "exp_intermediate_ros": float,
    "exp_surface_ros": float,
    "exp_crown_ros": float,
    "exp_cfb": float_or_nan,
}


@pytest.fixture
def c6_calc_test_data(load_csv):
    """Fixture to load test data from c6_calc_test_data.csv"""
    return load_csv(
        "cffdrs/tests/data/c6_calc_test_data.csv",
        c6_calc_csv_schema,
    )


###########

# fbp test data fixture
input_csv_schema = {
    "id": int,
    "FuelType": str,
    "LAT": float,
    "LONG": float,
    "ELV": float_or_nan,
    "FFMC": float,
    "BUI": float,
    "WS": float,
    "WD": float_or_nan,
    "GS": float,
    "Dj": float,
    "D0": float_or_nan,
    "hr": float,
    "PC": float_or_nan,
    "PDF": float_or_nan,
    "GFL": float_or_nan,
    "cc": float_or_nan,
    "theta": float,
    "Accel": float,
    "Aspect": float_or_nan,
    "BUIEff": float,
    "CBH": float_or_nan,
    "CFL": float_or_nan,
    "ISI": float,
}


@pytest.fixture
def fbp_input_data(load_csv) -> list[FBPInput]:
    """Fixture to load test data from fbp_input_data.csv"""
    data = load_csv(
        "cffdrs/tests/data/test_fbp.csv",
        input_csv_schema,
    )
    inputs = []
    for row in data:
        inp = FBPInput(
            id=row["id"],
            fuel_type=row["FuelType"],
            lat=row["LAT"],
            lon=row["LONG"],
            elv=row["ELV"],
            ffmc=row["FFMC"],
            bui=row["BUI"],
            ws=row["WS"],
            wd=row["WD"],
            gs=row["GS"],
            dj=row["Dj"],
            d0=row["D0"],
            hr=row["hr"],
            pc=row["PC"],
            pdf=row["PDF"],
            gfl=row["GFL"],
            cc=row["cc"],
            theta=row["theta"],
            accel=row["Accel"],
            aspect=row["Aspect"],
            bui_eff=row["BUIEff"],
            cbh=row["CBH"],
            cfl=row["CFL"],
            isi=row["ISI"],
        )
        inputs.append(inp)
    return inputs
