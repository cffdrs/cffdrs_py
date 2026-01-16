import csv
import math
import pytest

def float_or_nan(value: str) -> float:
    value = value.strip()

    if not value or value.upper() == "NA":
        return math.nan

    if value.upper() == "INF":
        return math.inf

    if value.upper() == "-INF":
        return -math.inf

    return float(value)

@pytest.fixture
def load_csv():
    def _load_csv(path, schema):
        rows = []
        with open(path, newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                rows.append({
                    key: schema[key](value)
                    for key, value in row.items()
                })
        return rows
    return _load_csv


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


@pytest.fixture
def hffmc_test_data():
   """Fixture to load test data from hffmc_test_data.csv"""
   with open("cffdrs/tests/data/hffmc_test_data.csv", "r") as f:
       reader = csv.DictReader(f)
       data = []
       for row in reader:
           # convert all numeric columns to float
           row["temp"] = float(row["temp"])
           row["rh"] = float(row["rh"])
           row["ws"] = float(row["ws"])
           row["prec"] = float(row["prec"])
           row["ffmc_old"] = float(row["ffmc_old"])
           row["time_step"] = float(row["time_step"])
           row["expected_hffmc"] = float(row["expected_hffmc"])
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
    "expected_critical_surface_rate_of_spread": float_or_nan
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