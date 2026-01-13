import csv
import pytest


@pytest.fixture
def fwi_test_data():
    """Fixture to load test data from fwi_test_data.csv"""
    with open("cffdrs/data/fwi_test_data.csv", "r") as f:
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
   with open("cffdrs/data/hffmc_test_data.csv", "r") as f:
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
