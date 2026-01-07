import csv
import pytest


@pytest.fixture
def fwi_test_data():
    """Fixture to load test data from fwi_test_data.csv"""
    with open("cffdrs/data/fwi_test_data.csv", "r") as f:
        reader = csv.DictReader(f)
        data = list(reader)
    return data
