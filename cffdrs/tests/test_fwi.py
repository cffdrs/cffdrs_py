import pytest
from cffdrs.fwi import ffmc, dmc, dc, isi, bui, fwi


'''
Tabular test data are provided in Van Wagner and Pickett (1985): https://ostrnrcan-dostrncan.canada.ca/entities/publication/29706108-2891-4e5d-a59a-a77c96bc507c. 
This test data set is available through the R cffdrs package. The following R script produces FWI System outputs that are within one decimal place of the test data published back in 1985, 
the results of this script were saved into `fwi_test_data.csv` and a fixture was created for loading that data into tests:

library(cffdrs)
data("test_fwi")
output <- fwi(
  input = test_fwi,
  init = data.frame(ffmc = 85, dmc = 6, dc = 15, lat = 40)
)
output[, c("FFMC", "DMC", "DC", "ISI", "BUI", "FWI")] <- round(output[, c("FFMC", "DMC", "DC", "ISI", "BUI", "FWI")], digits = 1)
output[, "DSR"] <- round(output[, "DSR"], digits = 2)

# Save to CSV
write.csv(output, file = "./test_fwi_output.csv", row.names = FALSE)

'''



def test_ffmc(fwi_test_data):
    """Test ffmc function using sequential test data from fwi_test_data.csv"""
    initial_ffmc = 85.0
    current_ffmc = initial_ffmc

    for row in fwi_test_data:
        temp = row["TEMP"]
        rh = row["RH"]
        ws = row["WS"]
        prec = row["PREC"]
        expected_ffmc = row["FFMC"]

        result = ffmc(current_ffmc, temp, rh, ws, prec)
        # Use a small tolerance for floating point comparison
        assert pytest.approx(result, abs=0.1) == expected_ffmc, (
            f"Failed for row: {row}, got {result}, expected {expected_ffmc}"
        )

        current_ffmc = result


def test_dmc(fwi_test_data):
    """Test dmc function using sequential test data from fwi_test_data.csv"""
    initial_dmc = 6.0
    current_dmc = initial_dmc

    for row in fwi_test_data:
        temp = row["TEMP"]
        rh = row["RH"]
        prec = row["PREC"]
        lat = row["LAT"]
        mon = row["MON"]
        expected_dmc = row["DMC"]

        result = dmc(current_dmc, temp, rh, prec, lat, mon)
        # Use a small tolerance for floating point comparison
        assert pytest.approx(result, abs=0.1) == expected_dmc, (
            f"Failed for row: {row}, got {result}, expected {expected_dmc}"
        )

        current_dmc = result


def test_dc(fwi_test_data):
    """Test dc function using sequential test data from fwi_test_data.csv"""
    initial_dc = 15.0
    current_dc = initial_dc

    for row in fwi_test_data:
        temp = row["TEMP"]
        rh = row["RH"]
        prec = row["PREC"]
        lat = row["LAT"]
        mon = row["MON"]
        expected_dc = row["DC"]

        result = dc(current_dc, temp, rh, prec, lat, mon)
        # Use a small tolerance for floating point comparison
        assert pytest.approx(result, abs=0.1) == expected_dc, (
            f"Failed for row: {row}, got {result}, expected {expected_dc}"
        )

        current_dc = result


def test_isi(fwi_test_data):
    """Test isi function using test data from fwi_test_data.csv"""
    current_ffmc = 85.0

    for row in fwi_test_data:
        temp = row["TEMP"]
        rh = row["RH"]
        ws = row["WS"]
        prec = row["PREC"]
        expected_isi = row["ISI"]

        result_ffmc = ffmc(current_ffmc, temp, rh, ws, prec)
        result = isi(result_ffmc, ws)
        # Use a small tolerance for floating point comparison
        assert pytest.approx(result, abs=0.1) == expected_isi, (
            f"Failed for row: {row}, got {result}, expected {expected_isi}"
        )

        current_ffmc = result_ffmc


def test_bui(fwi_test_data):
    """Test bui function using test data from fwi_test_data.csv"""
    current_ffmc = 85.0
    current_dmc = 6.0
    current_dc = 15.0

    for row in fwi_test_data:
        temp = row["TEMP"]
        rh = row["RH"]
        ws = row["WS"]
        prec = row["PREC"]
        lat = row["LAT"]
        mon = row["MON"]
        expected_bui = row["BUI"]

        result_ffmc = ffmc(current_ffmc, temp, rh, ws, prec)
        result_dmc = dmc(current_dmc, temp, rh, prec, lat, mon)
        result_dc = dc(current_dc, temp, rh, prec, lat, mon)
        result = bui(result_dmc, result_dc)
        # Use a small tolerance for floating point comparison
        assert pytest.approx(result, abs=0.1) == expected_bui, (
            f"Failed for row: {row}, got {result}, expected {expected_bui}"
        )

        current_ffmc = result_ffmc
        current_dmc = result_dmc
        current_dc = result_dc


def test_fwi(fwi_test_data):
    """Test fwi function using test data from fwi_test_data.csv"""
    current_ffmc = 85.0
    current_dmc = 6.0
    current_dc = 15.0

    for row in fwi_test_data:
        temp = row["TEMP"]
        rh = row["RH"]
        ws = row["WS"]
        prec = row["PREC"]
        lat = row["LAT"]
        mon = row["MON"]
        expected_fwi = row["FWI"]

        result_ffmc = ffmc(current_ffmc, temp, rh, ws, prec)
        result_dmc = dmc(current_dmc, temp, rh, prec, lat, mon)
        result_dc = dc(current_dc, temp, rh, prec, lat, mon)
        result_isi = isi(result_ffmc, ws)
        result_bui = bui(result_dmc, result_dc)
        result = fwi(result_isi, result_bui)
        # Use a small tolerance for floating point comparison
        assert pytest.approx(result, abs=0.1) == expected_fwi, (
            f"Failed for row: {row}, got {result}, expected {expected_fwi}"
        )

        current_ffmc = result_ffmc
        current_dmc = result_dmc
        current_dc = result_dc
