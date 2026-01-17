import math
import pytest

from cffdrs.cfb_calc import critical_surface_intensity, crown_fraction_burned, surface_fire_rate_of_spread



def test_critical_surface_intensity(cfb_calc_test_data):
    """Test the Crown Fraction Burned calculation function with test data from CSV."""

    for row in cfb_calc_test_data:
        fmc = row["fmc"]
        cbh = row["cbh"]


        # Test critical surface intensity
        expected_csi = row["expected_critical_surface_intensity"]
        csi = critical_surface_intensity(fmc, cbh)
        assert pytest.approx(expected_csi, abs=0.01) == csi, f"Failed for row: {row} - CSI: {csi}"


def test_surface_fire_rate_of_spread(cfb_calc_test_data):
    """Test the Crown Fraction Burned calculation function with test data from CSV."""

    for row in cfb_calc_test_data:
        fmc = row["fmc"]
        sfc = row["sfc"]
        cbh = row["cbh"]

        csi = critical_surface_intensity(fmc, cbh)

        expected_rso = row["expected_critical_surface_rate_of_spread"]
        rso = surface_fire_rate_of_spread(csi, sfc)
        assert pytest.approx(expected_rso, abs=0.01, nan_ok=True) == rso, f"Failed for row: {row} - RSO: {rso}"

def test_crown_fraction_burned(cfb_calc_test_data):
    """Test the Crown Fraction Burned calculation function with test data from CSV."""

    for row in cfb_calc_test_data:
        fmc = row["fmc"]
        sfc = row["sfc"]
        ros = row["ros"]
        cbh = row["cbh"]

        csi = critical_surface_intensity(fmc, cbh)
        rso = surface_fire_rate_of_spread(csi, sfc)

        expected_cfb = row["expected_cfb"]
        cfb = crown_fraction_burned(ros, rso)
        if math.isnan(expected_cfb):
            assert math.isnan(cfb) or cfb == 0, f"Failed for row: {row} - CSI: {csi}, RSO: {rso}, CFB: {cfb}"
        else:
            assert pytest.approx(expected_cfb, abs=0.01) == cfb, f"Failed for row: {row} - CSI: {csi}, RSO: {rso}, CFB: {cfb}"