import math
import pytest

from cffdrs.c6_calc import (
    intermediate_surface_rate_of_spread_c6,
    surface_rate_of_spread_c6,
    crown_rate_of_spread_c6,
    crown_fraction_burned_c6,
)
from cffdrs.cfb_calc import critical_surface_intensity, surface_fire_rate_of_spread


def test_intermediate_surface_rate_of_spread_c6(c6_calc_test_data):
    """Test the intermediate_surface_rate_of_spread_c6 function with test data from CSV."""

    for row in c6_calc_test_data:
        isi = row["isi"]

        expected_rsi = row["exp_intermediate_ros"]
        rsi = intermediate_surface_rate_of_spread_c6(isi)
        assert pytest.approx(expected_rsi, abs=0.01) == rsi, f"Failed for row: {row} - RSI: {rsi}"


def test_surface_rate_of_spread_c6(c6_calc_test_data):
    """Test the surface_rate_of_spread_c6 function with test data from CSV."""

    for row in c6_calc_test_data:
        isi = row["isi"]
        bui = row["bui"]

        rsi = intermediate_surface_rate_of_spread_c6(isi)

        expected_rss = row["exp_surface_ros"]
        rss = surface_rate_of_spread_c6(rsi, bui)
        assert pytest.approx(expected_rss, abs=0.01) == rss, f"Failed for row: {row} - RSS: {rss}"


def test_crown_rate_of_spread_c6(c6_calc_test_data):
    """Test the crown_rate_of_spread_c6 function with test data from CSV."""

    for row in c6_calc_test_data:
        isi = row["isi"]
        fmc = row["fmc"]

        expected_rsc = row["exp_crown_ros"]
        rsc = crown_rate_of_spread_c6(isi, fmc)
        assert pytest.approx(expected_rsc, abs=0.01) == rsc, f"Failed for row: {row} - RSC: {rsc}"


# @pytest.mark.skip(reason="Issue with NaN/0 handling in crown_fraction_burned_c6")
def test_crown_fraction_burned_c6(c6_calc_test_data):
    """Test the crown_fraction_burned_c6 function with test data from CSV."""

    for row in c6_calc_test_data:
        isi = row["isi"]
        bui = row["bui"]
        fmc = row["fmc"]
        cbh = row["cbh"]
        sfc = row["sfc"]

        rsi = intermediate_surface_rate_of_spread_c6(isi)
        rsc = crown_rate_of_spread_c6(isi, fmc)
        rss = surface_rate_of_spread_c6(rsi, bui)
        csi = critical_surface_intensity(fmc, cbh)
        rso = surface_fire_rate_of_spread(csi, sfc)

        expected_cfb = row["exp_cfb"]
        cfb = crown_fraction_burned_c6(rsc, rss, rso)
        if math.isnan(expected_cfb):
            assert math.isnan(cfb) or cfb == 0, f"Failed for row: {row} - CFB: {cfb}"
        else:
            assert pytest.approx(expected_cfb, 0.01, nan_ok=True) == cfb, (
                f"Failed for row: {row} - CFB: {cfb}"
            )
