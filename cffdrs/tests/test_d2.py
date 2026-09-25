import math

import pytest

from cffdrs.buildup_effect import buildup_effect
from cffdrs.constants import D2, FUEL_TYPE_CODES, FUEL_TYPE_NAMES
from cffdrs.crown_base_height import crown_base_height
from cffdrs.crown_fuel_load import crown_fuel_load
from cffdrs.distance_at_time import distance_at_time
from cffdrs.fire_behaviour_prediction import fire_behaviour_prediction
from cffdrs.length_to_breadth_at_time import length_to_breadth_at_time
from cffdrs.models import FBPInput
from cffdrs.rate_of_spread import rate_of_spread
from cffdrs.rate_of_spread_at_time import rate_of_spread_at_time
from cffdrs.slope_calc import slope_adjustment
from cffdrs.surface_fuel_consumption import surface_fuel_consumption


def test_d2_is_appended_without_changing_existing_fuel_type_codes():
    assert FUEL_TYPE_CODES["D1"] == 7
    assert FUEL_TYPE_CODES["M1"] == 8
    assert FUEL_TYPE_CODES["NF"] == 17
    assert FUEL_TYPE_CODES["WA"] == 18
    assert FUEL_TYPE_CODES["D2"] == D2 == len(FUEL_TYPE_NAMES) - 1


def test_d2_input_normalization_and_deciduous_defaults():
    assert FBPInput(fuel_type="D-2").fuel_type == "D2"
    assert crown_base_height("D2", 0, 0, 0) == 0
    assert crown_fuel_load("D2", 0) == 0


@pytest.mark.parametrize("bui", [0, 79, 79.999])
def test_d2_has_no_surface_fuel_consumption_or_buildup_effect_below_80(bui):
    assert surface_fuel_consumption("D2", ffmc=90, bui=bui, pc=50, gfl=0.35) == 0
    assert buildup_effect("D2", bui) == 0


@pytest.mark.parametrize("bui", [80, 100, 200])
def test_d2_uses_d1_surface_fuel_consumption_and_buildup_equations_at_80_or_above(bui):
    expected_sfc = 1.5 * (1 - math.exp(-0.0183 * bui))
    expected_be = math.exp(50 * math.log(0.9) * (1 / bui - 1 / 32))

    assert surface_fuel_consumption("D2", 90, bui, 50, 0.35) == pytest.approx(expected_sfc)
    assert buildup_effect("D2", bui) == pytest.approx(expected_be)


def test_d2_negative_bui_retains_internal_disable_buildup_effect_sentinel():
    assert buildup_effect("D2", -1) == 1


def test_d2_rate_of_spread_is_one_fifth_of_d1_at_or_above_threshold():
    args = {
        "isi": 10,
        "bui": 80,
        "fmc": 0,
        "sfc": 1,
        "pc": 50,
        "pdf": 35,
        "cc": 80,
        "cbh": 0,
    }

    d1_ros = rate_of_spread("D1", **args)
    d2_ros = rate_of_spread("D2", **args)

    assert d2_ros == pytest.approx(0.2 * d1_ros)


def test_d2_rate_of_spread_keeps_safety_floor_below_threshold():
    ros = rate_of_spread("D2", 10, 79.999, 0, 0, 50, 35, 80, 0)

    assert ros == 0.000001


def test_d2_slope_adjustment_accounts_for_reduced_rsi_coefficient():
    args = (90, 100, 15, math.pi, 25, 0, 0, 1, 50, 35, 80, 0, 0)

    d1 = slope_adjustment("D1", *args)
    d2 = slope_adjustment("D2", *args)

    assert d2.wsv == pytest.approx(d1.wsv)
    assert d2.raz == pytest.approx(d1.raz)


def test_d2_uses_d1_acceleration_equations():
    assert distance_at_time("D2", 5, 60, 0.5) == pytest.approx(
        distance_at_time("D1", 5, 60, 0.5)
    )
    assert rate_of_spread_at_time("D2", 5, 60, 0.5) == pytest.approx(
        rate_of_spread_at_time("D1", 5, 60, 0.5)
    )
    assert length_to_breadth_at_time("D2", 2, 60, 0.5) == pytest.approx(
        length_to_breadth_at_time("D1", 2, 60, 0.5)
    )


def test_d2_fbp_outputs_below_bui_threshold():
    result = fire_behaviour_prediction(FBPInput(fuel_type="D-2", bui=79), "All")

    assert result.sfc == 0
    assert result.be == 0
    assert result.ros == 0.000001
    assert result.tfc == 0
    assert result.hfi == 0
    assert result.fmc == 0
    assert result.cfb == 0


def test_d2_fbp_outputs_at_bui_threshold():
    input_data = FBPInput(fuel_type="D2", bui=80)
    result = fire_behaviour_prediction(input_data, "All")

    expected_sfc = surface_fuel_consumption("D2", input_data.ffmc, 80, input_data.pc, input_data.gfl)
    expected_be = buildup_effect("D2", 80)
    expected_ros = rate_of_spread(
        "D2",
        result.isi,
        80,
        result.fmc,
        result.sfc,
        input_data.pc,
        input_data.pdf,
        input_data.cc,
        0,
    )

    assert result.sfc == pytest.approx(expected_sfc)
    assert result.be == pytest.approx(expected_be)
    assert result.ros == pytest.approx(expected_ros)


def test_d2_can_disable_buildup_effect_without_disabling_spread():
    result = fire_behaviour_prediction(
        FBPInput(fuel_type="D2", bui=80, isi=10, bui_eff=0), "All"
    )

    expected_rsi = 0.2 * 30 * (1 - math.exp(-0.0232 * 10)) ** 1.6
    assert result.be == 1
    assert result.ros == pytest.approx(expected_rsi)
