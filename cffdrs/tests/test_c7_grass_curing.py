import math

import pytest

from cffdrs.back_rate_of_spread import back_rate_of_spread
from cffdrs.fire_behaviour_prediction import fire_behaviour_prediction
from cffdrs.fwi import initial_spread_index
from cffdrs.models import FBPInput
from cffdrs.rate_of_spread import _curing_factor, rate_of_spread
from cffdrs.slope_calc import slope_adjustment

# C7 parameters
A, B, C0 = 45, 0.0305, 2.0

ROS_ARGS = (10, 50, 100, 1, 50, 35, 80, 0)  # isi, bui, fmc, sfc, pc, pdf, cc, cbh
SLOPE_ARGS = (90, 50, 15, math.pi, 25, 0, 100, 1, 50, 35, 80, 0, 0)


@pytest.mark.parametrize(
    "c7_cc, expected",
    [(0, 0), (30, 0.005 * (math.exp(0.061 * 30) - 1)), (58.8, 0.176), (80, 0.6), (100, 1)],
)
def test_curing_factor(c7_cc, expected):
    assert _curing_factor(c7_cc) == pytest.approx(expected)


def test_c7_default_curing_matches_unmodified_c7():
    expected_rsi = A * (1 - math.exp(-B * 10)) ** C0
    assert rate_of_spread("C7", 10, -1, 100, 1, 50, 35, 80, 0) == pytest.approx(expected_rsi)
    assert FBPInput(fuel_type="C7").c7_cc == 100


def test_c7_cc_does_not_use_o1_cc():
    # The O1 cc argument must not affect C7
    args = list(ROS_ARGS)
    args[6] = 0
    assert rate_of_spread("C7", *args) == pytest.approx(rate_of_spread("C7", *ROS_ARGS))


@pytest.mark.parametrize("c7_cc", [30, 58.8, 80])
def test_c7_rate_of_spread_scaled_by_curing_factor(c7_cc):
    expected = rate_of_spread("C7", *ROS_ARGS) * _curing_factor(c7_cc)
    assert rate_of_spread("C7", *ROS_ARGS, c7_cc=c7_cc) == pytest.approx(expected)


@pytest.mark.parametrize("c7_cc", [30, 58.8, 80])
def test_c7_slope_adjustment_matches_pseudo_code(c7_cc):
    ffmc, _, _, _, gs, _ = SLOPE_ARGS[:6]
    cf = _curing_factor(c7_cc)
    # Pseudo code from bcgov/wps#5844
    sf = math.exp(3.533 * (gs / 100) ** 1.2)
    isz = initial_spread_index(ffmc, 0)
    rsz = A * (1 - math.exp(-B * isz)) ** C0 * cf
    rsf = rsz * sf
    isf = math.log(max(1 - (rsf / (cf * A)) ** (1 / C0), 0.01)) / -B
    # With curing applied to both RSZ and ISF, CF cancels out, so the slope
    # equivalent ISI (and hence WSV/RAZ) is independent of curing.
    rsf_uncured = A * (1 - math.exp(-B * isz)) ** C0 * sf
    expected_isf = math.log(max(1 - (rsf_uncured / A) ** (1 / C0), 0.01)) / -B
    assert isf == pytest.approx(expected_isf)

    cured = slope_adjustment("C7", *SLOPE_ARGS, c7_cc=c7_cc)
    uncured = slope_adjustment("C7", *SLOPE_ARGS)
    assert cured.wsv == pytest.approx(uncured.wsv)
    assert cured.raz == pytest.approx(uncured.raz)


def test_c7_back_rate_of_spread_scaled_by_curing_factor():
    args = (90, 50, 20, 100, 1, 50, 35, 80, 0)
    assert back_rate_of_spread("C7", *args, c7_cc=80) == pytest.approx(
        back_rate_of_spread("C7", *args) * _curing_factor(80)
    )


def test_c7_fbp_applies_curing():
    input = {
        "fuel_type": "C7",
        "ffmc": 92,
        "bui": 80,
        "ws": 20,
        "gs": 20,
        "aspect": 180,
        "cbh": 10,
        "cfl": 0.5,
    }
    default_cured_100 = fire_behaviour_prediction(FBPInput(**input), "All")
    cured_80 = fire_behaviour_prediction(FBPInput(**input, c7_cc=80), "All")

    assert cured_80.wsv == pytest.approx(default_cured_100.wsv)
    # ROS is RSI * BE for C7 (CFB is derived from ROS, not fed back), so it scales by CF
    assert cured_80.ros == pytest.approx(default_cured_100.ros * _curing_factor(80))
    assert cured_80.bros == pytest.approx(default_cured_100.bros * _curing_factor(80))


def test_c7_cc_does_not_affect_other_fuel_types():
    for fuel_type in ("C2", "O1A", "M1"):
        assert fire_behaviour_prediction(
            FBPInput(fuel_type=fuel_type, c7_cc=40), "All"
        ) == fire_behaviour_prediction(FBPInput(fuel_type=fuel_type), "All")


@pytest.mark.parametrize("c7_cc", [-1, 101])
def test_c7_cc_out_of_range_defaults_to_100(c7_cc):
    with pytest.warns(UserWarning):
        assert FBPInput(fuel_type="C7", c7_cc=c7_cc).c7_cc == 100


def test_c7_cc_nan_defaults_to_100():
    assert FBPInput(fuel_type="C7", c7_cc=math.nan).c7_cc == 100
