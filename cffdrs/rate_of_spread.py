from dataclasses import dataclass
from typing import NamedTuple
import math
from cffdrs.constants import (
    FuelType,
    FUEL_TYPE_CODES,
    ROS_A,
    ROS_B,
    ROS_C0,
    C1,
    C2,
    C3,
    C4,
    C5,
    C6,
    C7,
    D1,
    M1,
    M2,
    M3,
    M4,
    S1,
    S2,
    S3,
    O1A,
    O1B,
)
from cffdrs.c6_calc import (
    intermediate_surface_rate_of_spread_c6,
    crown_rate_of_spread_c6,
    _surface_rate_of_spread_c6,
    crown_fraction_burned_c6,
    rate_of_spread_c6,
)
from cffdrs.cfb_calc import (
    critical_surface_intensity,
    surface_fire_rate_of_spread,
    crown_fraction_burned,
)
from cffdrs.buildup_effect import _buildup_effect


class _RateOfSpreadOutput(NamedTuple):
    ros: float
    cfb: float
    csi: float
    rso: float


# The basic a*(1-exp(-b*isi))**c0 formula (Eq. 26) applies to these fuel
# types directly. M1-M4 and O1A/O1B are handled separately below; C6 is
# handled separately further down (its ROS depends on CFB).
_BASIC_ROS_FUEL_TYPES = (C1, C2, C3, C4, C5, C7, D1, S1, S2, S3)


def _floored_basic_rsi(fuel_type_code: int, isi: float) -> float:
    """
    Mirrors calling rate_of_spread(name, isi, bui=-1, ...).ros for one of the
    basic a/b/c0 fuel types (only ever C2 or D1 here). Passing bui=-1 makes
    buildup_effect short-circuit to 1.0, so that nested call reduces to the
    raw RSI formula with the same "floor to 0.000001" constraint that
    rate_of_spread_extended applies to its ros output. This is what lets the
    M1-M4 blends below be flat arithmetic instead of a recursive call into a
    different fuel type - the recursion FBP normally does here can't be
    traced/compiled by array-vectorization tools like numba or jax.
    """
    rsi = (
        ROS_A[fuel_type_code]
        * (1 - math.exp(-ROS_B[fuel_type_code] * isi)) ** ROS_C0[fuel_type_code]
    )
    return rsi if rsi > 0 else 0.000001


def _rate_of_spread_extended(
    fuel_type_code: int,
    isi: float,
    bui: float,
    fmc: float,
    sfc: float,
    pc: float,
    pdf: float,
    cc: float,
    cbh: float,
) -> _RateOfSpreadOutput:
    """
    Vectorization-ready Rate of Spread Calculation.

    Same as rate_of_spread_extended(), but takes an int fuel_type_code (see
    cffdrs.constants.FUEL_TYPE_CODES) instead of a fuel type string, and has
    the M1-M4 recursive sub-calls unrolled into flat arithmetic (see
    _floored_basic_rsi).

    :param fuel_type_code: The Fire Behaviour Prediction fuel type code
    :param isi: Initial Spread Index
    :param bui: Buildup Index
    :param fmc: Foliar Moisture Content
    :param sfc: Surface Fuel Consumption (kg/m^2)
    :param pc: Percent Conifer (%)
    :param pdf: Percent Dead Balsam Fir (%)
    :param cc: Constant
    :param cbh: Crown to base height(m)

    :returns: _RateOfSpreadOutput with ros, cfb, csi, rso
    """
    # Eq. 26 (FCFDG 1992) - Initial Rate of Spread for Conifer and Slash types
    rsi = -1
    if fuel_type_code in _BASIC_ROS_FUEL_TYPES:
        rsi = (
            ROS_A[fuel_type_code]
            * (1 - math.exp(-ROS_B[fuel_type_code] * isi)) ** ROS_C0[fuel_type_code]
        )
    # Eq. 27 (FCFDG 1992) - Initial Rate of Spread for M1 Mixedwood type
    if fuel_type_code == M1:
        rsi = pc / 100 * _floored_basic_rsi(C2, isi) + (100 - pc) / 100 * _floored_basic_rsi(
            D1, isi
        )
    # Eq. 27 (FCFDG 1992) - Initial Rate of Spread for M2 Mixedwood type
    if fuel_type_code == M2:
        rsi = pc / 100 * _floored_basic_rsi(C2, isi) + 0.2 * (100 - pc) / 100 * _floored_basic_rsi(
            D1, isi
        )
    # Initial Rate of Spread for M3 Mixedwood
    rsi_m3 = -99
    # Eq. 30 (Wotton et. al 2009)
    if fuel_type_code == M3:
        rsi_m3 = ROS_A[M3] * ((1 - math.exp(-ROS_B[M3] * isi)) ** ROS_C0[M3])
    # Eq. 29 (Wotton et. al 2009)
    if fuel_type_code == M3:
        rsi = pdf / 100 * rsi_m3 + (1 - pdf / 100) * _floored_basic_rsi(D1, isi)
    # Initial Rate of Spread for M4 Mixedwood
    rsi_m4 = -99
    # Eq. 30 (Wotton et. al 2009)
    if fuel_type_code == M4:
        rsi_m4 = ROS_A[M4] * ((1 - math.exp(-ROS_B[M4] * isi)) ** ROS_C0[M4])
    # Eq. 33 (Wotton et. al 2009)
    if fuel_type_code == M4:
        rsi = pdf / 100 * rsi_m4 + 0.2 * (1 - pdf / 100) * _floored_basic_rsi(D1, isi)
    # Eq. 35b (Wotton et. al. 2009) - Calculate Curing function for grass
    cf = -99
    if fuel_type_code in (O1A, O1B):
        if cc < 58.8:
            cf = 0.005 * (math.exp(0.061 * cc) - 1)
        else:
            cf = 0.176 + 0.02 * (cc - 58.8)
    # Eq. 36 (FCFDG 1992) - Calculate Initial Rate of Spread for Grass
    if fuel_type_code in (O1A, O1B):
        rsi = (
            ROS_A[fuel_type_code]
            * ((1 - math.exp(-ROS_B[fuel_type_code] * isi)) ** ROS_C0[fuel_type_code])
            * cf
        )
    # Calculate Critical Surface Intensity
    csi = critical_surface_intensity(fmc, cbh)
    # Calculate Surface fire rate of spread (m/min)
    rso = surface_fire_rate_of_spread(csi, sfc)
    # use ifelse for C6 because ROS depends on CFB (opposite of other fuels)
    if fuel_type_code == C6:
        rsi = intermediate_surface_rate_of_spread_c6(isi)
    if fuel_type_code == C6:
        rsc = crown_rate_of_spread_c6(isi, fmc)
    else:
        rsc = float("nan")
    # HACK: need ROS first for non-C6
    # this is RSS for C6 and ROS otherwise
    if fuel_type_code == C6:
        rss = _surface_rate_of_spread_c6(rsi, bui)
    else:
        rss = _buildup_effect(fuel_type_code, bui) * rsi
    # Calculate Crown Fraction Burned (CFB), C6 has different calculations
    if fuel_type_code == C6:
        cfb = crown_fraction_burned_c6(rsc, rss, rso)
    else:
        cfb = crown_fraction_burned(rss, rso)
    if fuel_type_code == C6:
        ros = rate_of_spread_c6(rsc, rss, cfb)
    else:
        ros = rss
    # add a constraint
    if ros <= 0:
        ros = 0.000001
    return _RateOfSpreadOutput(ros=ros, cfb=cfb, csi=csi, rso=rso)


def _rate_of_spread(
    fuel_type_code: int,
    isi: float,
    bui: float,
    fmc: float,
    sfc: float,
    pc: float,
    pdf: float,
    cc: float,
    cbh: float,
) -> float:
    """
    Vectorization-ready Rate of Spread Calculation (ros only).

    Same as rate_of_spread(), but takes an int fuel_type_code (see
    cffdrs.constants.FUEL_TYPE_CODES) instead of a fuel type string.
    """
    ros_vars = _rate_of_spread_extended(fuel_type_code, isi, bui, fmc, sfc, pc, pdf, cc, cbh)
    return ros_vars.ros


@dataclass
class RateOfSpreadOutput:
    ros: float
    cfb: float
    csi: float
    rso: float


def rate_of_spread_extended(fuel_type: FuelType, isi, bui, fmc, sfc, pc, pdf, cc, cbh):
    """
    Rate of Spread Calculation

    Computes the Rate of Spread prediction based on fuel type and
    FWI conditions. Equations are from listed FCFDG (1992) and Wotton et. al.
    (2009), and are marked as such.

    All variables names are laid out in the same manner as Forestry Canada
    Fire Danger Group (FCFDG) (1992). Development and Structure of the
    Canadian Forest Fire Behavior Prediction System." Technical Report
    ST-X-3, Forestry Canada, Ottawa, Ontario.

    Wotton, B.M., Alexander, M.E., Taylor, S.W. 2009. Updates and revisions to
    the 1992 Canadian forest fire behavior prediction system. Nat. Resour.
    Can., Can. For. Serv., Great Lakes For. Cent., Sault Ste. Marie, Ontario,
    Canada. Information Report GLC-X-10, 45p.

    :param fuel_type: The Fire Behaviour Prediction FuelType
    :param isi: Initial Spread Index
    :param bui: Buildup Index
    :param fmc: Foliar Moisture Content
    :param sfc: Surface Fuel Consumption (kg/m^2)
    :param pc: Percent Conifer (%)
    :param pdf: Percent Dead Balsam Fir (%)
    :param cc: Constant
    :param cbh: Crown to base height(m)

    :returns: dict with ROS - Rate of Spread (m/min), CFB, CSI, RSO
    """
    result = _rate_of_spread_extended(
        FUEL_TYPE_CODES.get(fuel_type, -1), isi, bui, fmc, sfc, pc, pdf, cc, cbh
    )
    return RateOfSpreadOutput(ros=result.ros, cfb=result.cfb, csi=result.csi, rso=result.rso)


def rate_of_spread(fuel_type: FuelType, isi, bui, fmc, sfc, pc, pdf, cc, cbh):
    # HACK: C6 ROS depends on CFB so do this to not repeat calculations
    return _rate_of_spread(FUEL_TYPE_CODES.get(fuel_type, -1), isi, bui, fmc, sfc, pc, pdf, cc, cbh)
