from dataclasses import dataclass
from typing import NamedTuple
import math
from cffdrs.constants import (
    FFMC_COEFFICIENT,
    FuelType,
    FUEL_TYPE_CODES,
    UNKNOWN_FUEL_TYPE_CODE,
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
    NF,
    WA,
)
from cffdrs.fwi import _initial_spread_index
from cffdrs.r_helpers import safe_div
from cffdrs.rate_of_spread import _rate_of_spread


class _SlopeAdjustmentOutput(NamedTuple):
    wsv: float
    raz: float


# Fuel types whose slope-equivalent ISI (isf) is derived directly from the
# basic a/b/c0 RSI formula (Eq. 41a/41b). M1-M4 and O1A/O1B are handled
# separately below.
_IS_BASIC_FUEL_TYPES = (C1, C2, C3, C4, C5, C6, C7, D1, S1, S2, S3)


def _slope_adjustment(
    fuel_type_code: int,
    ffmc: float,
    bui: float,
    ws: float,
    waz: float,
    gs: float,
    saz: float,
    fmc: float,
    sfc: float,
    pc: float,
    pdf: float,
    cc: float,
    cbh: float,
    isi: float,
) -> _SlopeAdjustmentOutput:
    """
    Vectorization-ready Slope Adjusted wind speed / direction calculation.

    Same as slope_adjustment(), but takes an int fuel_type_code (see
    cffdrs.constants.FUEL_TYPE_CODES) instead of a fuel type string. The
    M1-M4 weighted-average sub-calls, which normally recurse into
    rate_of_spread with a different, fixed fuel-type string ("C2", "D1",
    "M3", "M4"), instead call _rate_of_spread with the equivalent int
    constant - a plain function call with a compile-time-known argument,
    not a runtime-varying recursive dispatch.

    :returns: _SlopeAdjustmentOutput with wsv, raz
    """
    default_output = _SlopeAdjustmentOutput(wsv=math.nan, raz=math.nan)

    # Non-fuel or unknown fuel types -> no spread
    if fuel_type_code in (NF, WA):
        return default_output

    no_bui = -1
    # Eq. 39 (FCFDG 1992) - Calculate Spread Factor
    sf = 10 if gs >= 70 else math.exp(3.533 * (gs / 100) ** 1.2)
    # ISI with 0 wind on level grounds
    isz = _initial_spread_index(ffmc, 0)
    # Surface spread rate with 0 wind on level ground
    rsz = _rate_of_spread(fuel_type_code, isz, no_bui, fmc, sfc, pc, pdf, cc, cbh)
    # Eq. 40 (FCFDG 1992) - Surface spread rate with 0 wind upslope
    rsf = rsz * sf

    # initialize local vars
    isf = -99
    rsf_c2 = -99
    rsf_d1 = -99
    rsf_m3 = -99
    rsf_m4 = -99
    cf = -99
    isf_c2 = -99
    isf_d1 = -99
    isf_m3 = -99
    isf_m4 = -99

    pdf100 = 100

    # Eqs. 41a, 41b (Wotton 2009) - Calculate the slope equivalent ISI
    if fuel_type_code in _IS_BASIC_FUEL_TYPES:
        a_val = ROS_A[fuel_type_code]
        b_val = ROS_B[fuel_type_code]
        c0_val = ROS_C0[fuel_type_code]
        temp = 1 - (rsf / a_val) ** (1 / c0_val)
        isf = math.log(max(temp, 0.01)) / (-b_val)

    # M1/M2 weighted average
    if fuel_type_code in (M1, M2):
        rsz = _rate_of_spread(C2, isz, no_bui, fmc, sfc, pc, pdf, cc, cbh)
        rsf_c2 = rsz * sf
        rsz = _rate_of_spread(D1, isz, no_bui, fmc, sfc, pc, pdf, cc, cbh)
        rsf_d1 = rsz * sf

        isf_c2 = math.log(
            max(
                1 - (rsf_c2 / ROS_A[C2]) ** (1 / ROS_C0[C2]),
                0.01,
            )
        ) / (-ROS_B[C2])
        isf_d1 = math.log(
            max(
                1 - (rsf_d1 / ROS_A[D1]) ** (1 / ROS_C0[D1]),
                0.01,
            )
        ) / (-ROS_B[D1])
        isf = pc / 100 * isf_c2 + (1 - pc / 100) * isf_d1

    # M3 weighted average
    if fuel_type_code == M3:
        rsz = _rate_of_spread(M3, isz, no_bui, fmc, sfc, pc, pdf100, cc, cbh)
        rsf_m3 = rsz * sf
        rsz = _rate_of_spread(D1, isz, no_bui, fmc, sfc, pc, pdf100, cc, cbh)
        rsf_d1 = rsz * sf
        isf_m3 = math.log(
            max(
                1 - (rsf_m3 / ROS_A[M3]) ** (1 / ROS_C0[M3]),
                0.01,
            )
        ) / (-ROS_B[M3])
        isf_d1 = math.log(
            max(
                1 - (rsf_d1 / ROS_A[D1]) ** (1 / ROS_C0[D1]),
                0.01,
            )
        ) / (-ROS_B[D1])
        isf = pdf / 100 * isf_m3 + (1 - pdf / 100) * isf_d1

    # M4 weighted average
    if fuel_type_code == M4:
        rsz = _rate_of_spread(M4, isz, no_bui, fmc, sfc, pc, pdf100, cc, cbh)
        rsf_m4 = rsz * sf
        rsz = _rate_of_spread(D1, isz, no_bui, fmc, sfc, pc, pdf100, cc, cbh)
        rsf_d1 = rsz * sf
        isf_m4 = math.log(
            max(
                1 - (rsf_m4 / ROS_A[M4]) ** (1 / ROS_C0[M4]),
                0.01,
            )
        ) / (-ROS_B[M4])
        isf_d1 = math.log(
            max(
                1 - (rsf_d1 / ROS_A[D1]) ** (1 / ROS_C0[D1]),
                0.01,
            )
        ) / (-ROS_B[D1])
        isf = pdf / 100 * isf_m4 + (1 - pdf / 100) * isf_d1

    # Grass curing factor
    if fuel_type_code in (O1A, O1B):
        if cc < 58.8:
            cf = 0.005 * (math.exp(0.061 * cc) - 1)
        else:
            cf = 0.176 + 0.02 * (cc - 58.8)
        a_val = ROS_A[fuel_type_code]
        b_val = ROS_B[fuel_type_code]
        c0_val = ROS_C0[fuel_type_code]
        temp = 1 - safe_div(rsf, (cf * a_val)) ** (1 / c0_val)
        isf = math.log(max(temp, 0.01)) / (-b_val)

    # Only set WSV/RAZ to nan for non-spreading fuels
    if isf <= 0 or math.isnan(isf):
        return default_output

    # Eq. 46 (FCFDG 1992)
    m = FFMC_COEFFICIENT * (101 - ffmc) / (59.5 + ffmc)
    # Eq. 45 (FCFDG 1992) - FFMC function from the ISI equation
    ff = 91.9 * math.exp(-0.1386 * m) * (1 + (m**5.31) / 49300000)
    # Eqs. 44a, 44d (Wotton 2009) - Slope equivalent wind speed
    wse = 1 / 0.05039 * math.log(isf / (0.208 * ff))
    # Eqs. 44b, 44e (Wotton 2009) - Slope equivalent wind speed
    if wse > 40 and isf < (0.999 * 2.496 * ff):
        wse = 28 - (1 / 0.0818 * math.log(1 - isf / (2.496 * ff)))
    if wse > 40 and isf >= (0.999 * 2.496 * ff):
        wse = 112.45

    # Eq. 47 (FCFDG 1992) - x component
    wsx = ws * math.sin(waz) + wse * math.sin(saz)
    # Eq. 48 (FCFDG 1992) - y component
    wsy = ws * math.cos(waz) + wse * math.cos(saz)

    wsv = math.sqrt(wsx**2 + wsy**2)

    # Eq. 50 (FCFDG 1992) - the net effective wind direction (radians)
    raz = math.acos(safe_div(wsy, wsv))
    # Eq. 51 (FCFDG 1992) - convert possible negative RAZ into more understandable
    # directions
    if wsx < 0:
        raz = 2 * math.pi - raz

    return _SlopeAdjustmentOutput(wsv=wsv, raz=raz)


@dataclass
class SlopeAdjustmentOutput:
    wsv: float
    raz: float


def slope_adjustment(
    fuel_type: FuelType, ffmc, bui, ws, waz, gs, saz, fmc, sfc, pc, pdf, cc, cbh, isi
):
    """
    Slope Adjusted wind speed or slope direction of spread calculation

    Calculate the net effective windspeed (WSV), the net effective
    wind direction (RAZ) or the wind azimuth (WAZ).

    All variables names are laid out in the same manner as FCFDG (1992) and
    Wotton (2009).

    Forestry Canada Fire Danger Group (FCFDG) (1992). "Development and
    Structure of the Canadian Forest Fire Behavior Prediction System."
    Technical Report ST-X-3, Forestry Canada, Ottawa, Ontario.

    Wotton, B.M., Alexander, M.E., Taylor, S.W. 2009. Updates and revisions to
    the 1992 Canadian forest fire behavior prediction system. Nat. Resour.
    Can., Can. For. Serv., Great Lakes For. Cent., Sault Ste. Marie, Ontario,
    Canada. Information Report GLC-X-10, 45p.

    :param fuel_type: The Fire Behaviour Prediction FuelType
    :param ffmc: Fine Fuel Moisture Code
    :param bui: The Buildup Index value
    :param ws: Windspeed (km/h)
    :param waz: Wind Azimuth
    :param gs: Ground Slope (%)
    :param saz: Slope Azimuth
    :param fmc: Foliar Moisture Content
    :param sfc: Surface Fuel Consumption (kg/m^2)
    :param pc: Percent Conifer (%)
    :param pdf: Percent Dead Balsam Fir (%)
    :param cc: Constant
    :param cbh: Crown Base Height (m)
    :param isi: Initial Spread Index

    :returns: RAZ and WSV - Rate of spread azimuth (degrees) and Wind Slope speed (km/hr)
    """
    result = _slope_adjustment(
        FUEL_TYPE_CODES.get(fuel_type, UNKNOWN_FUEL_TYPE_CODE),
        ffmc,
        bui,
        ws,
        waz,
        gs,
        saz,
        fmc,
        sfc,
        pc,
        pdf,
        cc,
        cbh,
        isi,
    )
    return SlopeAdjustmentOutput(wsv=result.wsv, raz=result.raz)
