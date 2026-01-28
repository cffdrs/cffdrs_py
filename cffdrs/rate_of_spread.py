from dataclasses import dataclass
import math
from cffdrs.constants import FUEL_TYPE_ROS, FuelType
from cffdrs.c6_calc import (
    intermediate_surface_rate_of_spread_c6,
    crown_rate_of_spread_c6,
    surface_rate_of_spread_c6,
    crown_fraction_burned_c6,
    rate_of_spread_c6,
)
from cffdrs.cfb_calc import (
    critical_surface_intensity,
    surface_fire_rate_of_spread,
    crown_fraction_burned,
)
from cffdrs.buildup_effect import buildup_effect


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
    # HACK: C6 ROS depends on CFB so do this to not repeat calculations
    bo_bui = -1

    # Calculate RSI (set up data vectors first)
    # Eq. 26 (FCFDG 1992) - Initial Rate of Spread for Conifer and Slash types
    rsi = -1
    if fuel_type in ["C1", "C2", "C3", "C4", "C5", "C7", "D1", "S1", "S2", "S3"]:
        a_val = FUEL_TYPE_ROS[fuel_type]["a"]
        b_val = FUEL_TYPE_ROS[fuel_type]["b"]
        c0_val = FUEL_TYPE_ROS[fuel_type]["c0"]
        rsi = a_val * (1 - math.exp(-b_val * isi)) ** c0_val
    # Eq. 27 (FCFDG 1992) - Initial Rate of Spread for M1 Mixedwood type
    if fuel_type == "M1":
        rsi = pc / 100 * rate_of_spread("C2", isi, bo_bui, fmc, sfc, pc, pdf, cc, cbh) + (
            100 - pc
        ) / 100 * rate_of_spread("D1", isi, bo_bui, fmc, sfc, pc, pdf, cc, cbh)
    # Eq. 27 (FCFDG 1992) - Initial Rate of Spread for M2 Mixedwood type
    if fuel_type == "M2":
        rsi = pc / 100 * rate_of_spread("C2", isi, bo_bui, fmc, sfc, pc, pdf, cc, cbh) + 0.2 * (
            100 - pc
        ) / 100 * rate_of_spread("D1", isi, bo_bui, fmc, sfc, pc, pdf, cc, cbh)
    # Initial Rate of Spread for M3 Mixedwood
    rsi_m3 = -99
    # Eq. 30 (Wotton et. al 2009)
    if fuel_type == "M3":
        a_val = FUEL_TYPE_ROS["M3"]["a"]
        b_val = FUEL_TYPE_ROS["M3"]["b"]
        c0_val = FUEL_TYPE_ROS["M3"]["c0"]
        rsi_m3 = a_val * ((1 - math.exp(-b_val * isi)) ** c0_val)
    # Eq. 29 (Wotton et. al 2009)
    if fuel_type == "M3":
        rsi = pdf / 100 * rsi_m3 + (1 - pdf / 100) * rate_of_spread(
            "D1", isi, bo_bui, fmc, sfc, pc, pdf, cc, cbh
        )
    # Initial Rate of Spread for M4 Mixedwood
    rsi_m4 = -99
    # Eq. 30 (Wotton et. al 2009)
    if fuel_type == "M4":
        a_val = FUEL_TYPE_ROS["M4"]["a"]
        b_val = FUEL_TYPE_ROS["M4"]["b"]
        c0_val = FUEL_TYPE_ROS["M4"]["c0"]
        rsi_m4 = a_val * ((1 - math.exp(-b_val * isi)) ** c0_val)
    # Eq. 33 (Wotton et. al 2009)
    if fuel_type == "M4":
        rsi = pdf / 100 * rsi_m4 + 0.2 * (1 - pdf / 100) * rate_of_spread(
            "D1", isi, bo_bui, fmc, sfc, pc, pdf, cc, cbh
        )
    # Eq. 35b (Wotton et. al. 2009) - Calculate Curing function for grass
    cf = -99
    if fuel_type in ["O1A", "O1B"]:
        if cc < 58.8:
            cf = 0.005 * (math.exp(0.061 * cc) - 1)
        else:
            cf = 0.176 + 0.02 * (cc - 58.8)
    # Eq. 36 (FCFDG 1992) - Calculate Initial Rate of Spread for Grass
    if fuel_type in ["O1A", "O1B"]:
        a_val = FUEL_TYPE_ROS[fuel_type]["a"]
        b_val = FUEL_TYPE_ROS[fuel_type]["b"]
        c0_val = FUEL_TYPE_ROS[fuel_type]["c0"]
        rsi = a_val * ((1 - math.exp(-b_val * isi)) ** c0_val) * cf
    # used to be called like this and return ROS here
    # # Calculate the Rate of Spread (ROS)
    # ROS <- rate_of_spread(FUELTYPE, ISI, BUI, FMC, SFC, PC, PDF, CC, CBH)
    # Calculate Critical Surface Intensity
    csi = critical_surface_intensity(fmc, cbh)
    # Calculate Surface fire rate of spread (m/min)
    rso = surface_fire_rate_of_spread(csi, sfc)
    # use ifelse for C6 because ROS depends on CFB (opposite of other fuels)
    if fuel_type == "C6":
        rsi = intermediate_surface_rate_of_spread_c6(isi)
    if fuel_type == "C6":
        rsc = crown_rate_of_spread_c6(isi, fmc)
    else:
        rsc = float("nan")
    # HACK: need ROS first for non-C6
    # this is RSS for C6 and ROS otherwise
    if fuel_type == "C6":
        rss = surface_rate_of_spread_c6(rsi, bui)
    else:
        rss = buildup_effect(fuel_type, bui) * rsi
    # Calculate Crown Fraction Burned (CFB), C6 has different calculations
    if fuel_type == "C6":
        cfb = crown_fraction_burned_c6(rsc, rss, rso)
    else:
        cfb = crown_fraction_burned(rss, rso)
    if fuel_type == "C6":
        ros = rate_of_spread_c6(rsc, rss, cfb)
    else:
        ros = rss
    # add a constraint
    if ros <= 0:
        ros = 0.000001
    # CFB <- ros_and_cfb$CFB
    return RateOfSpreadOutput(ros=ros, cfb=cfb, csi=csi, rso=rso)


def rate_of_spread(fuel_type: FuelType, isi, bui, fmc, sfc, pc, pdf, cc, cbh):
    # HACK: C6 ROS depends on CFB so do this to not repeat calculations
    ros_vars = rate_of_spread_extended(fuel_type, isi, bui, fmc, sfc, pc, pdf, cc, cbh)
    return ros_vars.ros
