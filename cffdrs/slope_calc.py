import math
from cffdrs.constants import FUEL_TYPE_ROS, FFMC_COEFFICIENT, FuelType
from cffdrs.fwi import initial_spread_index
from cffdrs.rate_of_spread import rate_of_spread


def slope_adjustment(fuel_type: FuelType, ffmc, bui, ws, waz, gs, saz, fmc, sfc, pc, pdf, cc, cbh, isi):
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
    NoBUI = -1
    # Eq. 39 (FCFDG 1992) - Calculate Spread Factor
    SF = 10 if gs >= 70 else math.exp(3.533 * (gs / 100)**1.2)
    # ISI with 0 wind on level grounds
    ISZ = initial_spread_index(ffmc, 0)
    # Surface spread rate with 0 wind on level ground
    RSZ = rate_of_spread(fuel_type, ISZ, NoBUI, fmc, sfc, pc, pdf, cc, cbh)
    # Eq. 40 (FCFDG 1992) - Surface spread rate with 0 wind upslope
    RSF = RSZ * SF
    # initialize some local vars
    ISF = -99
    RSF_C2 = -99
    RSF_D1 = -99
    RSF_M3 = -99
    RSF_M4 = -99
    CF = -99
    ISF_C2 = -99
    ISF_D1 = -99
    ISF_M3 = -99
    ISF_M4 = -99

    # Eqs. 41a, 41b (Wotton 2009) - Calculate the slope equivalent ISI

    is_basic = fuel_type in [
        "C1", "C2", "C3", "C4", "C5", "C6", "C7", "D1", "S1", "S2", "S3"
    ]
    if is_basic:
        a_val = FUEL_TYPE_ROS[fuel_type]["a"]
        b_val = FUEL_TYPE_ROS[fuel_type]["b"]
        c0_val = FUEL_TYPE_ROS[fuel_type]["c0"]
        temp = 1 - (RSF / a_val)**(1 / c0_val)
        if temp >= 0.01:
            ISF = math.log(temp) / (-b_val)
        else:
            ISF = math.log(0.01) / (-b_val)

    # When calculating the M1/M2 types, we are going to calculate for both C2
    # and D1 types, and combine
    # Surface spread rate with 0 wind on level ground
    if fuel_type in ["M1", "M2"]:
        RSZ = rate_of_spread("C2", ISZ, NoBUI, fmc, sfc, pc, pdf, cc, cbh)
        # Eq. 40 (FCFDG 1992) - Surface spread rate with 0 wind upslope for C2
        RSF_C2 = RSZ * SF
        RSZ = rate_of_spread("D1", ISZ, NoBUI, fmc, sfc, pc, pdf, cc, cbh)
        # Eq. 40 (FCFDG 1992) - Surface spread rate with 0 wind upslope for D1
        RSF_D1 = RSZ * SF
        RSF0 = 1 - (RSF_C2 / FUEL_TYPE_ROS["C2"]["a"])**(1 / FUEL_TYPE_ROS["C2"]["c0"])
        # Eq. 41a (Wotton 2009) - Calculate the slope equivalent ISI
        if RSF0 >= 0.01:
            ISF_C2 = math.log(1 - (RSF_C2 / FUEL_TYPE_ROS["C2"]["a"])**(1 / FUEL_TYPE_ROS["C2"]["c0"])) / (-FUEL_TYPE_ROS["C2"]["b"])
        else:
            # Eq. 41b (Wotton 2009) - Calculate the slope equivalent ISI
            ISF_C2 = math.log(0.01) / (-FUEL_TYPE_ROS["C2"]["b"])
        RSF0 = 1 - (RSF_D1 / FUEL_TYPE_ROS["D1"]["a"])**(1 / FUEL_TYPE_ROS["D1"]["c0"])
        # Eq. 41a (Wotton 2009) - Calculate the slope equivalent ISI
        if RSF0 >= 0.01:
            ISF_D1 = math.log(1 - (RSF_D1 / FUEL_TYPE_ROS["D1"]["a"])**(1 / FUEL_TYPE_ROS["D1"]["c0"])) / (-FUEL_TYPE_ROS["D1"]["b"])
        else:
            # Eq. 41b (Wotton 2009) - Calculate the slope equivalent ISI
            ISF_D1 = math.log(0.01) / (-FUEL_TYPE_ROS["D1"]["b"])
        # Eq. 42a (Wotton 2009) - Calculate weighted average for the M1/M2 types
        ISF = pc / 100 * ISF_C2 + (1 - pc / 100) * ISF_D1

    # Set % Dead Balsam Fir to 100%
    PDF100 = 100
    # Surface spread rate with 0 wind on level ground
    if fuel_type == "M3":
        RSZ = rate_of_spread("M3", ISZ, NoBUI, fmc, sfc, pc, PDF100, cc, cbh)
        # Eq. 40 (FCFDG 1992) - Surface spread rate with 0 wind upslope for M3
        RSF_M3 = RSZ * SF
        # Surface spread rate with 0 wind on level ground, using D1
        RSZ = rate_of_spread("D1", ISZ, NoBUI, fmc, sfc, pc, PDF100, cc, cbh)
        # Eq. 40 (FCFDG 1992) - Surface spread rate with 0 wind upslope for M3
        RSF_D1 = RSZ * SF
        RSF0 = 1 - (RSF_M3 / FUEL_TYPE_ROS["M3"]["a"])**(1 / FUEL_TYPE_ROS["M3"]["c0"])
        # Eq. 41a (Wotton 2009) - Calculate the slope equivalent ISI
        if RSF0 >= 0.01:
            ISF_M3 = math.log(1 - (RSF_M3 / FUEL_TYPE_ROS["M3"]["a"])**(1 / FUEL_TYPE_ROS["M3"]["c0"])) / (-FUEL_TYPE_ROS["M3"]["b"])
        else:
            # Eq. 41b (Wotton 2009) - Calculate the slope equivalent ISI
            ISF_M3 = math.log(0.01) / (-FUEL_TYPE_ROS["M3"]["b"])
        # Eq. 40 (FCFDG 1992) - Surface spread rate with 0 wind upslope for D1
        RSF0 = 1 - (RSF_D1 / FUEL_TYPE_ROS["D1"]["a"])**(1 / FUEL_TYPE_ROS["D1"]["c0"])
        # Eq. 41a (Wotton 2009) - Calculate the slope equivalent ISI
        if RSF0 >= 0.01:
            ISF_D1 = math.log(1 - (RSF_D1 / FUEL_TYPE_ROS["D1"]["a"])**(1 / FUEL_TYPE_ROS["D1"]["c0"])) / (-FUEL_TYPE_ROS["D1"]["b"])
        else:
            # Eq. 41b (Wotton 2009) - Calculate the slope equivalent ISI
            ISF_D1 = math.log(0.01) / (-FUEL_TYPE_ROS["D1"]["b"])
        # Eq. 42b (Wotton 2009) - Calculate weighted average for the M3 type
        ISF = pdf / 100 * ISF_M3 + (1 - pdf / 100) * ISF_D1
    # Surface spread rate with 0 wind on level ground, using M4
    if fuel_type == "M4":
        RSZ = rate_of_spread("M4", ISZ, NoBUI, fmc, sfc, pc, PDF100, cc, cbh)
        # Eq. 40 (FCFDG 1992) - Surface spread rate with 0 wind upslope for M4
        RSF_M4 = RSZ * SF
        # Surface spread rate with 0 wind on level ground, using M4
        RSZ = rate_of_spread("D1", ISZ, NoBUI, fmc, sfc, pc, PDF100, cc, cbh)
        # Eq. 40 (FCFDG 1992) - Surface spread rate with 0 wind upslope for D1
        RSF_D1 = RSZ * SF
        # Eq. 40 (FCFDG 1992) - Surface spread rate with 0 wind upslope for D1
        RSF0 = 1 - (RSF_M4 / FUEL_TYPE_ROS["M4"]["a"])**(1 / FUEL_TYPE_ROS["M4"]["c0"])
        # Eq. 41a (Wotton 2009) - Calculate the slope equivalent ISI
        if RSF0 >= 0.01:
            ISF_M4 = math.log(1 - (RSF_M4 / FUEL_TYPE_ROS["M4"]["a"])**(1 / FUEL_TYPE_ROS["M4"]["c0"])) / (-FUEL_TYPE_ROS["M4"]["b"])
        else:
            # Eq. 41b (Wotton 2009) - Calculate the slope equivalent ISI
            ISF_M4 = math.log(0.01) / (-FUEL_TYPE_ROS["M4"]["b"])
        # Eq. 40 (FCFDG 1992) - Surface spread rate with 0 wind upslope for D1
        RSF0 = 1 - (RSF_D1 / FUEL_TYPE_ROS["D1"]["a"])**(1 / FUEL_TYPE_ROS["D1"]["c0"])
        # Eq. 41a (Wotton 2009) - Calculate the slope equivalent ISI (D1)
        if RSF0 >= 0.01:
            ISF_D1 = math.log(1 - (RSF_D1 / FUEL_TYPE_ROS["D1"]["a"])**(1 / FUEL_TYPE_ROS["D1"]["c0"])) / (-FUEL_TYPE_ROS["D1"]["b"])
        else:
            # Eq. 41b (Wotton 2009) - Calculate the slope equivalent ISI (D1)
            ISF_D1 = math.log(0.01) / (-FUEL_TYPE_ROS["D1"]["b"])
        # Eq. 42c (Wotton 2009) - Calculate weighted average for the M4 type
        ISF = pdf / 100 * ISF_M4 + (1 - pdf / 100) * ISF_D1
    # Eqs. 35a, 35b (Wotton 2009) - Curing Factor pivoting around % 58.8
    if fuel_type in ["O1A", "O1B"]:
        if cc < 58.8:
            CF = 0.005 * (math.exp(0.061 * cc) - 1)
        else:
            CF = 0.176 + 0.02 * (cc - 58.8)
        # Eqs. 43a, 43b (Wotton 2009) - slope equivilent ISI for Grass
        a_val = FUEL_TYPE_ROS[fuel_type]["a"]
        b_val = FUEL_TYPE_ROS[fuel_type]["b"]
        c0_val = FUEL_TYPE_ROS[fuel_type]["c0"]
        temp = 1 - (RSF / (CF * a_val))**(1 / c0_val)
        if temp >= 0.01:
            ISF = math.log(temp) / (-b_val)
        else:
            ISF = math.log(0.01) / (-b_val)
    # Initialize RAZ and WSV
    RAZ = -99
    WSV = -99
    # Eq. 46 (FCFDG 1992)
    m = FFMC_COEFFICIENT * (101 - ffmc) / (59.5 + ffmc)
    # Eq. 45 (FCFDG 1992) - FFMC function from the ISI equation
    fF = 91.9 * math.exp(-0.1386 * m) * (1 + (m**5.31) / 49300000)
    # Eqs. 44a, 44d (Wotton 2009) - Slope equivalent wind speed
    WSE = 1 / 0.05039 * math.log(ISF / (0.208 * fF))
    # Eqs. 44b, 44e (Wotton 2009) - Slope equivalent wind speed
    if WSE > 40 and ISF < (0.999 * 2.496 * fF):
        WSE = 28 - (1 / 0.0818 * math.log(1 - ISF / (2.496 * fF)))
    # Eqs. 44c (Wotton 2009) - Slope equivalent wind speed
    if WSE > 40 and ISF >= (0.999 * 2.496 * fF):
        WSE = 112.45
    # Eq. 47 (FCFDG 1992) - resultant vector magnitude in the x-direction
    WSX = ws * math.sin(math.radians(waz)) + WSE * math.sin(math.radians(saz))
    # Eq. 48 (FCFDG 1992) - resultant vector magnitude in the y-direction
    WSY = ws * math.cos(math.radians(waz)) + WSE * math.cos(math.radians(saz))
    # Eq. 49 (FCFDG 1992) - the net effective wind speed
    WSV = math.sqrt(WSX * WSX + WSY * WSY)
    if fuel_type in ["NF", "WA"]:
        WSV = float('nan')
    # Eq. 50 (FCFDG 1992) - the net effective wind direction (radians)
    RAZ = math.acos(WSY / WSV)
    # Eq. 51 (FCFDG 1992) - convert possible negative RAZ into more understandable
    # directions
    if WSX < 0:
        RAZ = 2 * math.pi - RAZ
    if fuel_type in ["NF", "WA"]:
        RAZ = float('nan')
    # Convert RAZ to degrees
    RAZ = math.degrees(RAZ)
    return {"WSV": WSV, "RAZ": RAZ}