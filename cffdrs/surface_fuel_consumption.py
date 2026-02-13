import math
from cffdrs.constants import FuelType


def surface_fuel_consumption(fuel_type: FuelType, ffmc, bui, pc, gfl):
    """
    Surface Fuel Consumption Calculator

    Computes the Surface Fuel Consumption by Fuel Type.
    All variables names are laid out in the same manner as FCFDG (1992) or
    Wotton et. al (2009)

    Forestry Canada Fire Danger Group (FCFDG) (1992). "Development and
    Structure of the Canadian Forest Fire Behavior Prediction System."
    Technical Report ST-X-3, Forestry Canada, Ottawa, Ontario.

    Wotton, B.M., Alexander, M.E., Taylor, S.W. 2009. Updates and revisions to
    the 1992 Canadian forest fire behavior prediction system. Nat. Resour.
    Can., Can. For. Serv., Great Lakes For. Cent., Sault Ste. Marie, Ontario,
    Canada. Information Report GLC-X-10, 45p.

    :param fuel_type: The Fire Behaviour Prediction FuelType
    :param BUI: Buildup Index
    :param FFMC: Fine Fuel Moisture Code
    :param PC: Percent Conifer (%)
    :param GFL: Grass Fuel Load (kg/m^2)

    :returns: SFC Surface Fuel Consumption (kg/m^2)
    """
    sfc = -999
    # Eqs. 9a, 9b (Wotton et. al. 2009) - Solving the lower bound of FFMC value
    # for the C1 fuel type SFC calculation
    if fuel_type == "C1":
        if ffmc > 84:
            sfc = 0.75 + 0.75 * (1 - math.exp(-0.23 * (ffmc - 84))) ** 0.5
        else:
            sfc = 0.75 - 0.75 * (1 - math.exp(-0.23 * (84 - ffmc))) ** 0.5
    # Eq. 10 (FCFDG 1992) - C2, M3, and M4 Fuel Types
    elif fuel_type == "C2" or fuel_type == "M3" or fuel_type == "M4":
        sfc = 5.0 * (1 - math.exp(-0.0115 * bui))
    # Eq. 11 (FCFDG 1992) - C3, C4 Fuel Types
    elif fuel_type == "C3" or fuel_type == "C4":
        sfc = 5.0 * (1 - math.exp(-0.0164 * bui)) ** 2.24
    # Eq. 12 (FCFDG 1992) - C5, C6 Fuel Types
    elif fuel_type == "C5" or fuel_type == "C6":
        sfc = 5.0 * (1 - math.exp(-0.0149 * bui)) ** 2.48
    # Eqs. 13, 14, 15 (FCFDG 1992) - C7 Fuel Types
    elif fuel_type == "C7":
        if ffmc > 70:
            sfc = 2 * (1 - math.exp(-0.104 * (ffmc - 70)))
        else:
            sfc = 0
        sfc += 1.5 * (1 - math.exp(-0.0201 * bui))
    # Eq. 16 (FCFDG 1992) - D1 Fuel Type
    elif fuel_type == "D1":
        sfc = 1.5 * (1 - math.exp(-0.0183 * bui))
    # Eq. 17 (FCFDG 1992) - M1 and M2 Fuel Types
    elif fuel_type == "M1" or fuel_type == "M2":
        sfc = pc / 100 * (5.0 * (1 - math.exp(-0.0115 * bui))) + (
            (100 - pc) / 100 * (1.5 * (1 - math.exp(-0.0183 * bui)))
        )
    # Eq. 18 (FCFDG 1992) - Grass Fuel Types
    elif fuel_type == "O1A" or fuel_type == "O1B":
        sfc = gfl
    # Eq. 19, 20, 25 (FCFDG 1992) - S1 Fuel Type
    elif fuel_type == "S1":
        sfc = 4.0 * (1 - math.exp(-0.025 * bui)) + 4.0 * (1 - math.exp(-0.034 * bui))
    # Eq. 21, 22, 25 (FCFDG 1992) - S2 Fuel Type
    elif fuel_type == "S2":
        sfc = 10.0 * (1 - math.exp(-0.013 * bui)) + 6.0 * (1 - math.exp(-0.060 * bui))
    # Eq. 23, 24, 25 (FCFDG 1992) - S3 Fuel Type
    elif fuel_type == "S3":
        sfc = 12.0 * (1 - math.exp(-0.0166 * bui)) + 20.0 * (1 - math.exp(-0.0210 * bui))
    # Constrain SFC value
    if sfc <= 0:
        sfc = 0.000001
    return sfc
