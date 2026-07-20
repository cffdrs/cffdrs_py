import math
from cffdrs.constants import (
    FuelType,
    FUEL_TYPE_CODES,
    UNKNOWN_FUEL_TYPE_CODE,
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
    O1A,
    O1B,
    S1,
    S2,
    S3,
)


def _surface_fuel_consumption(
    fuel_type_code: int, ffmc: float, bui: float, pc: float, gfl: float
) -> float:
    """
    Vectorization-ready Surface Fuel Consumption Calculator.

    Same as surface_fuel_consumption(), but takes an int fuel_type_code (see
    cffdrs.constants.FUEL_TYPE_CODES) instead of a fuel type string.

    :param fuel_type_code: The Fire Behaviour Prediction fuel type code
    :param BUI: Buildup Index
    :param FFMC: Fine Fuel Moisture Code
    :param PC: Percent Conifer (%)
    :param GFL: Grass Fuel Load (kg/m^2)

    :returns: SFC Surface Fuel Consumption (kg/m^2)
    """
    sfc = -999
    if fuel_type_code == C1:
        if ffmc > 84:
            sfc = 0.75 + 0.75 * (1 - math.exp(-0.23 * (ffmc - 84))) ** 0.5
        else:
            sfc = 0.75 - 0.75 * (1 - math.exp(-0.23 * (84 - ffmc))) ** 0.5
    elif fuel_type_code in (C2, M3, M4):
        sfc = 5.0 * (1 - math.exp(-0.0115 * bui))
    elif fuel_type_code in (C3, C4):
        sfc = 5.0 * (1 - math.exp(-0.0164 * bui)) ** 2.24
    elif fuel_type_code in (C5, C6):
        sfc = 5.0 * (1 - math.exp(-0.0149 * bui)) ** 2.48
    elif fuel_type_code == C7:
        if ffmc > 70:
            sfc = 2 * (1 - math.exp(-0.104 * (ffmc - 70)))
        else:
            sfc = 0
        sfc += 1.5 * (1 - math.exp(-0.0201 * bui))
    elif fuel_type_code == D1:
        sfc = 1.5 * (1 - math.exp(-0.0183 * bui))
    elif fuel_type_code in (M1, M2):
        sfc = pc / 100 * (5.0 * (1 - math.exp(-0.0115 * bui))) + (
            (100 - pc) / 100 * (1.5 * (1 - math.exp(-0.0183 * bui)))
        )
    elif fuel_type_code in (O1A, O1B):
        sfc = gfl
    elif fuel_type_code == S1:
        sfc = 4.0 * (1 - math.exp(-0.025 * bui)) + 4.0 * (1 - math.exp(-0.034 * bui))
    elif fuel_type_code == S2:
        sfc = 10.0 * (1 - math.exp(-0.013 * bui)) + 6.0 * (1 - math.exp(-0.060 * bui))
    elif fuel_type_code == S3:
        sfc = 12.0 * (1 - math.exp(-0.0166 * bui)) + 20.0 * (1 - math.exp(-0.0210 * bui))
    if sfc <= 0:
        sfc = 0.000001
    return sfc


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
    # UNKNOWN_FUEL_TYPE_CODE never matches any branch below, which mirrors
    # the old behavior of falling through to the sfc<=0 constraint.
    return _surface_fuel_consumption(
        FUEL_TYPE_CODES.get(fuel_type, UNKNOWN_FUEL_TYPE_CODE), ffmc, bui, pc, gfl
    )
