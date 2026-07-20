import math
from cffdrs.constants import FuelType, FUEL_TYPE_CODES, O1A, O1B


def _length_to_breadth(fuel_type_code: int, wsv: float) -> float:
    """
    Vectorization-ready Length-to-Breadth ratio function.

    Same as length_to_breadth(), but takes an int fuel_type_code (see
    cffdrs.constants.FUEL_TYPE_CODES) instead of a fuel type string.

    :param fuel_type_code: The Fire Behaviour Prediction fuel type code
    :param wsv: The Wind Speed (km/h)

    :returns: Length to Breadth ratio value
    """
    if fuel_type_code in (O1A, O1B):
        if wsv >= 1.0:
            lb = 1.1 * (wsv**0.464)
        else:
            lb = 1.0
    else:
        base = 1 - math.exp(-0.030 * wsv)
        if base < 0:
            lb = math.nan
        else:
            lb = 1.0 + 8.729 * (base**2.155)

    return lb


def length_to_breadth(fuel_type: FuelType, wsv):
    """
    Length-to-Breadth ratio

    Computes the Length to Breadth ratio of an elliptically shaped
    fire. Equations are from listed FCFDG (1992) except for errata 80 from
    Wotton et. al. (2009).

    All variables names are laid out in the same manner as Forestry Canada
    Fire Danger Group (FCFDG) (1992). Development and Structure of the
    Canadian Forest Fire Behavior Prediction System." Technical Report
    ST-X-3, Forestry Canada, Ottawa, Ontario.

    Wotton, B.M., Alexander, M.E., Taylor, S.W. 2009. Updates and revisions to
    the 1992 Canadian forest fire behavior prediction system. Nat. Resour.
    Can., Can. For. Serv., Great Lakes For. Cent., Sault Ste. Marie, Ontario,
    Canada. Information Report GLC-X-10, 45p.

    :param fuel_type: The Fire Behaviour Prediction FuelType
    :param wsv: The Wind Speed (km/h)

    :returns: Length to Breadth ratio value
    """
    # -1 for an unrecognized fuel type never matches (O1A, O1B), which
    # mirrors the old fall-through to the non-grass "else" formula.
    return _length_to_breadth(FUEL_TYPE_CODES.get(fuel_type, -1), wsv)
