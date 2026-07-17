import math
from cffdrs.constants import FuelType, O1A, O1B


def length_to_breadth_core(fuel_type_code: int, wsv):
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
    # calculation is depending on if fuel type is grass (O1) or other fueltype
    if fuel_type in ["O1A", "O1B"]:
        # Correction to original Equation 80 is made here
        # Eq. 80a / 80b from Wotton 2009
        if wsv >= 1.0:
            lb = 1.1 * (wsv**0.464)
        else:
            lb = 1.0  # Eq. 80/81
    else:
        # Eq. 79
        base = 1 - math.exp(-0.030 * wsv)
        if base < 0:
            lb = math.nan
        else:
            lb = 1.0 + 8.729 * (base**2.155)

    return lb
