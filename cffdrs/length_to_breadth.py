import math
from cffdrs.constants import FuelType


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
            LB = 1.1 * (wsv**0.464)
        else:
            LB = 1.0  # Eq. 80/81
    else:
        LB = 1.0 + 8.729 * (1 - math.exp(-0.030 * wsv))**2.155  # Eq. 79
    return LB