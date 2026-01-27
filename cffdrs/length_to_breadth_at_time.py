import math
from cffdrs.constants import FuelType


def length_to_breadth_at_time(fuel_type: FuelType, lb, hr, cfb):
    """
    Length-to-Breadth ratio at time t

    Computes the Length to Breadth ratio of an elliptically shaped
    fire at elapsed time since ignition. Equations are from listed FCFDG (1992)
    and Wotton et. al. (2009), and are marked as such.

    All variables names are laid out in the same manner as Forestry Canada
    Fire Danger Group (FCFDG) (1992). Development and Structure of the
    Canadian Forest Fire Behavior Prediction System." Technical Report
    ST-X-3, Forestry Canada, Ottawa, Ontario.

    Wotton, B.M., Alexander, M.E., Taylor, S.W. 2009. Updates and revisions to
    the 1992 Canadian forest fire behavior prediction system. Nat. Resour.
    Can., Can. For. Serv., Great Lakes For. Cent., Sault Ste. Marie, Ontario,
    Canada. Information Report GLC-X-10, 45p.

    :param fuel_type: The Fire Behaviour Prediction FuelType
    :param lb: Length to Breadth ratio
    :param hr: Time since ignition (hours)
    :param cfb: Crown Fraction Burned

    :returns: Length to Breadth ratio at time since ignition
    """
    # Eq. 72 (FCFDG 1992) - alpha constant value, dependent on fuel type

    if fuel_type in ["C1", "O1A", "O1B", "S1", "S2", "S3", "D1"]:
        alpha = 0.115
    else:
        # In R, negative base ** non-integer exponent → NaN
        # Replicate that behavior for test consistency when cfb < 0
        if cfb < 0:
            alpha = math.nan
        else:
            alpha = 0.115 - 18.8 * (cfb**2.5) * math.exp(-8 * cfb)

    # Eq. 81 (Wotton et.al. 2009) - LB at time since ignition
    LBt = (lb - 1) * (1 - math.exp(-alpha * hr)) + 1
    return LBt
