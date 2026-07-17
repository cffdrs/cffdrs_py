import math
from cffdrs.constants import FuelType, C1, O1A, O1B, S1, S2, S3, D1


def rate_of_spread_at_time_core(fuel_type_code: int, roseq, hr, cfb):
    """
    Vectorization-ready Rate of spread at time t calculation.

    Same as rate_of_spread_at_time(), but takes an int fuel_type_code (see
    cffdrs.constants.FUEL_TYPE_CODES) instead of a fuel type string.

    :param fuel_type_code: The Fire Behaviour Prediction fuel type code
    :param roseq: Equilibrium Rate of Spread (m/min)
    :param hr: Time since ignition (hours)
    :param cfb: Crown Fraction Burned

    :returns: ROSt Rate of Spread at time since ignition value
    """
    # Eq. 72 - alpha constant value, dependent on fuel type
    if fuel_type_code in (C1, O1A, O1B, S1, S2, S3, D1):
        alpha = 0.115
    else:
        if cfb < 0:
            alpha = math.nan
        else:
            alpha = 0.115 - 18.8 * (cfb**2.5) * math.exp(-8 * cfb)

    # Eq. 70 - Rate of Spread at time since ignition
    ros_t = roseq * (1 - math.exp(-alpha * hr))
    return ros_t


def rate_of_spread_at_time(fuel_type: FuelType, roseq, hr, cfb):
    """
    Rate of spread at time t calculation

    Computes the Rate of Spread prediction based on fuel type and
    FWI conditions at elapsed time since ignition. Equations are from listed
    FCFDG (1992).

    All variables names are laid out in the same manner as Forestry Canada
    Fire Danger Group (FCFDG) (1992). Development and Structure of the
    Canadian Forest Fire Behavior Prediction System." Technical Report
    ST-X-3, Forestry Canada, Ottawa, Ontario.

    :param fuel_type: The Fire Behaviour Prediction FuelType
    :param roseq: Equilibrium Rate of Spread (m/min)
    :param hr: Time since ignition (hours)
    :param cfb: Crown Fraction Burned

    :returns: ROSt Rate of Spread at time since ignition value
    """
    # Eq. 72 - alpha constant value, dependent on fuel type

    if fuel_type in ("C1", "O1A", "O1B", "S1", "S2", "S3", "D1"):
        alpha = 0.115
    else:
        # In R, negative base ** non-integer exponent → NaN
        # replicate that behavior here for test consistency
        if cfb < 0:
            alpha = math.nan
        else:
            alpha = 0.115 - 18.8 * (cfb**2.5) * math.exp(-8 * cfb)

    # Eq. 70 - Rate of Spread at time since ignition
    ros_t = roseq * (1 - math.exp(-alpha * hr))
    return ros_t
