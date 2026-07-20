import math
from cffdrs.constants import FuelType, FUEL_TYPE_CODES, C1, O1A, O1B, S1, S2, S3, D1


def _distance_at_time(fuel_type_code: int, roseq: float, hr: float, cfb: float) -> float:
    """
    Vectorization-ready Distance at time t calculator.

    Same as distance_at_time(), but takes an int fuel_type_code (see
    cffdrs.constants.FUEL_TYPE_CODES) instead of a fuel type string.

    :param fuel_type_code: The Fire Behaviour Prediction fuel type code
    :param roseq: The predicted equilibrium rate of spread (m/min)
    :param hr: The elapsed time (min)
    :param cfb: Crown Fraction Burned

    :returns: DISTt Head fire spread distance at time t
    """
    # Eq. 72 (FCFDG 1992) - alpha constant for the DISTt calculation
    if fuel_type_code in (C1, O1A, O1B, S1, S2, S3, D1):
        alpha = 0.115
    else:
        if cfb < 0:
            alpha = math.nan
        else:
            alpha = 0.115 - 18.8 * (cfb**2.5) * math.exp(-8 * cfb)

    # Eq. 71 (FCFDG 1992) Calculate Head fire spread distance
    if math.isnan(alpha) or alpha == 0:
        dist_t = math.nan
    else:
        dist_t = roseq * (hr + math.exp(-alpha * hr) / alpha - 1 / alpha)

    return dist_t


def distance_at_time(fuel_type: FuelType, roseq, hr, cfb):
    """
    Distance at time t calculator

    Calculate the Head fire spread distance at time t. In the
    documentation this variable is just "D".

    All variables names are laid out in the same manner as Forestry Canada Fire
    Danger Group (FCFDG) (1992). Development and Structure of the  Canadian
    Forest Fire Behavior Prediction System." Technical Report ST-X-3,
    Forestry Canada, Ottawa, Ontario.

    :param fuel_type: The Fire Behaviour Prediction FuelType
    :param roseq: The predicted equilibrium rate of spread (m/min)
    :param hr: The elapsed time (min)
    :param cfb: Crown Fraction Burned

    :returns: DISTt Head fire spread distance at time t
    """
    # -1 for an unrecognized fuel type never matches the alpha=0.115 list,
    # which mirrors the old fall-through to the cfb-based "else" formula.
    return _distance_at_time(FUEL_TYPE_CODES.get(fuel_type, -1), roseq, hr, cfb)
