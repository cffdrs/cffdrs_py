import math
from cffdrs.constants import FuelType


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
    # Eq. 72 (FCFDG 1992)
    # Calculate the alpha constant for the DISTt calculation
    if fuel_type in ["C1", "O1A", "O1B", "S1", "S2", "S3", "D1"]:
        alpha = 0.115
    else:
        alpha = 0.115 - 18.8 * (cfb**2.5) * math.exp(-8 * cfb)
    # Eq. 71 (FCFDG 1992) Calculate Head fire spread distance
    DISTt = roseq * (hr + math.exp(-alpha * hr) / alpha - 1 / alpha)
    return DISTt