import math
from cffdrs.constants import FuelType
from cffdrs.r_helpers import safe_power


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
    if fuel_type in ["C1", "O1A", "O1B", "S1", "S2", "S3", "D1"]:
        alpha = 0.115
    else:
        alpha = 0.115 - 18.8 * (cfb**2.5) * math.exp(-8 * cfb)
    # Eq. 70 - Rate of Spread at time since ignition
    ROSt = roseq * (1 - safe_power(-alpha, hr))
    return ROSt