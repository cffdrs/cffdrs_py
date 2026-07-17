import math
from cffdrs.constants import FuelType, FUEL_TYPE_DEFAULTS, BUI_O, BUI_Q


def buildup_effect_core(fuel_type_code: int, bui):
    """
    Vectorization-ready Build Up Effect Calculator.

    Same as buildup_effect(), but takes an int fuel_type_code (see
    cffdrs.constants.FUEL_TYPE_CODES) instead of a fuel type string.

    :param fuel_type_code: The Fire Behaviour Prediction fuel type code
    :param bui: The Buildup Index value

    :returns: BE Build up effect
    """
    # Eq. 54 (FCFDG 1992) The Buildup Effect
    buio = BUI_O[fuel_type_code]
    q = BUI_Q[fuel_type_code]

    if math.isnan(buio) or math.isnan(q):
        return math.nan if bui > 0 else 1.0

    if bui > 0 and buio > 0:
        return math.exp(50 * math.log(q) * (1 / bui - 1 / buio))

    return 1.0


def buildup_effect(fuel_type: FuelType, bui):
    """
    Build Up Effect Calculator

    Computes the Buildup Effect on Fire Spread Rate. All variables
    names are laid out in the same manner as Forestry Canada Fire Danger Group
    (FCFDG)(1992).

    R behavior: if fuel unknown and BUI > 0 -> NA, else -> 1

    :param fuel_type: The Fire Behaviour Prediction FuelType
    :param bui: The Buildup Index value

    :returns: BE Build up effect
    """
    # Eq. 54 (FCFDG 1992) The Buildup Effect
    fuel = FUEL_TYPE_DEFAULTS.get(fuel_type)

    if fuel is None:
        return math.nan if bui > 0 else 1.0

    buio = fuel["BUIo"]
    q = fuel["Q"]

    if bui > 0 and buio > 0:
        return math.exp(50 * math.log(q) * (1 / bui - 1 / buio))

    return 1.0
