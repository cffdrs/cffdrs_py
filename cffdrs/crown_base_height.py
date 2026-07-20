import math
from cffdrs.constants import FuelType, FUEL_TYPE_NAMES, FUEL_TYPE_CODES, C6

# Default crown base height by fuel type (Table 7, FCFDG 1992), indexed by
# fuel_type_code instead of fuel_type string. math.nan where undefined (NF, WA).
_CBH_DEFAULT_BY_NAME = {
    "C1": 2,
    "C2": 3,
    "C3": 8,
    "C4": 4,
    "C5": 18,
    "C6": 7,
    "C7": 10,
    "D1": 0,
    "M1": 6,
    "M2": 6,
    "M3": 6,
    "M4": 6,
    "S1": 0,
    "S2": 0,
    "S3": 0,
    "O1A": 0,
    "O1B": 0,
}
CBH_DEFAULT = tuple(_CBH_DEFAULT_BY_NAME.get(name, math.nan) for name in FUEL_TYPE_NAMES)


def _crown_base_height(fuel_type_code: int, cbh: float | None, sd: float, sh: float) -> float:
    """
    Vectorization-ready Crown Base Height function.

    Same as crown_base_height(), but takes an int fuel_type_code (see
    cffdrs.constants.FUEL_TYPE_CODES) instead of a fuel type string.

    :param fuel_type_code: The Fire Behaviour Prediction fuel type code
    :param cbh: Crown Base Height (m)
    :param sd: Stand density
    :param sh: Stand height

    :returns: Crown Base Height (m)
    """
    invalid_cbh = cbh is None or math.isnan(cbh) or cbh <= 0 or cbh > 50

    if invalid_cbh:
        if fuel_type_code == C6 and sd > 0 and sh > 0:
            cbh = -11.2 + 1.06 * sh + 0.0017 * sd
        elif 0 <= fuel_type_code < len(CBH_DEFAULT):
            cbh = CBH_DEFAULT[fuel_type_code]
        else:
            # Out-of-range/unrecognized fuel_type_code: same fallback as a
            # fuel type with no CBH default (e.g. NF, WA).
            cbh = math.nan

    return cbh


def crown_base_height(fuel_type: FuelType, cbh, sd, sh):
    """
    Crown Base Height function

    :param fuel_type: The Fire Behaviour Prediction FuelType
    :param cbh: Crown Base Height (m)
    :param sd: Stand density
    :param sh: Stand height

    :returns: Crown Base Height (m)
    """
    return _crown_base_height(FUEL_TYPE_CODES.get(fuel_type, -1), cbh, sd, sh)
