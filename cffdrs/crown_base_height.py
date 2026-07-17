import math
from cffdrs.constants import FuelType, FUEL_TYPE_NAMES, C6

# CBH_DEFAULT[code] mirrors the CBHs dict in crown_base_height() below, indexed
# by fuel_type_code instead of fuel_type string. math.nan where undefined (NF, WA).
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


def crown_base_height_core(fuel_type_code: int, cbh, sd, sh):
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
        else:
            cbh = CBH_DEFAULT[fuel_type_code]

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
    CBHs = {
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
    invalid_cbh = cbh is None or math.isnan(cbh) or cbh <= 0 or cbh > 50

    if invalid_cbh:
        # --- R: FUELTYPE == "C6" & SD > 0 & SH > 0
        if fuel_type == "C6" and sd > 0 and sh > 0:
            cbh = -11.2 + 1.06 * sh + 0.0017 * sd
        else:
            cbh = CBHs.get(fuel_type, math.nan)

    return cbh
