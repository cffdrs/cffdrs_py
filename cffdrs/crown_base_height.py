import math
from cffdrs.constants import FuelType


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
