import math
from cffdrs.constants import FuelType


def crown_fuel_load(fuel_type: FuelType, cfl):
    """
    Crown Fuel Load function

    :param fuel_type: The Fire Behaviour Prediction FuelType
    :param cfl: Crown Fuel Load

    :returns: Crown Fuel Load
    """
    CFLs = {
        "C1": 0.75, "C2": 0.8, "C3": 1.15, "C4": 1.2, "C5": 1.2, "C6": 1.8, "C7": 0.5,
        "D1": 0, "M1": 0.8, "M2": 0.8, "M3": 0.8, "M4": 0.8, "S1": 0, "S2": 0, "S3": 0, "O1A": 0, "O1B": 0
    }
    if cfl <= 0 or cfl > 2 or math.isnan(cfl):
        cfl = CFLs[fuel_type]
    return cfl