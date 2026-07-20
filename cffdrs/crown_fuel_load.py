import math
from cffdrs.constants import FuelType, FUEL_TYPE_NAMES, FUEL_TYPE_CODES, UNKNOWN_FUEL_TYPE_CODE

# Default crown fuel load by fuel type (Table 7, FCFDG 1992), indexed by
# fuel_type_code instead of fuel_type string. math.nan where undefined (NF, WA).
_CFL_DEFAULT_BY_NAME = {
    "C1": 0.75,
    "C2": 0.8,
    "C3": 1.15,
    "C4": 1.2,
    "C5": 1.2,
    "C6": 1.8,
    "C7": 0.5,
    "D1": 0,
    "M1": 0.8,
    "M2": 0.8,
    "M3": 0.8,
    "M4": 0.8,
    "S1": 0,
    "S2": 0,
    "S3": 0,
    "O1A": 0,
    "O1B": 0,
}
CFL_DEFAULT = tuple(_CFL_DEFAULT_BY_NAME.get(name, math.nan) for name in FUEL_TYPE_NAMES)


def _crown_fuel_load(fuel_type_code: int, cfl: float) -> float:
    """
    Vectorization-ready Crown Fuel Load function.

    Same as crown_fuel_load(), but takes an int fuel_type_code (see
    cffdrs.constants.FUEL_TYPE_CODES) instead of a fuel type string.

    :param fuel_type_code: The Fire Behaviour Prediction fuel type code
    :param cfl: Crown Fuel Load

    :returns: Crown Fuel Load
    """
    if cfl <= 0 or cfl > 2 or math.isnan(cfl):
        if 0 <= fuel_type_code < len(CFL_DEFAULT):
            cfl = CFL_DEFAULT[fuel_type_code]
        else:
            # Out-of-range/unrecognized fuel_type_code: same fallback as a
            # fuel type with no CFL default (e.g. NF, WA).
            cfl = math.nan
    return cfl


def crown_fuel_load(fuel_type: FuelType, cfl):
    """
    Crown Fuel Load function

    :param fuel_type: The Fire Behaviour Prediction FuelType
    :param cfl: Crown Fuel Load

    :returns: Crown Fuel Load
    """
    return _crown_fuel_load(FUEL_TYPE_CODES.get(fuel_type, UNKNOWN_FUEL_TYPE_CODE), cfl)
