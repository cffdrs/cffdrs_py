import math
from cffdrs.constants import FuelType, FFMC_COEFFICIENT, FUEL_TYPE_CODES, UNKNOWN_FUEL_TYPE_CODE
from cffdrs.rate_of_spread import _rate_of_spread


def _back_rate_of_spread(
    fuel_type_code: int,
    ffmc: float,
    bui: float,
    wsv: float,
    fmc: float,
    sfc: float,
    pc: float,
    pdf: float,
    cc: float,
    cbh: float,
) -> float:
    """
    Vectorization-ready Back Fire Rate of Spread Calculator.

    Same as back_rate_of_spread(), but takes an int fuel_type_code (see
    cffdrs.constants.FUEL_TYPE_CODES) instead of a fuel type string.

    :returns: BROS Back Fire Rate of Spread
    """
    # Eq. 46 (FCFDG 1992)
    # Calculate the FFMC function from the ISI equation
    m = FFMC_COEFFICIENT * (101 - ffmc) / (59.5 + ffmc)
    # Eq. 45 (FCFDG 1992)
    fF = 91.9 * math.exp(-0.1386 * m) * (1.0 + (m**5.31) / 49300000)
    # Eq. 75 (FCFDG 1992)
    # Calculate the Back fire wind function
    bf_w = math.exp(-0.05039 * wsv)
    # Calculate the ISI associated with the back fire spread rate
    # Eq. 76 (FCFDG 1992)
    bisi = 0.208 * bf_w * fF
    # Eq. 77 (FCFDG 1992)
    # Calculate final Back fire spread rate
    bros = _rate_of_spread(fuel_type_code, bisi, bui, fmc, sfc, pc, pdf, cc, cbh)
    return bros


def back_rate_of_spread(fuel_type: FuelType, ffmc, bui, wsv, fmc, sfc, pc, pdf, cc, cbh):
    """
    Back Fire Rate of Spread Calculator

    Calculate the Back Fire Spread Rate. All variables names are
    laid out in the same manner as Forestry Canada Fire Danger Group (FCFDG)
    (1992).

    :param fuel_type: The Fire Behaviour Prediction FuelType
    :param ffmc: Fine Fuel Moisture Code
    :param bui: Buildup Index
    :param wsv: Wind Speed Vector
    :param fmc: Foliar Moisture Content
    :param sfc: Surface Fuel Consumption
    :param pc: Percent Conifer
    :param pdf: Percent Dead Balsam Fir
    :param cc: Degree of Curing (just "C" in FCFDG 1992)
    :param cbh: Crown Base Height

    :returns: BROS Back Fire Rate of Spread
    """
    return _back_rate_of_spread(
        FUEL_TYPE_CODES.get(fuel_type, UNKNOWN_FUEL_TYPE_CODE),
        ffmc,
        bui,
        wsv,
        fmc,
        sfc,
        pc,
        pdf,
        cc,
        cbh,
    )
