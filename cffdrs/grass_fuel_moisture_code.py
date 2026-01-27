from cffdrs.constants import FFMC_COEFFICIENT


def grass_fuel_moisture_code(mc0):
    """
    Grass Fuel Moisture Calculation

    This is the actual calculation for grass fuel moisture

    :param mc0: An output from the mcCalc functions

    :returns: gfmc
    """
    # Eq. 12 - Calculate GFMC
    GFMC0 = 59.5 * ((250 - mc0) / (FFMC_COEFFICIENT + mc0))
    return GFMC0
