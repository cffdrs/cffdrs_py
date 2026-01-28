import math
from cffdrs.constants import FFMC_COEFFICIENT


def grass_fuel_moisture(temp, rh, ws, prec, isol, gfmc_old, rofl=0.3, time_step=1):
    """
    Moisture content Calculation

    Calculation of moisture content for use in the GFMC calculation

    :param temp: Temperature
    :param rh: Relative Humidity
    :param ws: Wind Speed
    :param prec: Precipitation
    :param isol: Insolation
    :param gfmcold: Yesterdays Grass Foliar Moisture Content
    :param rofl: The nominal fuel load of the fine fuel layer, default is 0.3 kg/m^2
    :param time_step: Time step (hour) [default 1 hour]

    :returns: MC0
    """
    # Eq. 13 - Calculate previous moisture code
    mc_old = FFMC_COEFFICIENT * ((101 - gfmc_old) / (59.5 + gfmc_old))
    # Eq. 11 - Calculate the moisture content of the layer in % after rainfall
    mc_r = mc_old + 100 * (prec / rofl) if prec > 0 else mc_old
    # Constrain to 250
    mc_r = min(mc_r, 250)
    mc_old = mc_r
    # Eq. 2 - Calculate Fuel temperature
    tf = temp + 35.07 * isol * math.exp(-0.06215 * ws)
    # Eq. 3 - Calculate Saturation Vapour Pressure (Baumgartner et a. 1982)
    es_t = 6.107 * 10 ** (7.5 * temp / (237 + temp))
    # Eq. 3 for Fuel temperature
    es_tf = 6.107 * 10 ** (7.5 * tf / (237 + tf))
    # Eq. 4 - Calculate Fuel Level Relative Humidity
    rh_f = rh * (es_t / es_tf)
    # Eq. 7 - Calculate Equilibrium Moisture Content for Drying phase
    emc_d = (
        1.62 * rh_f**0.532
        + 13.7 * math.exp((rh_f - 100) / 13.0)
        + 0.27 * (26.7 - tf) * (1 - math.exp(-0.115 * rh_f))
    )
    # Eq. 7 - Calculate Equilibrium Moisture Content for Wetting phase
    emc_w = (
        1.42 * rh_f**0.512
        + 12.0 * math.exp((rh_f - 100) / 18.0)
        + 0.27 * (26.7 - tf) * (1 - math.exp(-0.115 * rh_f))
    )
    # RH in terms of RH/100 for desorption
    rf = rh_f / 100 if mc_old > emc_d else rh
    # RH in terms of 1-RH/100 for absorption
    rf = (100 - rh_f) / 100 if mc_old < emc_w else rf

    # Eq. 10 - Calculate Inverse Response time of grass (hours)
    try:
        exp_term = math.exp(0.0365 * tf)
    except OverflowError:  # this can create a very large number and error out
        exp_term = math.inf  # tests pass by setting to infinity
    k_grass = 0.389633 * exp_term * (0.424 * (1 - rf**1.7) + 0.0694 * math.sqrt(ws) * (1 - rf**8))

    # Fuel is drying, calculate Moisture Content
    mc0 = (
        emc_d + (mc_old - emc_d) * math.exp(-1.0 * math.log(10.0) * k_grass * time_step)
        if mc_old > emc_d
        else mc_old
    )
    # Fuel is wetting, calculate moisture content
    mc0 = (
        emc_w + (mc_old - emc_w) * math.exp(-1.0 * math.log(10.0) * k_grass * time_step)
        if mc_old < emc_w
        else mc0
    )
    return mc0
