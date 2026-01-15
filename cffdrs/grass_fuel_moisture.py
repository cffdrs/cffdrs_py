import math
from cffdrs.constants import FFMC_COEFFICIENT


def grass_fuel_moisture(temp, rh, ws, prec, isol, gfmcold, rofl=0.3, time_step=1):
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
    MCold = FFMC_COEFFICIENT * ((101 - gfmcold) / (59.5 + gfmcold))
    # Eq. 11 - Calculate the moisture content of the layer in % after rainfall
    MCr = MCold + 100 * (prec / rofl) if prec > 0 else MCold
    # Constrain to 250
    MCr = min(MCr, 250)
    MCold = MCr
    # Eq. 2 - Calculate Fuel temperature
    Tf = temp + 35.07 * isol * math.exp(-0.06215 * ws)
    # Eq. 3 - Calculate Saturation Vapour Pressure (Baumgartner et a. 1982)
    eS_T = 6.107 * 10**(7.5 * temp / (237 + temp))
    # Eq. 3 for Fuel temperature
    eS_Tf = 6.107 * 10**(7.5 * Tf / (237 + Tf))
    # Eq. 4 - Calculate Fuel Level Relative Humidity
    RH_f = rh * (eS_T / eS_Tf)
    # Eq. 7 - Calculate Equilibrium Moisture Content for Drying phase
    EMC_D = (1.62 * RH_f**0.532 + 13.7 * math.exp((RH_f - 100) / 13.0)
             + 0.27 * (26.7 - Tf) * (1 - math.exp(-0.115 * RH_f)))
    # Eq. 7 - Calculate Equilibrium Moisture Content for Wetting phase
    EMC_W = (1.42 * RH_f**0.512 + 12.0 * math.exp((RH_f - 100) / 18.0)
             + 0.27 * (26.7 - Tf) * (1 - math.exp(-0.115 * RH_f)))
    # RH in terms of RH/100 for desorption
    Rf = RH_f / 100 if MCold > EMC_D else rh
    # RH in terms of 1-RH/100 for absorption
    Rf = (100 - RH_f) / 100 if MCold < EMC_W else Rf
    # Eq. 10 - Calculate Inverse Response time of grass (hours)
    K_GRASS = 0.389633 * math.exp(0.0365 * Tf) * (0.424 * (1 - Rf**1.7) + 0.0694 *
                math.sqrt(ws) * (1 - Rf**8))
    # Fuel is drying, calculate Moisture Content
    MC0 = (EMC_D + (MCold - EMC_D) * math.exp(-1.0 * math.log(10.0) * K_GRASS * time_step)
           if MCold > EMC_D else MCold)
    # Fuel is wetting, calculate moisture content
    MC0 = (EMC_W + (MCold - EMC_W) * math.exp(-1.0 * math.log(10.0) * K_GRASS * time_step)
           if MCold < EMC_W else MC0)
    return MC0