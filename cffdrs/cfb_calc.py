import math

from cffdrs.r_helpers import safe_div


def critical_surface_intensity(fmc, cbh):
    """
    Critical Surface Intensity Calculator

    :param fmc: Foliar Moisture Content
    :param cbh: Crown Base Height

    :returns: CSI Critical surface intensity
    """
    # Eq. 56 (FCFDG 1992) Critical surface intensity
    CSI = 0.001 * (cbh**1.5) * (460 + 25.9 * fmc) ** 1.5
    return CSI


def surface_fire_rate_of_spread(csi, sfc):
    """
    Surface Fire Rate of Spread Calculator

    Follows R rules:
    - x / 0     ->  Inf or -Inf
    - 0 / 0     ->  NaN
    - NaN propagates

    :param csi: Critical surface intensity
    :param sfc: Surface Fuel Consumption

    :returns: RSO Surface fire rate of spread (m/min)
    """
    # Eq. 57 (FCFDG 1992) Surface fire rate of spread (m/min)
    rso = safe_div(csi, (300 * sfc))
    return rso


def crown_fraction_burned(ros, rso):
    """
    Crown Fraction Burned Calculator

    :param ros: Rate of Spread
    :param rso: Surface fire rate of spread

    :returns: CFB Crown fraction burned
    """
    # Eq. 58 (FCFDG 1992) Crown fraction burned
    if math.isnan(rso):
        return math.nan

    if ros > rso:
        return 1 - math.exp(-0.23 * (ros - rso))
    return 0.0
