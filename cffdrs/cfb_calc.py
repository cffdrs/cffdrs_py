import math


def critical_surface_intensity(fmc, cbh):
    """
    Critical Surface Intensity Calculator

    :param fmc: Foliar Moisture Content
    :param cbh: Crown Base Height

    :returns: CSI Critical surface intensity
    """
    # Eq. 56 (FCFDG 1992) Critical surface intensity
    CSI = 0.001 * (cbh**1.5) * (460 + 25.9 * fmc)**1.5
    return CSI


def surface_fire_rate_of_spread(csi, sfc):
    """
    Surface Fire Rate of Spread Calculator

    :param csi: Critical surface intensity
    :param sfc: Surface Fuel Consumption

    :returns: RSO Surface fire rate of spread (m/min)
    """
    try:
        if sfc == 0:
            if csi == 0:
                return math.nan  # 0 / 0 -> NA
            else:
                return math.inf  # x / 0 -> Inf
        RSO = csi / (300 * sfc)
        return RSO
    except Exception:
        return math.nan


def crown_fraction_burned(ros, rso):
    """
    Crown Fraction Burned Calculator

    :param ros: Rate of Spread
    :param rso: Surface fire rate of spread

    :returns: CFB Crown fraction burned
    """
    # Eq. 58 (FCFDG 1992) Crown fraction burned
    CFB = 1 - math.exp(-0.23 * (ros - rso)) if ros > rso else 0
    return CFB