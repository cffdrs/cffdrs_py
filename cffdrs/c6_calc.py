import math
from cffdrs.buildup_effect import buildup_effect
from cffdrs.cfb_calc import crown_fraction_burned


def intermediate_surface_rate_of_spread_c6(isi):
    """
    Intermediate Surface Rate of Spread for C6 Calculator

    :param isi: Initial Spread Index

    :returns: RSI Intermediate surface fire spread rate
    """
    # Eq. 62 (FCFDG 1992) Intermediate surface fire spread rate
    rsi = 30 * (1 - math.exp(-0.08 * isi)) ** 3.0
    return rsi


def surface_rate_of_spread_c6(rsi, bui):
    """
    Surface Rate of Spread for C6 Calculator

    :param rsi: Intermediate surface fire spread rate
    :param bui: Buildup Index

    :returns: RSS Surface fire spread rate (m/min)
    """
    # Eq. 63 (FCFDG 1992) Surface fire spread rate (m/min)
    rss = rsi * buildup_effect("C6", bui)
    return rss


def crown_rate_of_spread_c6(isi, fmc):
    """
    Crown Rate of Spread for C6 Calculator

    :param isi: Initial Spread Index
    :param fmc: Foliar Moisture Content

    :returns: RSC Crown fire spread rate (m/min)
    """
    # Average foliar moisture effect
    fme_avg = 0.778
    # Eq. 59 (FCFDG 1992) Crown flame temperature (degrees K)
    tt = 1500 - 2.75 * fmc  # this comes from the source R code but is not used
    # Eq. 60 (FCFDG 1992) Head of ignition (kJ/kg)
    H = 460 + 25.9 * fmc  # this comes from the source R code but is not used
    # Eq. 61 (FCFDG 1992) Average foliar moisture effect
    fme = ((1.5 - 0.00275 * fmc) ** 4.0) / (460 + 25.9 * fmc) * 1000
    # Eq. 64 (FCFDG 1992) Crown fire spread rate (m/min)
    rsc = 60 * (1 - math.exp(-0.0497 * isi)) * fme / fme_avg
    return rsc


def crown_fraction_burned_c6(rsc, rss, rso):
    """
    Crown Fraction Burned for C6 Calculator

    :param rsc: Crown fire spread rate
    :param rss: Surface fire spread rate
    :param rso: Surface fire rate of spread

    :returns: CFB Crown fraction burned
    """

    cond1 = rsc > rss

    if math.isnan(rso):
        cond2 = math.nan
    else:
        cond2 = rss > rso

    if cond1 is False or cond2 is False:
        return 0.0

    if math.isnan(cond2):
        return math.nan

    return crown_fraction_burned(rss, rso)


def rate_of_spread_c6(rsc, rss, cfb):
    """
    Rate of Spread for C6 Calculator

    :param rsc: Crown fire spread rate
    :param rss: Surface fire spread rate
    :param cfb: Crown fraction burned

    :returns: ROS Rate of spread (m/min)
    """
    # Eq. 65 (FCFDG 1992) Calculate Rate of spread (m/min)
    ros = rss + cfb * (rsc - rss) if rsc > rss else rss
    return ros
