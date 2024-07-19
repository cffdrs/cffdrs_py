from math import exp
from cffdrs.fwi import FFMC_COEFFICIENT


def hourly_fine_fuel_moisture_code(temp, rh, ws, prec, fo=85, t0=1):
    """
    Hourly Fine Fuel Moisture Code Calculation

    Parameters
    ----------
    temp: float
        Temperature (centigrade)
    rh: float
        Relative Humidity (%)
    ws: float
        Wind speed (km/h)
    fo: float, default=85
        FFMC at the previous time step
    t0: float, default=1
        Time (in hours) between the previous value of FFMC and the current time at which we want to
        calculate a new value of the FFMC.

    Returns
    -------
    float
        Fine Fuel Moisture Code

    Notes
    -----
    Hourly Fine Fuel Moisture Code is based on a calculation routine first described in detail by Van
    Wagner (1977) and which has been updated in minor ways by the Canadian Forest
    Service to have it agree with the calculation methodology for the daily FFMC.
    In its simplest typical use this current routine calculates a value of FFMC based on a series of
    uninterrupted hourly weather observations of screen level (~1.4 m) temperature, relative humidity, 10 m
    wind speed, and 1-hour rainfall. This implementation of the function
    includes an optional time step input which is defaulted to one hour, but can
    be reduced if sub-hourly calculation of the code is needed.  The FFMC is in
    essence a bookkeeping system for moisture content and thus it needs to use
    the last time step's value of FFMC in its calculation as well.

    The hourly FFMC is very similar in its structure and calculation to the
    Canadian Forest Fire Weather Index System's daily FFMC
    but has an altered drying and wetting rate which more realistically reflects
    the drying and wetting of a pine needle litter layer sitting on a decaying
    organic layer.  This particular implementation of the Canadian Forest Fire
    Danger Rating System's hourly FFMC provides for a flexible time step; that
    is, the data need not necessarily be in time increments of one hour.  This
    flexibility has been added for some users who use this method with data
    sampled more frequently that one hour.  We do not recommend using a time
    step much greater than one hour. An important and implicit assumption in
    this calculation is that the input weather is constant over the time step of
    each calculation (e.g., typically over the previous hour).  This is a
    reasonable assumption for an hour; however it can become problematic for
    longer periods.  For brevity we have referred to this routine throughout
    this description as the hourly FFMC.

    """
    if rh < 0 or rh > 100:
        raise ValueError(f"Invalid rh: {rh}")
    if prec < 0:
        raise ValueError(f"Invalid prec: {prec}")
    if ws < 0:
        raise ValueError(f"Invalid ws: {ws}")
    # Eq. 1 (with a more precise multiplier than the daily)
    mo = FFMC_COEFFICIENT * (101 - fo) / (59.5 + fo)
    rf = prec
    # Eqs. 3a & 3b (Van Wagner & Pickett 1985)
    if mo <= 150:
        mr = mo + 42.5 * rf * exp(-100 / (251 - mo)) * (1 - exp(-6.93 / rf))
    else:
        mr = (
            mo
            + 42.5 * rf * exp(-100 / (251 - mo)) * (1 - exp(-6.93 / rf))
            + 0.0015 * ((mo - 150) ** 2) * (rf**0.5)
        )
    # The real moisture content of pine litter ranges up to about 250 percent,
    # so we cap it at 250
    mr = min(mr, 250)
    mo = mr if prec > 0.0 else mo
    # Eq. 2a Equilibrium moisture content from drying
    ed = (
        0.942 * (rh**0.679)
        + 11 * exp((rh - 100) / 10)
        + 0.18 * (21.1 - temp) * (1 - exp(-0.115 * rh))
    )
    # Eq. 3a Log drying rate at the normal temperature of 21.1C
    ko = 0.424 * (1 - (rh / 100) ** 1.7) + 0.0694 * (ws**0.5) * (1 - (rh / 100) ** 8)
    # Eq. 3b
    kd = ko * 0.0579 * exp(0.0365 * temp)
    # Eq. 8 (Van Wagner & Pickett 1985)
    md = ed + (mo - ed) * (10 ** (-kd * t0))
    # Eq. 2b Equilibrium moisture content from wetting
    ew = (
        0.618 * (rh**0.753)
        + 10 * exp((rh - 100) / 10)
        + 0.18 * (21.1 - temp) * (1 - exp(-0.115 * rh))
    )
    # Eq. 7a Log wetting rate at the normal temperature of 21.1 C
    k1 = 0.424 * (1 - ((100 - rh) / 100) ** 1.7) + 0.0694 * (ws**0.5) * (
        1 - ((100 - rh) / 100) ** 8
    )
    # Eq. 4b
    kw = k1 * 0.0579 * exp(0.0365 * temp)
    # Eq. 8 (Van Wagner & Pickett 1985)
    mw = ew - (ew - mo) * (10 ** (-kw * t0))
    # Constraints
    if mo > ed:
        m = md
    elif ed >= mo >= ew:
        m = mo
    else:
        m = mw
    # Eq. 6 - Final hffmc calculation
    fo = 59.5 * (250 - m) / (FFMC_COEFFICIENT + m)
    fo = max(fo, 0)
    return fo
