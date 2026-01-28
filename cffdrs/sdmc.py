import math


def sdmc(temp, rh, ws, prec, mon, dmc, sdmc_old=None):
    """
    Sheltered Duff Moisture Code

    sdmc is used to calculate sheltered DMC (sDMC, Wotton et al., 2005) based on daily noon weather
    observations of temperature, relative humidity, wind speed, 24-hour rainfall,
    and a previous day's calculated or estimated value of sDMC.

    The Duff Moisture Code (DMC) component of the Canadian Forest Fire Weather
    Index (FWI) System tracks moisture content of the forest floor away from the
    sheltering influences of overstory trees.  This sheltered Duff Moisture Code
    (sDMC) was developed to track moisture in the upper 5 cm of the organic
    layer in the rain sheltered areas near (<0.5 m) the boles of overstory trees
    (Wotton et al. 2005), an area where lightning strikes usually ignite the
    forest floor when they run to ground. The sDMC is very similar in structure
    (and identical in data requirements) to the DMC.  The sDMC, like all the FWI
    System moisture codes, is a bookkeeping system that tracks gain and loss of
    moisture from day-to-day; thus an estimate of the previous day's sDMC value
    is needed to provide a starting point for each day's moisture calculation.
    Like the other moisture codes in the FWI System the sDMC is converted from a
    moisture content value to an outputted CODE value which increases in value
    with decreasing moisture content.

    Args:
        temp: Temperature (centigrade)
        rh: Relative humidity (%)
        ws: 10-m height wind speed (km/h)
        prec: 1-hour rainfall (mm)
        mon: Month of the observations (integer 1-12)
        dmc: Duff Moisture Code
        sdmc_old: Previous day's value of SDMC. At the start of calculations,
        when there is no calculated previous day's SDMC value to use, the user must
        specify an estimate of this value. Where sdmc_old=None, the function
        will calculate the initial SDMC value based on the initial DMC.

    Returns:
        sdmc returns the SDMC value.

    References:
        Wotton, B.M., B.J. Stocks, and D.L. Martell. 2005. An index for
        tracking sheltered forest floor moisture within the Canadian Forest Fire
        Weather Index System. International Journal of Wildland Fire, 14, 169-182.

        Equations and FORTRAN program for the Canadian Forest Fire
        Weather Index System. 1985. Van Wagner, C.E.; Pickett, T.L.
        Canadian Forestry Service, Petawawa National Forestry
        Institute, Chalk River, Ontario. Forestry Technical Report 33.
        18 p.
    """
    # Constrain rh, ws and precipitation
    rh = min(99.9, max(10.0, rh))
    ws = max(0.0, ws)
    prec = max(0.0, prec)

    # Reference latitude for DMC day length adjustment
    el = [6.5, 7.5, 9.0, 12.8, 13.9, 13.9, 12.4, 10.9, 9.4, 8.0, 7.0, 6.0]

    # Initialize sdmc if it does not exist
    if sdmc_old is None:
        sdmc_old = 2.6 + (1.7 * dmc)
        sdmc_old -= 6.0
        sdmc_old = max(12.0, sdmc_old)

    # Constrain temperature
    t0 = max(-1.1, temp)

    # This is a modification multiplier at front
    rk = 4.91 / 3.57 * 1.894 * (t0 + 1.1) * (100 - rh) * el[mon - 1] * 0.0001

    # Eq.7 (Wotton et. al. 2005) calculates rain throughfall.
    rw = 0.218 * prec - 0.094 if prec < 7.69 else 0.83 * prec - 4.8
    rw = max(0.0, rw)  # prevent negative throughfall (matches R)

    # Alteration to Eq. 12 (Van Wagner & Pickett 1985)
    wmi = 20.0 + 280.0 / math.exp(0.023 * sdmc_old)

    # Eqs. 13a, 13b, 13c (Van Wagner & Pickett 1985)
    if sdmc_old <= 33:
        b = 100.0 / (0.5 + 0.3 * sdmc_old)
    elif sdmc_old <= 65:
        b = 14.0 - 1.3 * math.log(sdmc_old)
    else:
        b = 6.2 * math.log(sdmc_old) - 17.2

    # Eq. 14 (Van Wagner & Pickett 1985) - Moisture content after rain
    wmr = wmi + 1000.0 * rw / (48.77 + b * rw)

    # Alteration to Eq. 15 (Van Wagner & Pickett 1985)
    if prec <= 0.44:
        pr = sdmc_old
    else:
        # Ensure argument for log is > 0 to avoid math error
        pr = 43.43 * (5.6348 - math.log(max(wmr - 20, 0.0001)))

    # Constrain p
    pr = max(0.0, pr)

    # Calculate final SDMC
    sdmc = pr + rk

    # Constrain result
    sdmc = max(0.0, sdmc)

    return sdmc
