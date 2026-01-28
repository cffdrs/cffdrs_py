import math


def foliar_moisture_content(lat, long, elv, dj, d0):
    """
    Foliar Moisture Content Calculator

    Calculate Foliar Moisture Content on a specified day.
    All variables names are laid out in the same manner as Forestry Canada
    Fire Danger Group (FCFDG) (1992). Development and Structure of the
    Canadian Forest Fire Behavior Prediction System." Technical Report
    ST-X-3, Forestry Canada, Ottawa, Ontario.

    :param lat: Latitude (decimal degrees)
    :param long: Longitude (decimal degrees)
    :param elv: Elevation (metres)
    :param dj: Day of year (often referred to as julian date)
    :param d0: Date of minimum foliar moisture content. If D0, date of min FMC, is not known then D0 = NULL.

    :returns: FMC Foliar Moisture Content value
    """
    # Initialize
    fmc = -1
    latn = 0
    # Calculate Normalized Latitude
    # Eqs. 1 & 3 (FCFDG 1992)
    if d0 <= 0:
        if elv <= 0:
            latn = 46 + 23.4 * math.exp(-0.0360 * (150 - long))
        else:
            latn = 43 + 33.7 * math.exp(-0.0351 * (150 - long))
    # Calculate Date of minimum foliar moisture content
    # Eqs. 2 & 4 (FCFDG 1992)
    if d0 <= 0:
        if elv <= 0:
            d0 = 151 * (lat / latn)
        else:
            d0 = 142.1 * (lat / latn) + 0.0172 * elv
    # Round D0 to the nearest integer because it is a date
    d0 = round(d0, 0)
    # Number of days between day of year and date of min FMC
    # Eq. 5 (FCFDG 1992)
    nd = abs(dj - d0)
    # Calculate final FMC
    # Eqs. 6, 7, & 8 (FCFDG 1992)
    if nd < 30:
        fmc = 85 + 0.0189 * nd**2
    elif nd >= 30 and nd < 50:
        fmc = 32.9 + 3.17 * nd - 0.0288 * nd**2
    else:
        fmc = 120
    return fmc
