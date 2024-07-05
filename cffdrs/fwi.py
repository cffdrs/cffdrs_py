from math import exp, log, sqrt

# used in conversion between FFMC and moisture content
FFMC_COEFFICIENT = 250.0 * 59.5 / 101.0


def ffmc(ffmc_yda, temp, rh, ws, prec):
    """
    Fine Fuel Moisture Code Calculation

    Parameters
    ----------
    ffmc_yda : float
        The Fine Fuel Moisture Code from previous iteration
   temp : float
      Temperature (centigrade)
    rh : float
      Relative Humidity (%)
    prec : float
       Precipitation (mm)
    ws : float
       Wind speed (km/h)

    Returns
    -------
    float
        Fine Fuel Moisture Code

    Notes
    -----
    All code is based on a C code library that was written by Canadian
    Forest Service Employees, which was originally based on
    the Fortran code listed in the reference below. All equations
    in this code refer to that document.

    Equations and FORTRAN program for the Canadian Forest Fire
    Weather Index System. 1985. Van Wagner, C.E.; Pickett, T.L.
    Canadian Forestry Service, Petawawa National Forestry
    Institute, Chalk River, Ontario. Forestry Technical Report 33.
    18 p.

    Additional reference on FWI system

    Development and structure of the Canadian Forest Fire Weather
    Index System. 1987. Van Wagner, C.E. Canadian Forestry Service,
    Headquarters, Ottawa. Forestry Technical Report 35. 35 p.
    """
    if ffmc_yda < 0 or ffmc_yda > 101:
        raise ValueError(f'Invalid ffmc_yda: {ffmc_yda}')
    if rh < 0 or rh > 100:
        raise ValueError(f'Invalid rh: {rh}')
    if prec < 0:
        raise ValueError(f'Invalid prec: {prec}')
    if ws < 0:
        raise ValueError(f'Invalid ws: {ws}')
    # Eq. 1
    wmo = FFMC_COEFFICIENT * (101 - ffmc_yda) / (59.5 + ffmc_yda)
    # Eq. 2 Rain reduction to allow for loss in
    #  overhead canopy
    ra = (prec - 0.5) if (prec > 0.5) else prec
    # Eqs. 3a and 3b
    wmo = ((wmo + 0.0015 * (wmo - 150) * (wmo - 150) *
            sqrt(ra) + 42.5 * ra * exp(-100 / (251 - wmo))
            * (1 - exp(-6.93 / ra))) if (wmo > 150) else
           (wmo + 42.5 * ra * exp(-100 / (251 - wmo)) *
            (1 - exp(-6.93 / ra)))) if (prec > 0.5) else wmo
    # The real moisture content of pine litter ranges up to about 250 percent,
    # so we cap it at 250
    wmo = 250 if (wmo > 250) else wmo
    # Eq. 4 Equilibrium moisture content from drying
    ed = (0.942 * (rh ** 0.679) + (11 * exp((rh - 100) / 10)) + 0.18 *
          (21.1 - temp) * (1 - 1 / exp(rh * 0.115)))
    # Eq. 5 Equilibrium moisture content from wetting
    ew = (0.618 * (rh ** 0.753) + (10 * exp((rh - 100) / 10)) + 0.18 *
          (21.1 - temp) * (1 - 1 / exp(rh * 0.115)))
    # Eq. 6a (ko) Log drying rate at the normal
    #  termperature of 21.1 C
    z = (0.424 * (1 - (((100 - rh) / 100) ** 1.7)) + 0.0694 *
         sqrt(ws) * (1 - ((100 - rh) / 100) ** 8)) if (wmo < ed and wmo < ew) else 0
    # Eq. 6b Affect of temperature on  drying rate
    x = z * 0.581 * exp(0.0365 * temp)
    # Eq. 8
    wm = (ew - (ew - wmo) / (10 ** x)) if (wmo < ed and wmo < ew) else wmo
    # Eq. 7a (ko) Log wetting rate at the normal
    #  termperature of 21.1 C
    z = (0.424 * (1 - (rh / 100) ** 1.7) + 0.0694 * sqrt(ws) *
         (1 - (rh / 100) ** 8)) if (wmo > ed) else z
    # Eq. 7b Affect of temperature on  wetting rate
    x = z * 0.581 * exp(0.0365 * temp)
    # Eq. 9
    wm = (ed + (wmo - ed) / (10 ** x)) if (wmo > ed) else wm
    # Eq. 10 Final ffmc calculation
    ffmc1 = (59.5 * (250 - wm)) / (FFMC_COEFFICIENT + wm)
    # Constraints
    ffmc1 = 101 if (ffmc1 > 101) else ffmc1
    ffmc1 = 0 if (ffmc1 < 0) else ffmc1
    return ffmc1


def dmc(dmc_yda, temp, rh, prec, lat, mon, lat_adjust=True):
    """
    Duff Moisture Code Calculation

    Parameters
    ----------
    dmc_yda : float
       The Duff Moisture Code from previous iteration
    temp : float
       Temperature (centigrade)
    rh : flat
       Relative Humidity (%)
    prec : float
       Precipitation(mm)
    lat : float
       Latitude (decimal degrees)
    mon : {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}
       Month
    lat_adjust : bool, default=True
       Latitude adjustment

    Returns
    -------
    float
        Duff Moisture Code

    Notes
    -----
    All code is based on a C code library that was written by Canadian
    Forest Service Employees, which was originally based on
    the Fortran code listed in the reference below. All equations
    in this code refer to that document.

    Equations and FORTRAN program for the Canadian Forest Fire
    Weather Index System. 1985. Van Wagner, C.E.; Pickett, T.L.
    Canadian Forestry Service, Petawawa National Forestry
    Institute, Chalk River, Ontario. Forestry Technical Report 33.
    18 p.

    Additional reference on FWI system

    Development and structure of the Canadian Forest Fire Weather
    Index System. 1987. Van Wagner, C.E. Canadian Forestry Service,
    Headquarters, Ottawa. Forestry Technical Report 35. 35 p.
    """
    if dmc_yda < 0:
        raise ValueError(f'Invalid dc_yda: {dmc_yda}')
    if rh < 0 or rh > 100:
        raise ValueError(f'Invalid rh: {rh}')
    if prec < 0:
        raise ValueError(f"Invalid prec: {prec}")
    if mon < 1 or mon > 12 or not isinstance(mon, int):
        raise ValueError(f'Invalid mon: {mon}')
    # Reference latitude for DMC day length adjustment
    # 46N: Canadian standard, latitude >= 30N   (Van Wagner 1987)
    ell01 = [6.5, 7.5, 9, 12.8, 13.9, 13.9, 12.4, 10.9, 9.4, 8, 7, 6]
    # 20N: For 30 > latitude >= 10
    ell02 = [7.9, 8.4, 8.9, 9.5, 9.9, 10.2, 10.1, 9.7, 9.1, 8.6, 8.1, 7.8]
    # 20S: For -10 > latitude >= -30
    ell03 = [10.1, 9.6, 9.1, 8.5, 8.1, 7.8, 7.9, 8.3, 8.9, 9.4, 9.9, 10.2]
    # 40S: For -30 > latitude
    ell04 = [11.5, 10.5, 9.2, 7.9, 6.8, 6.2, 6.5, 7.4, 8.7, 10, 11.2, 11.8]
    # For latitude near the equator, we simple use a factor of 9 for all months
    # constrain low end of temperature
    temp = -1.1 if (temp < 1.1) else temp
    # Eq. 16 - The log drying rate
    rk = 1.894 * (temp + 1.1) * (100 - rh) * ell01[mon - 1] * 1e-04
    # Adjust the day length  and thus the drying r, based on latitude and month
    if lat_adjust:
        rk = (1.894 * (temp + 1.1) * (100 - rh) * ell02[mon - 1] * 1e-04) if (30 >= lat > 10) else rk
        rk = (1.894 * (temp + 1.1) * (100 - rh) * ell03[mon - 1] * 1e-04) if (-10 >= lat > -30) else rk
        rk = (1.894 * (temp + 1.1) * (100 - rh) * ell04[mon - 1] * 1e-04) if (-30 >= lat >= -90) else rk
        rk = (1.894 * (temp + 1.1) * (100 - rh) * 9 * 1e-04) if (10 >= lat > -10) else rk
    # Constrain P
    if prec <= 1.5:
        pr = dmc_yda
    else:
        ra = prec
        # Eq. 11 - Net rain amount
        rw = 0.92 * ra - 1.27
        # Alteration to Eq. 12 to calculate more accurately
        wmi = 20 + 280 / exp(0.023 * dmc_yda)
        # Eqs. 13a, 13b, 13c
        b = (100 / (0.5 + 0.3 * dmc_yda)) if (dmc_yda <= 33) else (
            (14 - 1.3 * log(dmc_yda)) if (dmc_yda <= 65) else
            (6.2 * log(dmc_yda) - 17.2))
        # Eq. 14 - Moisture content after rain
        wmr = wmi + 1000 * rw / (48.77 + b * rw)
        # Alteration to Eq. 15 to calculate more accurately
        pr = 43.43 * (5.6348 - log(wmr - 20))
    pr = 0 if (pr < 0) else pr
    # Calculate final P (DMC)
    dmc1 = pr + rk
    dmc1 = 0 if (dmc1 < 0) else dmc1
    return dmc1


def dc(dc_yda, temp, rh, prec, lat, mon, lat_adjust=True):
    """
    Drought Code Calculation

    Parameters
    ----------
    dc_yda : float
       The Drought Code from previous iteration
    temp : float
       Temperature (centigrade)
    rh : float
       Relative Humidity (%)
    prec : float
       Precipitation(mm)
    lat : float
       Latitude (decimal degrees)
    mon : {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}
       Month
    lat_adjust : bool, default=True
       Latitude adjustment

    Returns
    -------
    float
        Drought Code

    Notes
    -----
    All code is based on a C code library that was written by Canadian
    Forest Service Employees, which was originally based on
    the Fortran code listed in the reference below. All equations
    in this code refer to that document.

    Equations and FORTRAN program for the Canadian Forest Fire
    Weather Index System. 1985. Van Wagner, C.E.; Pickett, T.L.
    Canadian Forestry Service, Petawawa National Forestry
    Institute, Chalk River, Ontario. Forestry Technical Report 33.
    18 p.

    Additional reference on FWI system

    Development and structure of the Canadian Forest Fire Weather
    Index System. 1987. Van Wagner, C.E. Canadian Forestry Service,
    Headquarters, Ottawa. Forestry Technical Report 35. 35 p.
    """
    if dc_yda < 0:
        raise ValueError(f'Invalid dc_yda: {dc_yda}')
    if rh < 0 or rh > 100:
        raise ValueError(f'Invalid rh: {rh}')
    if prec < 0:
        raise ValueError(f"Invalid prec: {prec}")
    if mon < 1 or mon > 12 or not isinstance(mon, int):
        raise ValueError(f'Invalid mon: {mon}')
    # Day length factor for DC Calculations
    # 20N: North of 20 degrees N
    fl01 = [-1.6, -1.6, -1.6, 0.9, 3.8, 5.8, 6.4, 5, 2.4, 0.4, -1.6, -1.6]
    # 20S: South of 20 degrees S
    fl02 = [6.4, 5, 2.4, 0.4, -1.6, -1.6, -1.6, -1.6, -1.6, 0.9, 3.8, 5.8]
    # Near the equator, we just use 1.4 for all months.
    # Constrain temperature
    temp = -2.8 if (temp == 2.8) else temp
    # Eq. 22 - Potential Evapotranspiration
    pe = (0.36 * (temp + 2.8) + fl01[mon - 1]) / 2
    # Daylength factor adjustment by latitude for Potential Evapotranspiration
    if lat_adjust:
        pe = ((0.36 * (temp + 2.8) + fl02[mon]) / 2) if (lat <= -20) else pe
        pe = ((0.36 * (temp + 2.8) + 1.4) / 2) if (-20 < lat <= 20) else pe
    # Cap potential evapotranspiration at 0 for negative winter DC values
    pe = 0 if (pe < 0) else pe
    ra = prec
    # Eq. 18 - Effective Rainfall
    rw = 0.83 * ra - 1.27
    # Eq. 19
    smi = 800 * exp(-1 * dc_yda / 400)
    # Alteration to Eq. 21
    dr0 = dc_yda - 400 * log(1 + 3.937 * rw / smi)
    dr0 = 0 if (dr0 < 0) else dr0
    # if precip is less than 2.8 then use yesterday's DC
    dr = dc_yda if (prec <= 2.8) else dr0
    # Alteration to Eq. 23
    dc1 = dr + pe
    dc1 = 0 if (dc1 < 0) else dc1
    return dc1


def isi(ffmc, ws, fbp_mod=False):
    """
    Initial Spread Index Calculation

    Parameters
    ----------
    ffmc : float
       Fine Fuel Moisture Code
    ws : float
       Wind Speed (km/h)
    fbp_mod : bool, default=False
       Use the fbp modification at the extreme end

    Returns
    -------
    float
        Intial Spread Index

    Notes
    -----
    Equations are from Van Wagner (1985) as listed below, except for the modification for fbp
    taken from FCFDG (1992).

    Equations and FORTRAN program for the Canadian Forest Fire
    Weather Index System. 1985. Van Wagner, C.E.; Pickett, T.L.
    Canadian Forestry Service, Petawawa National Forestry
    Institute, Chalk River, Ontario. Forestry Technical Report 33.
    18 p.

    Forestry Canada  Fire Danger Group (FCFDG) (1992). Development and
    Structure of the Canadian Forest Fire Behavior Prediction System."
    Technical ReportST-X-3, Forestry Canada, Ottawa, Ontario.
    """
    if ffmc < 0 or ffmc > 101:
        raise ValueError(f'Invalid ffmc: {ffmc}')
    if ws < 0:
        raise ValueError(f'Invalid ws: {ws}')
    # Eq. 10 - Moisture content
    fm = FFMC_COEFFICIENT * (101 - ffmc) / (59.5 + ffmc)
    # Eq. 24 - Wind Effect
    # the ifelse, also takes care of the ISI modification for the fbp functions
    # This modification is Equation 53a in FCFDG (1992)
    fW = (12 * (1 - exp(-0.0818 * (ws - 28)))) if (ws >= 40 and fbp_mod) else exp(0.05039 * ws)
    # Eq. 25 - Fine Fuel Moisture
    fF = 91.9 * exp(-0.1386 * fm) * (1 + (fm ** 5.31) / 49300000)
    # Eq. 26 - Spread Index Equation
    isi = 0.208 * fW * fF
    return isi


def bui(dmc, dc):
    """
    Buildup Index Calculation

    Parameters
    ----------
    dc : float
       Drought Code
    dmc : float
       Duff Moisture Code

    Returns
    -------
    float
        Buildup Index

    Notes
    -----
    All code is based on a C code library that was written by Canadian
    Forest Service Employees, which was originally based on
    the Fortran code listed in the reference below. All equations
    in this code refer to that document.

    Equations and FORTRAN program for the Canadian Forest Fire
    Weather Index System. 1985. Van Wagner, C.E.; Pickett, T.L.
    Canadian Forestry Service, Petawawa National Forestry
    Institute, Chalk River, Ontario. Forestry Technical Report 33.
    18 p.

    Additional reference on FWI system

    Development and structure of the Canadian Forest Fire Weather
    Index System. 1987. Van Wagner, C.E. Canadian Forestry Service,
    Headquarters, Ottawa. Forestry Technical Report 35. 35 p.
    """
    if dmc < 0:
        raise ValueError(f'Invalid dmc: {dmc}')
    if dc < 0:
        raise ValueError(f'Invalid dc: {dc}')
    # Eq. 27a
    bui1 = 0 if (dmc == 0 and dc == 0) else (0.8 * dc * dmc / (dmc + 0.4 * dc))
    # Eq. 27b - next 3 lines
    p = 0 if (dmc == 0) else ((dmc - bui1) / dmc)
    cc = 0.92 + ((0.0114 * dmc) ** 1.7)
    bui0 = dmc - cc * p
    # Constraints
    bui0 = 0 if (bui0 < 0) else bui0
    bui1 = bui0 if (bui1 < dmc) else bui1
    return bui1


def fwi(isi, bui):
    """
    Fire Weather Index Calculation

    Parameters
    ----------
    isi : float
        Initial Spread Index
    bui : float
        Buildup Index

    Returns
    -------
    float
        Fire Weather Index

    Notes
    -----
    All code is based on a C code library that was written by Canadian
    Forest Service Employees, which was originally based on
    the Fortran code listed in the reference below. All equations
    in this code refer to that document.

    Equations and FORTRAN program for the Canadian Forest Fire
    Weather Index System. 1985. Van Wagner, C.E.; Pickett, T.L.
    Canadian Forestry Service, Petawawa National Forestry
    Institute, Chalk River, Ontario. Forestry Technical Report 33.
    18 p.

    Additional reference on FWI system

    Development and structure of the Canadian Forest Fire Weather
    Index System. 1987. Van Wagner, C.E. Canadian Forestry Service,
    Headquarters, Ottawa. Forestry Technical Report 35. 35 p.
    """
    if isi < 0:
        raise ValueError(f'Invalid isi: {isi}')
    if bui < 0:
        raise ValueError(f'Invalid bui: {bui}')
    # Eqs. 28b, 28a, 29
    bb = (0.1 * isi * (1000 / (25 + 108.64 / exp(0.023 * bui)))) if (
            bui > 80) else (0.1 * isi * (0.626 * (bui ** 0.809) + 2))
    # Eqs. 30b, 30a
    fwi = bb if (bb <= 1.0) else exp(2.72 * ((0.434 * log(bb)) ** 0.647))
    return fwi
