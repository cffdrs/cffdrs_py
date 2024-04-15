from math import exp, log, sqrt
import numpy.ma as np


def ffmc_arr(ffmc_yda, temp, rh, ws, prec):
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
    if np.any(ffmc_yda < 0) or np.any(ffmc_yda > 101):
        raise ValueError(f'Invalid ffmc_yda: {ffmc_yda}')
    if np.any(rh < 0) or np.any(rh > 100):
        raise ValueError(f'Invalid rh: {rh}')
    if np.any(prec < 0):
        raise ValueError(f'Invalid prec: {prec}')
    if np.any(ws < 0):
        raise ValueError(f'Invalid ws: {ws}')
    # Eq. 1
    wmo = 147.2 * (101 - ffmc_yda) / (59.5 + ffmc_yda)
    # Eq. 2 Rain reduction to allow for loss in
    #  overhead canopy
    ###
    # Non-array operations
    ###
    # ra = (prec - 0.5) if (prec > 0.5) else prec
    ra = np.where(prec > 0.5, prec - 0.5, prec)
    # Eqs. 3a and 3b
    ###
    # Non-array operations
    ###    
    # wmo = ((wmo + 0.0015 * (wmo - 150) * (wmo - 150) *
    #         sqrt(ra) + 42.5 * ra * exp(-100 / (251 - wmo))
    #         * (1 - exp(-6.93 / ra))) if (wmo > 150) else
    #        (wmo + 42.5 * ra * exp(-100 / (251 - wmo)) *
    #         (1 - exp(-6.93 / ra)))) if (prec > 0.5) else wmo
    
    wmo = eq3a3b(wmo, ra, prec)
    
    # The real moisture content of pine litter ranges up to about 250 percent,
    # so we cap it at 250
    
    ###
    # Non-array operations
    ###
    # wmo = 250 if (wmo > 250) else wmo
    
    wmo = np.where(wmo > 250, 250, wmo)
    
    # Eq. 4 Equilibrium moisture content from drying
    ed = (0.942 * (rh ** 0.679) + (11 * np.exp((rh - 100) / 10)) + 0.18 *
          (21.1 - temp) * (1 - 1 / np.exp(rh * 0.115)))
    # Eq. 5 Equilibrium moisture content from wetting
    ew = (0.618 * (rh ** 0.753) + (10 * np.exp((rh - 100) / 10)) + 0.18 *
          (21.1 - temp) * (1 - 1 / np.exp(rh * 0.115)))
    # Eq. 6a (ko) Log drying rate at the normal
    #  termperature of 21.1 C
    ###
    # Non-array operations
    ###
    # z = (0.424 * (1 - (((100 - rh) / 100) ** 1.7)) + 0.0694 *
    #      sqrt(ws) * (1 - ((100 - rh) / 100) ** 8)) if (wmo < ed and wmo < ew) else 0
    
    z = np.where(wmo < ed,
                 np.where(wmo < ew,
                    (0.424 * (1 - (((100 - rh) / 100) ** 1.7)) + 0.0694 * np.sqrt(ws) * (1 - ((100 - rh) / 100) ** 8)),
                    0), 
                 0)
    
    # Eq. 6b Affect of temperature on  drying rate
    x = z * 0.581 * np.exp(0.0365 * temp)
    # Eq. 8
    ###
    # Non-array operations
    ###
    # wm = (ew - (ew - wmo) / (10 ** x)) if (wmo < ed and wmo < ew) else wmo
    
    wm = np.where(wmo < ed,
                  np.where(wmo < ew, 
                    (ew - (ew - wmo) / (10 ** x)), 
                    wmo),
                  wmo)
    
    # Eq. 7a (ko) Log wetting rate at the normal
    #  termperature of 21.1 C
    ###
    # Non-array operations
    ###
    # z = (0.424 * (1 - (rh / 100) ** 1.7) + 0.0694 * sqrt(ws) *
    #      (1 - (rh / 100) ** 8)) if (wmo > ed) else z
    
    z = np.where(wmo > ed, 
            (0.424 * (1 - (rh / 100) ** 1.7) + 0.0694 * np.sqrt(ws) * (1 - (rh / 100) ** 8)), 
            z)
    
    # Eq. 7b Affect of temperature on  wetting rate
    x = z * 0.581 * np.exp(0.0365 * temp)
    # Eq. 9
    ###
    # Non-array operations
    ###
    # wm = (ed + (wmo - ed) / (10 ** x)) if (wmo > ed) else wm
    
    wm = np.where(wmo > ed, 
                  (ed + (wmo - ed) / (10 ** x)),
                  wm)
    
    # Eq. 10 Final ffmc calculation
    ffmc1 = (59.5 * (250 - wm)) / (147.2 + wm)
    # Constraints
    ###
    # Non-array operations
    ###
    # ffmc1 = 101 if (ffmc1 > 101) else ffmc1
    
    ffmc1 = np.where(ffmc1 > 101, 101, ffmc1)
    
    ###
    # Non-array operations
    ###
    # ffmc1 = 0 if (ffmc1 < 0) else ffmc1
    
    ffmc1 = np.where(ffmc1 < 0, 0, ffmc1)
    
    return ffmc1

def eq3a3b(wmo, ra, prec):
    # test = np.sqrt(ra)
    # test2 = np.exp(-100 / (251 - wmo))
    # test3 = (1 - np.exp(-6.93 / ra))
    # test4 = np.exp(-6.93 / ra)
    wmo1 = (wmo + 0.0015 * (wmo - 150) * (wmo - 150) *
            np.sqrt(ra) + 42.5 * ra * np.exp(-100 / (251 - wmo))
            * (1 - np.exp(-6.93 / ra))) 
    wmo2 = (wmo + 42.5 * ra * np.exp(-100 / (251 - wmo)) *
            (1 - np.exp(-6.93 / ra)))
    wmo = np.where(prec > 0.5,
                   np.where(wmo > 150, 
                        wmo1, 
                        wmo2),
                   wmo
    )
    return wmo

def dmc_arr(dmc_yda, temp, rh, prec, lat, mon, lat_adjust=True):
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
    mon : {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}
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
    if np.any(dmc_yda < 0):
        raise ValueError(f'Invalid dc_yda: {dmc_yda}')
    if np.any(rh < 0) or np.any(rh > 100):
        raise ValueError(f'Invalid rh: {rh}')
    if np.any(prec < 0):
        raise ValueError(f'Invalid prec: {prec}')
    # if ws < 0:
    #     raise ValueError(f'Invalid ws: {ws}')
    if mon < 0 or mon > 11 or not isinstance(mon, int):
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
    ###
    # Non-array operations
    ###
    # temp = -1.1 if (temp == 1.1) else temp 

    temp = np.where(temp == 1.1, -1.1, temp)
    
    # Eq. 16 - The log drying rate
    rk = 1.894 * (temp + 1.1) * (100 - rh) * ell01[mon] * 1e-04
    # Adjust the day length  and thus the drying r, based on latitude and month
    if lat_adjust:
        ###
        # Non-array operations
        ###
        
        # rk = (1.894 * (temp + 1.1) * (100 - rh) * ell02[mon] * 1e-04) if (30 >= lat > 10) else rk
        # rk = (1.894 * (temp + 1.1) * (100 - rh) * ell03[mon] * 1e-04) if (-10 >= lat > -30) else rk
        # rk = (1.894 * (temp + 1.1) * (100 - rh) * ell04[mon] * 1e-04) if (-30 >= lat >= -90) else rk
        # rk = (1.894 * (temp + 1.1) * (100 - rh) * 9 * 1e-04) if (10 >= lat > -10) else rk
        # rk = np.where(lat > 10 and lat <= 30, (1.894 * (temp + 1.1) * (100 - rh) * ell02[mon] * 1e-04), rk)
        rk = np.where(lat > 10, 
            np.where(lat <= 30, (1.894 * (temp + 1.1) * (100 - rh) * ell02[mon] * 1e-04), rk),
            rk) 
        rk = np.where(lat > -30, 
            np.where(lat <= -10, (1.894 * (temp + 1.1) * (100 - rh) * ell03[mon] * 1e-04), rk),
            rk) 
        rk = np.where(lat >= -90, 
            np.where(lat <= -30, (1.894 * (temp + 1.1) * (100 - rh) * ell04[mon] * 1e-04), rk),
            rk) 
        rk = np.where(lat > -10, 
            np.where(lat <= 10, (1.894 * (temp + 1.1) * (100 - rh) * 9 * 1e-04), rk),
            rk) 
   
    ###
    # Non-array operations
    ###
    # Constrain P
    # if prec <= 1.5:
    #     pr = dmc_yda
    # else:
    #     ra = prec
    #     # Eq. 11 - Net rain amount
    #     rw = 0.92 * ra - 1.27
    #     # Alteration to Eq. 12 to calculate more accurately
    #     wmi = 20 + 280 / exp(0.023 * dmc_yda)
    #     # Eqs. 13a, 13b, 13c
    #     b = (100 / (0.5 + 0.3 * dmc_yda)) if (dmc_yda <= 33) else (
    #         (14 - 1.3 * log(dmc_yda)) if (dmc_yda <= 65) else
    #         (6.2 * log(dmc_yda) - 17.2))
    #     # Eq. 14 - Moisture content after rain
    #     wmr = wmi + 1000 * rw / (48.77 + b * rw)
    #     # Alteration to Eq. 15 to calculate more accurately
    #     pr = 43.43 * (5.6348 - log(wmr - 20))
        
    pr = np.where((prec <= 1.5), dmc_yda, calc_rain(prec, dmc_yda))
    
    ###
    # Non-array operations
    ###
    # pr = 0 if (pr < 0) else pr
    
    pr = np.where((pr < 0), 0, pr)
    
    # Calculate final P (DMC)
    dmc1 = pr + rk
    
    ###
    # Non-array operations
    ###
    # dmc1 = 0 if (dmc1 < 0) else dmc1
    
    dmc1 = np.where((dmc1 < 0), 0, dmc1)
    return dmc1

def calc_rain(prec, dmc_yda):
    ra = prec
    # Eq. 11 - Net rain amount
    rw = 0.92 * ra - 1.27
    # Alteration to Eq. 12 to calculate more accurately
    wmi = 20 + 280 / np.exp(0.023 * dmc_yda)
    # Eqs. 13a, 13b, 13c
    # b = (100 / (0.5 + 0.3 * dmc_yda)) if (dmc_yda <= 33) else (
    #         (14 - 1.3 * log(dmc_yda)) if (dmc_yda <= 65) else
    #         (6.2 * log(dmc_yda) - 17.2))
    
    term1 = (100 / (0.5 + 0.3 * dmc_yda))
    term2 = (14 - 1.3 * np.log(dmc_yda))
    term3 = (6.2 * np.log(dmc_yda) - 17.2)
    
    b = np.where((dmc_yda <= 65), 
                 np.where(dmc_yda <= 33, term1, term2),
                 term3)
    wmr = wmi + 1000 * rw / (48.77 + b * rw)
    # pr = 43.43 * (5.6348 - np.log(wmr - 20))
    # Decomposed the previous equation.
    wmr = wmr - 20
    pr = 43.43 * (5.6348 - np.log(wmr))
    return pr

def dc_arr(dc_yda, temp, rh, prec, lat, mon, lat_adjust=True):
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
    mon : {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}
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
    if np.any(dc_yda < 0):
        raise ValueError(f'Invalid dc_yda: {dc_yda}')
    if np.any(rh < 0) or np.any(rh > 100):
        raise ValueError(f'Invalid rh: {rh}')
    if np.any(prec < 0):
        raise ValueError(f'Invalid prec: {prec}')
    # if ws < 0:
    #     raise ValueError(f'Invalid ws: {ws}')
    if mon < 0 or mon > 11 or not isinstance(mon, int):
        raise ValueError(f'Invalid mon: {mon}')
    # Day length factor for DC Calculations
    # 20N: North of 20 degrees N
    fl01 = [-1.6, -1.6, -1.6, 0.9, 3.8, 5.8, 6.4, 5, 2.4, 0.4, -1.6, -1.6]
    # 20S: South of 20 degrees S
    fl02 = [6.4, 5, 2.4, 0.4, -1.6, -1.6, -1.6, -1.6, -1.6, 0.9, 3.8, 5.8]
    # Near the equator, we just use 1.4 for all months.
    # Constrain temperature
    
    ###
    # Non-array operations
    ###
    # temp = -2.8 if (temp == 2.8) else temp
    
    temp = np.where(temp == 2.8, -2.8, temp)
    
    # Eq. 22 - Potential Evapotranspiration
    pe = (0.36 * (temp + 2.8) + fl01[mon]) / 2
    # Daylength factor adjustment by latitude for Potential Evapotranspiration
    if lat_adjust:
        
        ###
        # Non-array operations
        ###
        # pe = ((0.36 * (temp + 2.8) + fl02[mon]) / 2) if (lat <= -20) else pe
        # pe = ((0.36 * (temp + 2.8) + 1.4) / 2) if (-20 < lat <= 20) else pe
        
        pe = np.where(lat <= -20, ((0.36 * (temp + 2.8) + fl02[mon]) / 2) ,pe)
        pe = np.where(lat > -20,
                      np.where(lat <= 20, ((0.36 * (temp + 2.8) + 1.4) / 2), pe),
                      pe)
    # Cap potential evapotranspiration at 0 for negative winter DC values
    
    ###
    # Non-array operations
    ###
    # pe = 0 if (pe < 0) else pe
    
    pe = np.where(pe < 0, 0, pe)
    
    ra = prec
    # Eq. 18 - Effective Rainfall
    rw = 0.83 * ra - 1.27
    # Eq. 19
    smi = 800 * np.exp(-1 * dc_yda / 400)
    # Alteration to Eq. 21
    dr0 = dc_yda - 400 * np.log(1 + 3.937 * rw / smi)
    
    ###
    # Non-array operations
    ###    
    # dr0 = 0 if (dr0 < 0) else dr0
    dr0 = np.where(dr0 < 0, 0, dr0)
    
    # if precip is less than 2.8 then use yesterday's DC
    ###
    # Non-array operations
    ### 
    # dr = dc_yda if (prec <= 2.8) else dr0
    
    dr = np.where(prec <= 2.8, dc_yda, dr0)
    # Alteration to Eq. 23
    dc1 = dr + pe
    
    ###
    # Non-array operations
    ### 
    # dc1 = 0 if (dc1 < 0) else dc1
    
    dc1 = np.where(dc1 < 0, 0, dc1)
    return dc1

def isi_arr(ffmc, ws, fbp_mod=False):
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
    if np.any(ffmc < 0) or np.any(ffmc > 101):
        raise ValueError(f'Invalid ffmc: {ffmc}')
    if np.any(ws < 0):
        raise ValueError(f'Invalid ws: {ws}')
    # Eq. 10 - Moisture content
    fm = 147.2 * (101 - ffmc) / (59.5 + ffmc)
    # Eq. 24 - Wind Effect
    # the ifelse, also takes care of the ISI modification for the fbp functions
    # This modification is Equation 53a in FCFDG (1992)
    ###
    # Non-array operations
    ### 
    # fW = (12 * (1 - exp(-0.0818 * (ws - 28)))) if (ws >= 40 and fbp_mod) else exp(0.05039 * ws)
    if fbp_mod:
        fW = np.where(ws >= 40, (12 * (1 - np.exp(-0.0818 * (ws - 28)))), np.exp(0.05039 * ws))
    else:
        fW = np.exp(0.05039 * ws)
    # Eq. 25 - Fine Fuel Moisture
    fF = 91.9 * np.exp(-0.1386 * fm) * (1 + (fm ** 5.31) / 49300000)
    # Eq. 26 - Spread Index Equation
    isi = 0.208 * fW * fF
    return isi

def bui_arr(dmc_loc, dc_loc):
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
    if np.any(dmc_loc < 0):
        raise ValueError(f'Invalid dmc: {dmc}')
    if np.any(dc_loc < 0):
        raise ValueError(f'Invalid dc: {dc}')
    # Eq. 27a
    
    ###
    # Non-array operations
    ### 
    
    # bui1 = 0 if (dmc == 0 and dc == 0) else (0.8 * dc * dmc / (dmc + 0.4 * dc))
    
    bui1 = np.where(dmc_loc == 0,
                    # np.where(dc_loc == 0, 0, (0.8 * dc_loc * dmc_loc / (dmc_loc + 0.4 * dc_loc)))
                    # (0.8 * dc_loc * dmc_loc / (dmc_loc + 0.4 * dc_loc)))
                    np.where(dc_loc == 0, 0, eq27a(dc_loc,dmc_loc)),
                    eq27a(dc_loc,dmc_loc))
    
    # Eq. 27b - next 3 lines
    ###
    # Non-array operations
    ### 
    # p = 0 if (dmc == 0) else ((dmc - bui1) / dmc)
    
    p = np.where(dmc_loc == 0, 0, ((dmc_loc - bui1) / dmc_loc))
    cc = 0.92 + ((0.0114 * dmc_loc) ** 1.7)
    bui0 = dmc_loc - cc * p
    # Constraints
    
    ###
    # Non-array operations
    ###     
    # bui0 = 0 if (bui0 < 0) else bui0
    
    bui0 = np.where(bui0 < 0, 0, bui0)

    ###
    # Non-array operations
    ###    
    # bui1 = bui0 if (bui1 < dmc) else bui1
    bui1 = np.where(bui1 < dmc_loc, bui0, bui1)
    
    return bui1

def eq27a(dc_loc, dmc_loc):
    # Original equation
    # (0.8 * dc_loc * dmc_loc / (dmc_loc + 0.4 * dc_loc))
    num = 0.8 * dc_loc * dmc_loc
    den = dmc_loc + 0.4 * dc_loc
    bui1 = num/den
    return bui1

def fwi_arr(isi, bui):
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
    if np.any(isi < 0):
        raise ValueError(f'Invalid isi: {isi}')
    if np.any(bui < 0):
        raise ValueError(f'Invalid bui: {bui}')
    # Eqs. 28b, 28a, 29
    ###
    # Non array operations
    ###
    #bb = (0.1 * isi * (1000 / (25 + 108.64 / np.exp(0.023 * bui)))) if (
    #        bui > 80) else (0.1 * isi * (0.626 * (bui ** 0.809) + 2))
    
    bb = np.where(bui > 80, 
                  (0.1 * isi * (1000 / (25 + 108.64 / np.exp(0.023 * bui)))), 
                  (0.1 * isi * (0.626 * (bui ** 0.809) + 2)))
    
    # Eqs. 30b, 30a
    ###
    # Non-array operations
    ###
    # fwi = bb if (bb <= 1.0) else exp(2.72 * ((0.434 * log(bb)) ** 0.647))
    
    fwi = np.where(bb <= 1.0, bb, np.exp(2.72 * ((0.434 * np.log(bb)) ** 0.647)))   
    return fwi

# if __name__ == "__main__":
    
#     ffmc_yda = 8.633961
#     wind = 2.775034475120417*3.6
#     temp = 10.882793908945514
#     rh = 79.0543043718733
#     prec = 0
    
#     ffmc_yda = np.array([ffmc_yda])
#     wind = np.array([wind])
#     temp = np.array([temp])
#     rh = np.array([rh])
#     prec = np.array([prec])
    
#     ffmc = ffmc_arr(ffmc_yda, temp,rh,wind,prec)
#     print("FFMC ", str(ffmc))