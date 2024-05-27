import numpy as np
import numpy.typing as npt

def buildup_index(dc: npt.NDArray[np.float64], dmc: npt.NDArray[np.float64]):
    """
    Buildup Index Raster Calculation

    Parameters
    ----------
    dc : numpy array of float
       The Drought Code raster
    dmc : numpy array of float
       The Duff Moisture Code raster

    Returns
    -------
    numpy array of float
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
    np.array(1)
        
    pass