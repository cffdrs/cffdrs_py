import cProfile
import os
import pstats
from osgeo import gdal
from cffdrs.fwi import bui

def get_raster_array(path: str):
    source = gdal.Open(path, gdal.GA_ReadOnly)
    source_band = source.GetRasterBand(1)
    nodata_value = source_band.GetNoDataValue()
    source_data = source.ReadAsArray()
    source_data[source_data == nodata_value] = 0
    del source
    return source_data

def profile_numba_2000m_tif(bui_tif_dir):
    """
    Run profiler against the 2000m resolution tifs using the numpy impl and capture profiles
    """
    profile_dir = os.path.join(os.path.dirname(__file__), 'profiles')
    dc_array = get_raster_array(os.path.join(bui_tif_dir, 'dc20240528.tif'))
    dmc_array = get_raster_array(os.path.join(bui_tif_dir, 'dmc20240528.tif'))

    for i in range(10):
        with cProfile.Profile() as pr:
            bui(dmc_array, dc_array)
            with open(f"{profile_dir}/numba_2000m_profile_iter_{i}.txt", 'w') as f:
                pstats.Stats( pr, stream=f ).strip_dirs().sort_stats("cumtime").print_stats()

if __name__ == "__main__":
    """
    Reads in 2000m dc and dmc rasters to profile the BUI calculation using numba's 
    vectorized annotation.
    """
    parent_dir = os.path.dirname(__file__)
    bui_tif_dir = os.path.join(parent_dir, 'fixtures')
    profile_numba_2000m_tif(bui_tif_dir)