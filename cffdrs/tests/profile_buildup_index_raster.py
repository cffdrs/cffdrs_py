import cProfile
import os
import pstats
from osgeo import gdal
from cffdrs.raster.buildup_index_raster import buildup_index, buildup_index_vectorized

def get_raster_array(path: str):
    source = gdal.Open(path, gdal.GA_ReadOnly)
    source_band = source.GetRasterBand(1)
    nodata_value = source_band.GetNoDataValue()
    source_data = source.ReadAsArray()
    source_data[source_data == nodata_value] = 0
    del source
    return source_data


def profile_2000m_tif_old(bui_tif_dir):
    """
    Run profiler against the 2000m resolution tiffs using the numpy impl and capture profiles
    """
    profile_dir = os.path.join(os.path.dirname(__file__), 'profiles')
    for i in range(10):
        dc_array = get_raster_array(os.path.join(bui_tif_dir, 'dc20240528.tif'))
        dmc_array = get_raster_array(os.path.join(bui_tif_dir, 'dmc20240528.tif'))
        cProfile.runctx('buildup_index(dc_array, dmc_array)', 
                        {'buildup_index': buildup_index, 'dc_array': dc_array, 'dmc_array': dmc_array }, 
                        { }, 
                        f"{profile_dir}/2000m_profile_iter_{i}.prof")


def profile_25m_tif(bui_tif_dir):
    """
    Run profiler against the 2000m resolution tifs using the numpy impl and capture profiles
    """
    profile_dir = os.path.join(os.path.dirname(__file__), 'profiles')
    dc_array = get_raster_array(os.path.join(bui_tif_dir, 'dc20240528_25m.tif'))
    dmc_array = get_raster_array(os.path.join(bui_tif_dir, 'dmc20240528_25m.tif'))

    for i in range(10):
        with cProfile.Profile() as pr:
            buildup_index(dc_array, dmc_array)
            with open(f"{profile_dir}/25m_profile_iter_{i}.txt", 'w') as f:
                pstats.Stats( pr, stream=f ).strip_dirs().sort_stats("cumtime").print_stats()

def profile_vectorized_25m_tif(bui_tif_dir):
    """
    Run profiler against the 2000m resolution tifs using the numpy impl and capture profiles
    """
    profile_dir = os.path.join(os.path.dirname(__file__), 'profiles')
    dc_array = get_raster_array(os.path.join(bui_tif_dir, 'dc20240528.tif'))
    dmc_array = get_raster_array(os.path.join(bui_tif_dir, 'dmc20240528.tif'))

    for i in range(10):
        with cProfile.Profile() as pr:
            buildup_index_vectorized(dmc_array, dc_array)
            with open(f"{profile_dir}/vec_2000m_profile_iter_{i}.txt", 'w') as f:
                pstats.Stats( pr, stream=f ).strip_dirs().sort_stats("cumtime").print_stats()

def profile_numba_25m_tif(bui_tif_dir):
    """
    Run profiler against the 2000m resolution tifs using the numpy impl and capture profiles
    """
    profile_dir = os.path.join(os.path.dirname(__file__), 'profiles')
    dc_array = get_raster_array(os.path.join(bui_tif_dir, 'dc20240528.tif'))
    dmc_array = get_raster_array(os.path.join(bui_tif_dir, 'dmc20240528.tif'))

    for i in range(10):
        with cProfile.Profile() as pr:
            buildup_index_vectorized(dmc_array, dc_array)
            with open(f"{profile_dir}/numba_2000m_profile_iter_{i}.txt", 'w') as f:
                pstats.Stats( pr, stream=f ).strip_dirs().sort_stats("cumtime").print_stats()

if __name__ == "__main__":
    parent_dir = os.path.dirname(__file__)
    bui_tif_dir = os.path.join(parent_dir, 'fixtures')
    profile_numba_25m_tif(bui_tif_dir)