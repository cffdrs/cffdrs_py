
import numpy as np
from cffdrs.buildup_index_raster import buildup_index


def test_zero_bui_input():
    buildup_index(np.array(0, 0), np.array(0, 0))