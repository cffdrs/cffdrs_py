
import numpy as np
from cffdrs.buildup_index_raster import buildup_index


def test_zero_bui_input():
    assert np.array([0, 0]).all() == buildup_index(np.array([0, 0]), np.array([0, 0])).all()