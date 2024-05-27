
import pytest
import numpy as np
from cffdrs.buildup_index_raster import buildup_index

@pytest.mark.parametrize(
    "dc,dmc",
    [
        (np.array([0, 0]), np.array([0, 0])),
        (np.array([0, 100]), np.array([0, 100])),
        (np.array([-1, 0]), np.array([-1, 0])),
    ],
)
def test_bui_zero_cases(dc, dmc):
    assert np.array([0, 0]).all() == buildup_index(dc, dmc).all()