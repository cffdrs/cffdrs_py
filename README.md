# cffdrs

This is a python translation of the R code available from [cran](https://cran.r-project.org/web/packages/cffdrs/index.html)

This project provides a group of new functions to calculate the outputs of the two main components of the [Canadian Forest Fire Danger Rating System (CFFDRS) Van Wagner and Pickett (1985)](https://cfs.nrcan.gc.ca/publications?id=19973) at various time scales: the [Fire Weather Index (FWI) System Wan Wagner (1985)](https://cfs.nrcan.gc.ca/publications?id=19927) and the [Fire Behaviour Prediction (FBP) System Forestry Canada Fire Danger Group (1992)](http://cfs.nrcan.gc.ca/pubwarehouse/pdfs/10068.pdf). Some functions have two versions, table and raster based.

## Usage

### Fire Weather Index (FWI) System

The FWI system consists of several components that can be calculated individually or chained together.

#### Fine Fuel Moisture Code (FFMC)

```python
from cffdrs.fwi import fine_fuel_moisture_code

# Calculate FFMC
ffmc = fine_fuel_moisture_code(ffmc_yda=85.0, temp=20.0, rh=60.0, ws=15.0, prec=0.0)
```

#### Duff Moisture Code (DMC)

```python
from cffdrs.fwi import duff_moisture_code

# Calculate DMC
dmc = duff_moisture_code(dmc_yda=20.0, temp=25.0, rh=45.0, prec=0.0, lat=55.0, mon=6)
```

#### Drought Code (DC)

```python
from cffdrs.fwi import drought_code

# Calculate DC
dc = drought_code(dc_yda=150.0, temp=30.0, rh=30.0, prec=0.0, lat=50.0, mon=7)
```

#### Initial Spread Index (ISI)

```python
from cffdrs.fwi import initial_spread_index

# Calculate ISI
isi = initial_spread_index(ffmc=85.0, ws=20.0)
```

#### Buildup Index (BUI)

```python
from cffdrs.fwi import buildup_index

# Calculate BUI
bui = buildup_index(dmc=25.0, dc=160.0)
```

#### Fire Weather Index (FWI)

```python
from cffdrs.fwi import fire_weather_index

# Calculate FWI
fwi = fire_weather_index(isi=10.0, bui=50.0)
```

### Fire Behaviour Prediction (FBP) System

The FBP system predicts fire behavior based on fuel type and weather conditions.

```python
from cffdrs.fbp import fbp
from cffdrs.models import FBPInput

# Create input parameters
input_data = FBPInput(
    fuel_type="C2",
    ffmc=85.0,
    bui=40.0,
    ws=15.0,
    wd=45.0,
    gs=0.0,
    lat=55.0,
    lon=-120.0,
    elv=100.0,
    dj=180.0,
    d0=0.0,
    hr=1.0,
    pc=50.0,
    pdf=35.0,
    cc=80.0,
    gfl=0.35,
    cbh=2.0,
    cfl=1.0,
    isi=8.0,
    fmc=0.0,
    theta=0.0,
    accel=0,
    aspect=0.0,
    bui_eff=1
)

# Calculate primary outputs
primary_results = fbp(input=input_data, output="Primary")


# Calculate all outputs
all_results = fbp(input=input_data, output="All")
print(all_results[0])
# Output: FBPAllOutput(id='1', cfb=..., cfc=..., fd='S', hfi=..., raz=..., ros=..., sfc=..., tfc=..., be=..., sf=..., isi=..., ffmc=..., fmc=..., d0=..., rso=..., csi=..., fros=..., bros=..., hrost=..., frost=..., brost=..., fcfb=..., bcfb=..., ffi=..., bfi=..., ftfc=..., btfc=..., ti=..., fti=..., bti=..., lb=..., lbt=..., wsv=..., dh=..., db=..., df=..., tros=..., trost=..., tcfb=..., tfi=..., ttfc=..., tti=...)
```

#### Batch Processing

You can process multiple inputs at once:

```python
inputs = [FBPInput(fuel_type="C2", ffmc=85.0, bui=40.0, ws=15.0), FBPInput(fuel_type="C3", ffmc=90.0, bui=50.0, ws=20.0)]
results = fbp(input=inputs, output="Primary")
for result in results:
    print(result)
```

#### Vectorization support

`cffdrs` does not vectorize anything itself and has no runtime dependencies, so it does not depend
on numpy, numba, or jax. This section is for callers who want to do that themselves: `fbp()`'s
Python-level loop over `FBPInput` objects can be a bottleneck over large batches, and two things
about the normal API get in the way of vectorizing it directly. Fuel type is a string (`"C2"`),
and the M1-M4 mixedwood fuel types compute part of their rate of spread by calling back into
`rate_of_spread`/`slope_adjustment` with a *different* fuel type, a data-dependent recursive call
that array-vectorization tools like `numba` or `jax` can't trace or compile.

To make that possible, every FBP function that branches on fuel type also has a leading-underscore
counterpart that a caller can wrap themselves:

- it takes an int `fuel_type_code` instead of a fuel type string (see
  `cffdrs.constants.FUEL_TYPE_CODES` for the mapping)
- the M1-M4 recursive sub-calls are unrolled into flat arithmetic instead
- it has no Python object construction or input validation in the hot path (plain scalars in,
  plain scalars/tuples out)

These functions are prefixed with `_` because they're a lower-level, less stable surface than the rest of the public API. They're still importable directly by name. You bring the vectorization tool (`numpy.vectorize`, `numba.guvectorize`, etc.).

The one you'll usually want is `_fire_behaviour_prediction` in `cffdrs.fire_behaviour_prediction`.
It's what `fire_behaviour_prediction(input, "All")` delegates to internally, so it mirrors that
computation exactly, but takes already-validated scalar inputs (the range-clamping and unit
conversion `FBPInput.__post_init__` normally does; a caller vectorizing over arrays of inputs is
expected to do the equivalent once, up front, over the whole array) and an int `fuel_type_code`,
and always returns the full set of fields as a `NamedTuple` (`_FBPOutput`, with `fd_code` in place
of the `fd` string; see `FD_SURFACE`/`FD_INTERMITTENT`/`FD_CROWN`/`FD_NONE`).

```python
import numpy as np
from cffdrs.constants import FUEL_TYPE_CODES
from cffdrs.fire_behaviour_prediction import _fire_behaviour_prediction

# Map string fuel types to int codes once, up front.
fuel_type_codes = np.array([FUEL_TYPE_CODES[ft] for ft in ["C2", "C3", "D1"]])
ffmc = np.array([85.0, 90.0, 88.0])
bui = np.array([40.0, 50.0, 45.0])
ws = np.array([15.0, 20.0, 10.0])
# ... one array per _fire_behaviour_prediction parameter, pre-clamped/validated
# the same way FBPInput.__post_init__ does for the scalar API.

vectorized = np.vectorize(_fire_behaviour_prediction, otypes=[object])
results = vectorized(fuel_type_codes, ffmc, bui, ws, ...)  # one _FBPOutput per row
ros = np.array([r.ros for r in results])
```

A `numba.guvectorize`'d version of the same call would compile down to native code and avoid the
Python-level loop `numpy.vectorize` still does under the hood. That's useful once you're processing
rasters or large tables, and it's the lack of recursion/string dispatch in `_fire_behaviour_prediction`
that makes it eligible for `nopython` mode in the first place.
