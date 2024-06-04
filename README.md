# cffdrs

This is a python translation of the R code available from [cran](https://cran.r-project.org/web/packages/cffdrs/index.html) - this does not include all the functionality available in that package yet.

This project provides a group of new functions to calculate the outputs of the two main components of the [Canadian Forest Fire Danger Rating System (CFFDRS) Van Wagner and Pickett (1985)](https://cfs.nrcan.gc.ca/publications?id=19973) at various time scales: the [Fire Weather Index (FWI) System Wan Wagner (1985)](https://cfs.nrcan.gc.ca/publications?id=19927) and the [Fire Behaviour Prediction (FBP) System Forestry Canada Fire Danger Group (1992)](http://cfs.nrcan.gc.ca/pubwarehouse/pdfs/10068.pdf). Some functions have two versions, table and raster based.

## Local Development Setup

1. Install poetry for `venv` and dependency management: https://python-poetry.org/docs/#installation
2. Verify your Python version is version 3.10 or greater with `python --version`
3. Activate the poetry shell with `poetry shell`; this will create a `venv` if not run before, or else activate an existing `venv`.
4. Install dependencies with `poetry install`; at the current time of writing the only runtime dependency is `numba` used for vectorizing scalar FWI calculations, and the only dev dependency is `pytest` for running tests.
5. Verify tests run by running `pytest` in the root `cffdrs` directory.
