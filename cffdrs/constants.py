# Constants for cffdrs_py

import math
from typing import Literal, get_args


FFMC_COEFFICIENT = 250.0 * 59.5 / 101.0

FuelType = Literal[
    "C1",
    "C2",
    "C3",
    "C4",
    "C5",
    "C6",
    "C7",
    "D1",
    "M1",
    "M2",
    "M3",
    "M4",
    "S1",
    "S2",
    "S3",
    "O1A",
    "O1B",
    "NF",
    "WA",
]

# Fuel type mappings
# BUIo = average BUI for the fuel type
# q = proportion of maximum possible spread rate that is reached at a standard BUI
# Table 7 (FCFDG 1992)
FUEL_TYPE_DEFAULTS = {
    "C1": {"BUIo": 72, "Q": 0.9},
    "C2": {"BUIo": 64, "Q": 0.7},
    "C3": {"BUIo": 62, "Q": 0.75},
    "C4": {"BUIo": 66, "Q": 0.8},
    "C5": {"BUIo": 56, "Q": 0.8},
    "C6": {"BUIo": 62, "Q": 0.8},
    "C7": {"BUIo": 106, "Q": 0.85},
    "D1": {"BUIo": 32, "Q": 0.9},
    "M1": {"BUIo": 50, "Q": 0.8},
    "M2": {"BUIo": 50, "Q": 0.8},
    "M3": {"BUIo": 50, "Q": 0.8},
    "M4": {"BUIo": 50, "Q": 0.8},
    "S1": {"BUIo": 38, "Q": 0.75},
    "S2": {"BUIo": 63, "Q": 0.75},
    "S3": {"BUIo": 31, "Q": 0.75},
    "O1A": {"BUIo": 1, "Q": 1.0},
    "O1B": {"BUIo": 1, "Q": 1.0},
}

# fuel type rate of spread regression coefficients
# Table 6 (FCFDG 1992)
FUEL_TYPE_ROS = {
    "C1": {"a": 90, "b": 0.0649, "c0": 4.5},
    "C2": {"a": 110, "b": 0.0282, "c0": 1.5},
    "C3": {"a": 110, "b": 0.0444, "c0": 3.0},
    "C4": {"a": 110, "b": 0.0293, "c0": 1.5},
    "C5": {"a": 30, "b": 0.0697, "c0": 4.0},
    "C6": {"a": 30, "b": 0.0800, "c0": 3.0},
    "C7": {"a": 45, "b": 0.0305, "c0": 2.0},
    "D1": {"a": 30, "b": 0.0232, "c0": 1.6},
    "M1": {"a": 0, "b": 0, "c0": 0},
    "M2": {"a": 0, "b": 0, "c0": 0},
    "M3": {"a": 120, "b": 0.0572, "c0": 1.4},
    "M4": {"a": 100, "b": 0.0404, "c0": 1.48},
    "S1": {"a": 75, "b": 0.0297, "c0": 1.3},
    "S2": {"a": 40, "b": 0.0438, "c0": 1.7},
    "S3": {"a": 55, "b": 0.0829, "c0": 3.2},
    "O1A": {"a": 190, "b": 0.0310, "c0": 1.4},
    "O1B": {"a": 250, "b": 0.0350, "c0": 1.7},
}

# --- Vectorization support -------------------------------------------------
#
# The functions above dispatch on fuel_type strings, and a few (rate_of_spread
# for M1-M4, slope_adjustment) recurse into themselves with a different fixed
# fuel_type string to compute a sub-formula. Neither pattern can be compiled
# or traced by array-vectorization tools like numba or jax. Every module in
# this package that branches on fuel type also exposes a "_core" sibling
# function that takes an int fuel_type_code (looked up here) instead of a
# string, with any such recursive calls unrolled into flat arithmetic. This
# package does not depend on numba or jax itself - the *_core functions are
# just plain, vectorization-friendly Python that a caller can wrap themselves.
#
# FUEL_TYPE_NAMES[code] == name, and FUEL_TYPE_CODES[name] == code.
FUEL_TYPE_NAMES: tuple = get_args(FuelType)
FUEL_TYPE_CODES: dict = {name: code for code, name in enumerate(FUEL_TYPE_NAMES)}

(
    C1,
    C2,
    C3,
    C4,
    C5,
    C6,
    C7,
    D1,
    M1,
    M2,
    M3,
    M4,
    S1,
    S2,
    S3,
    O1A,
    O1B,
    NF,
    WA,
) = range(len(FUEL_TYPE_NAMES))

# ROS_A/ROS_B/ROS_C0[code] mirror FUEL_TYPE_ROS[name]["a"/"b"/"c0"], math.nan
# where a fuel type has no entry (NF, WA).
ROS_A = tuple(FUEL_TYPE_ROS.get(name, {}).get("a", math.nan) for name in FUEL_TYPE_NAMES)
ROS_B = tuple(FUEL_TYPE_ROS.get(name, {}).get("b", math.nan) for name in FUEL_TYPE_NAMES)
ROS_C0 = tuple(FUEL_TYPE_ROS.get(name, {}).get("c0", math.nan) for name in FUEL_TYPE_NAMES)

# BUI_O/BUI_Q[code] mirror FUEL_TYPE_DEFAULTS[name]["BUIo"/"Q"], math.nan
# where a fuel type has no entry (NF, WA).
BUI_O = tuple(FUEL_TYPE_DEFAULTS.get(name, {}).get("BUIo", math.nan) for name in FUEL_TYPE_NAMES)
BUI_Q = tuple(FUEL_TYPE_DEFAULTS.get(name, {}).get("Q", math.nan) for name in FUEL_TYPE_NAMES)
