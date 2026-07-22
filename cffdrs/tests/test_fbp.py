import pytest
from cffdrs.fbp import fbp
from cffdrs.models import FBPInput
from cffdrs.constants import FUEL_TYPE_CODES
from cffdrs.fire_behaviour_prediction import (
    fire_behaviour_prediction,
    _fire_behaviour_prediction,
    FD_SURFACE,
    FD_INTERMITTENT,
    FD_CROWN,
    FD_NONE,
)
from cffdrs.tests.conftest import run_csv_test, string_or_none


primary_output_csv_schema = {
    "ID": int,
    "CFB": float,
    "CFC": float,
    "FD": string_or_none,
    "HFI": float,
    "RAZ": float,
    "ROS": float,
    "SFC": float,
    "TFC": float,
}


secondary_output_csv_schema = {
    "ID": int,
    "BE": float,
    "SF": float,
    "ISI": float,
    "FFMC": float,
    "FMC": float,
    "D0": float,
    "RSO": float,
    "CSI": float,
    "FROS": float,
    "BROS": float,
    "HROSt": float,
    "FROSt": float,
    "BROSt": float,
    "FCFB": float,
    "BCFB": float,
    "FFI": float,
    "BFI": float,
    "FTFC": float,
    "BTFC": float,
    "TI": float,
    "FTI": float,
    "BTI": float,
    "LB": float,
    "LBt": float,
    "WSV": float,
    "DH": float,
    "DB": float,
    "DF": float,
    "TROS": float,
    "TROSt": float,
    "TCFB": float,
    "TFI": float,
    "TTFC": float,
    "TTI": float,
}

all_output_csv_schema = {**primary_output_csv_schema, **secondary_output_csv_schema}


def test_fbp_01(fbp_input_data, load_csv):
    """Test FBP function with default primary option"""

    # should default to primary output
    results = fbp(input=fbp_input_data)

    # convert results to dicts
    results_dicts = [res.__dict__ for res in results]

    expected_primary = load_csv(
        "cffdrs/tests/data/fbp_01.csv",
        primary_output_csv_schema,
    )

    run_csv_test(results_dicts, expected_primary)


def test_fbp_02(fbp_input_data, load_csv):
    """Test FBP function with declared primary option"""

    results = fbp(input=fbp_input_data, output="Primary")

    # convert results to dicts
    results_dicts = [res.__dict__ for res in results]

    expected_primary = load_csv(
        "cffdrs/tests/data/fbp_01.csv",  # fbp_02.csv not included as identical to fbp_01.csv
        primary_output_csv_schema,
    )

    run_csv_test(results_dicts, expected_primary)


def test_fbp_03(fbp_input_data, load_csv):
    """Test FBP function with shorthand declared primary option"""

    results = fbp(input=fbp_input_data, output="P")

    # convert results to dicts
    results_dicts = [res.__dict__ for res in results]

    expected_primary = load_csv(
        "cffdrs/tests/data/fbp_01.csv",  # fbp_03.csv not included as identical to fbp_01.csv
        primary_output_csv_schema,
    )

    for result_dict, expected in zip(results_dicts, expected_primary):
        for key in expected.keys():
            assert pytest.approx(result_dict[key.lower()], abs=0.01) == expected[key]


def test_fbp_04(fbp_input_data, load_csv):
    """Test FBP function with declared secondary option"""

    results = fbp(input=fbp_input_data, output="Secondary")

    # convert results to dicts
    results_dicts = [res.__dict__ for res in results]

    expected_secondary = load_csv(
        "cffdrs/tests/data/fbp_04.csv",
        secondary_output_csv_schema,
    )

    run_csv_test(results_dicts, expected_secondary)


def test_fbp_05(fbp_input_data, load_csv):
    """Test FBP function with declared shorthand secondary option"""

    results = fbp(input=fbp_input_data, output="S")

    # convert results to dicts
    results_dicts = [res.__dict__ for res in results]

    expected_secondary = load_csv(
        "cffdrs/tests/data/fbp_04.csv",
        secondary_output_csv_schema,
    )

    run_csv_test(results_dicts, expected_secondary)


def test_fbp_06(fbp_input_data, load_csv):
    """Test FBP function with declared All option"""

    results = fbp(input=fbp_input_data, output="All")

    # convert results to dicts
    results_dicts = [res.__dict__ for res in results]

    expected = load_csv(
        "cffdrs/tests/data/fbp_06.csv",
        all_output_csv_schema,
    )

    run_csv_test(results_dicts, expected)


def test_fbp_07(fbp_input_data, load_csv):
    """Test FBP function with declared shorthand All option"""

    results = fbp(input=fbp_input_data, output="A")

    # convert results to dicts
    results_dicts = [res.__dict__ for res in results]

    expected = load_csv(
        "cffdrs/tests/data/fbp_06.csv",
        all_output_csv_schema,
    )

    run_csv_test(results_dicts, expected)


############
# fbp_08 and fbp_09 not included as the same data exists in other csv tests
############


@pytest.mark.xfail(
    reason="I believe default CFL should be 0, not 1 as it's currently written - https://github.com/cffdrs/cffdrs_r/blob/4e40dd3af841f3a708abd83d7f2d43fefb08649a/R/fire_behaviour_prediction.r#L38"
)
def test_fbp_10(load_csv):
    """Test FBP function with no input (uses default FBPInput)"""

    results = fbp()

    # convert results to dicts
    results_dicts = [res.__dict__ for res in results]

    expected = load_csv(
        "cffdrs/tests/data/fbp_10.csv",
        all_output_csv_schema,
    )

    run_csv_test(results_dicts, expected)


def test_fbp_11(fbp_input_data, load_csv):
    """Test FBP function with all NF fuel type"""

    modified_inputs = [FBPInput(**{**row.__dict__, "fuel_type": "NF"}) for row in fbp_input_data]
    results = fbp(input=modified_inputs, output="All")

    # convert results to dicts
    results_dicts = [res.__dict__ for res in results]

    expected = load_csv(
        "cffdrs/tests/data/fbp_11.csv",
        all_output_csv_schema,
    )

    run_csv_test(results_dicts, expected)


def test_fbp_12(fbp_input_data, load_csv):
    """Test FBP function with all WA fuel type"""

    modified_inputs = [FBPInput(**{**row.__dict__, "fuel_type": "WA"}) for row in fbp_input_data]
    results = fbp(input=modified_inputs, output="All")

    # convert results to dicts
    results_dicts = [res.__dict__ for res in results]

    expected = load_csv(
        "cffdrs/tests/data/fbp_12.csv",
        all_output_csv_schema,
    )

    run_csv_test(results_dicts, expected)


_FD_CODE_BY_STRING = {"S": FD_SURFACE, "I": FD_INTERMITTENT, "C": FD_CROWN, None: FD_NONE}


def _assert_core_matches_scalar(fbp_input: FBPInput):
    scalar = fire_behaviour_prediction(fbp_input, "All")
    core = _fire_behaviour_prediction(
        FUEL_TYPE_CODES[fbp_input.fuel_type],
        fbp_input.ffmc,
        fbp_input.bui,
        fbp_input.ws,
        fbp_input.wd,
        fbp_input.gs,
        fbp_input.aspect,
        fbp_input.pc,
        fbp_input.pdf,
        fbp_input.cc,
        fbp_input.gfl,
        fbp_input.cbh,
        fbp_input.cfl,
        fbp_input.fmc,
        fbp_input.isi,
        fbp_input.lat,
        fbp_input.lon,
        fbp_input.elv,
        fbp_input.dj,
        fbp_input.d0,
        fbp_input.sd,
        fbp_input.sh,
        fbp_input.hr,
        fbp_input.theta,
        fbp_input.accel,
        fbp_input.bui_eff,
    )

    for name in scalar.__dataclass_fields__:
        if name == "id":
            continue
        scalar_value = getattr(scalar, name)
        if name == "fd":
            assert _FD_CODE_BY_STRING[scalar_value] == core.fd_code, (
                f"fd mismatch for {fbp_input}: scalar={scalar_value}, core fd_code={core.fd_code}"
            )
            continue
        core_value = getattr(core, name)
        assert pytest.approx(scalar_value, abs=1e-6, nan_ok=True) == core_value, (
            f"{name} mismatch for {fbp_input}: scalar={scalar_value}, core={core_value}"
        )


def test_fbp_core_equivalence(fbp_input_data):
    """
    _fire_behaviour_prediction (int fuel_type_code, pre-validated scalar
    inputs, no recursion) must match fire_behaviour_prediction(..., "All")
    field-for-field, for every row of test_fbp.csv.
    """
    for fbp_input in fbp_input_data:
        _assert_core_matches_scalar(fbp_input)


def test_fbp_core_equivalence_nf(fbp_input_data):
    """Same equivalence check, forcing every row to the NF (non-fuel) fuel type."""
    for row in fbp_input_data:
        _assert_core_matches_scalar(FBPInput(**{**row.__dict__, "fuel_type": "NF"}))


def test_fbp_core_equivalence_wa(fbp_input_data):
    """Same equivalence check, forcing every row to the WA (water) fuel type."""
    for row in fbp_input_data:
        _assert_core_matches_scalar(FBPInput(**{**row.__dict__, "fuel_type": "WA"}))
