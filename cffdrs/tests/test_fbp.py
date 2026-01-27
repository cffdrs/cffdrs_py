import pytest
from cffdrs.fbp import fbp
from cffdrs.models import FBPInput
from cffdrs.tests.conftest import run_csv_test


primary_output_csv_schema = {
    "ID": int,
    "CFB": float,
    "CFC": float,
    "FD": str,
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


@pytest.mark.xfail(reason="I believe default CFL should be 0, not 1 as it's currently written")
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
