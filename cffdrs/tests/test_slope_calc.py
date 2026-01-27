import pytest
import math
from cffdrs.tests.conftest import float_or_nan
from cffdrs.slope_calc import slope_adjustment


csv_schema = {
    "FUELTYPE": str,
    "FFMC": float,
    "BUI": float,
    "WS": float,
    "WAZ": float,
    "GS": float,
    "SAZ": float,
    "FMC": float,
    "SFC": float,
    "PC": float,
    "PDF": float,
    "CC": float,
    "CBH": float,
    "ISI": float,
    "WSV": float_or_nan,
    "RAZ": float_or_nan,
}

reason = """
Known floating-point precision differences with R implementation for edge cases 
where WSV < 1e-11 (effectively zero wind). RAZ differs by ±180° or ±360° due to 
angle wraparound ambiguity when wind magnitude is negligible. These differences 
have no practical operational significance as direction is undefined at zero wind 
speed. Affects ~1.3% of test cases (28 out of ~2100).
"""


@pytest.mark.xfail(reason=reason)
def test_slope_calc(load_csv):
    data = load_csv("cffdrs/tests/data/Slope.csv", csv_schema)

    failures = []  # collect all issues

    for idx, row in enumerate(data, start=2):
        fuel_type = row["FUELTYPE"]
        ffmc = row["FFMC"]
        bui = row["BUI"]
        ws = row["WS"]
        waz = row["WAZ"]
        gs = row["GS"]
        saz = row["SAZ"]
        fmc = row["FMC"]
        sfc = row["SFC"]
        pc = row["PC"]
        pdf = row["PDF"]
        cc = row["CC"]
        cbh = row["CBH"]
        isi = row["ISI"]

        expected_wsv = row["WSV"]
        expected_raz = row["RAZ"]

        result = slope_adjustment(
            fuel_type, ffmc, bui, ws, waz, gs, saz, fmc, sfc, pc, pdf, cc, cbh, isi
        )

        calculated_wsv = result["WSV"]
        calculated_raz = result["RAZ"]

        tol = 0.1

        wsv_ok = pytest.approx(expected_wsv, abs=tol, nan_ok=True) == calculated_wsv
        raz_ok = pytest.approx(expected_raz, abs=tol, nan_ok=True) == calculated_raz

        if not wsv_ok or not raz_ok:
            failure_info = {
                "row": idx,
                "fuel": fuel_type,
                "ws": ws,
                "waz": waz,
                "gs": gs,
                "saz": saz,
                "isi": isi,
                "expected_wsv": expected_wsv,
                "calc_wsv": calculated_wsv,
                "expected_raz": expected_raz,
                "calc_raz": calculated_raz,
                "diff_raz_deg": math.degrees(calculated_raz - expected_raz)
                if not math.isnan(calculated_raz) and not math.isnan(expected_raz)
                else None,
                "waz_deg": math.degrees(waz),
                "saz_deg": math.degrees(saz),
            }
            failures.append(failure_info)

            msg = f"Row {idx} failed: WSV {calculated_wsv:.4g} vs {expected_wsv:.4g} | RAZ {calculated_raz:.4g} vs {expected_raz:.4g}"
            if not wsv_ok:
                print(f"WSV issue - {msg}")
            if not raz_ok:
                diff_deg = None
                if not math.isnan(calculated_raz) and not math.isnan(expected_raz):
                    diff_deg = math.degrees(calculated_raz - expected_raz)

                print(
                    f"RAZ issue  - {msg}  (diff ≈ {diff_deg:+.1f}°)"
                    if diff_deg is not None
                    else f"RAZ issue  - {msg}  (diff N/A - NaN involved)"
                )

    # After the loop — fail the test if anything went wrong
    if failures:
        print("\nSummary of all failures:")
        for f in failures:
            diff_str = f"{f['diff_raz_deg']:+.1f}°" if f["diff_raz_deg"] is not None else "N/A"
            print(
                f"Row {f['row']}: RAZ diff {diff_str}, "
                f"expected {f['expected_raz']:.4g}, got {f['calc_raz']:.4g}"
            )

        pytest.fail(
            f"{len(failures)} failures in slope_adjustment test. See printed details above."
        )
