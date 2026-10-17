import json
import math
from pathlib import Path

import pytest

from app.core.reactor import mean_residence_time, u0, uavg
from app.core.sweep import build_sweep_table, fit_forced_exponential

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture
def legacy():
    return json.loads((FIXTURES / "spreadsheet_known_values.json").read_text())


def test_build_sweep_table_full_range_not_truncated(legacy):
    """Deliberate fix vs. the legacy chart (which only plotted 5 of 12 rows):
    the table itself must return every requested row."""
    ri = legacy["reactor_inputs"]
    L, Uavg, U0 = ri["L"], uavg(ri["Q"], ri["a"]), u0(ri["Q"], ri["a"])
    q = mean_residence_time(L, Uavg)
    eavg_values = [row["Eavg"] for row in legacy["sweep"]["rows"]]

    rows = build_sweep_table(L, U0, q, legacy["sweep"]["k"], eavg_values)
    assert len(rows) == 12
    assert [r["Eavg"] for r in rows] == eavg_values


def test_build_sweep_table_davg_matches_legacy_exactly(legacy):
    """Davg = Eavg*q is pure arithmetic (no integration involved) so it
    should match the legacy cached values exactly, unlike C_C0_rtd."""
    ri = legacy["reactor_inputs"]
    L, Uavg, U0 = ri["L"], uavg(ri["Q"], ri["a"]), u0(ri["Q"], ri["a"])
    q = mean_residence_time(L, Uavg)
    eavg_values = [row["Eavg"] for row in legacy["sweep"]["rows"]]

    rows = build_sweep_table(L, U0, q, legacy["sweep"]["k"], eavg_values)
    for row, legacy_row in zip(rows, legacy["sweep"]["rows"]):
        assert row["Davg"] == pytest.approx(legacy_row["Davg"], rel=1e-9)
        assert row["C_C0_well_mixed"] == pytest.approx(legacy_row["C_C0_well_mixed"], rel=1e-9)


def test_build_sweep_table_rtd_survival_deviates_within_expected_bound(legacy):
    ri = legacy["reactor_inputs"]
    L, Uavg, U0 = ri["L"], uavg(ri["Q"], ri["a"]), u0(ri["Q"], ri["a"])
    q = mean_residence_time(L, Uavg)
    eavg_values = [row["Eavg"] for row in legacy["sweep"]["rows"]]

    rows = build_sweep_table(L, U0, q, legacy["sweep"]["k"], eavg_values)
    for row, legacy_row in zip(rows, legacy["sweep"]["rows"]):
        rel_diff = (legacy_row["C_C0_rtd"] - row["C_C0_rtd"]) / legacy_row["C_C0_rtd"]
        assert 0.0 < rel_diff < 0.15


def test_fit_forced_exponential_recovers_known_slope():
    b_true = -0.75
    x = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5]
    # skip x=0 (ln(1)=0 contributes nothing but avoid log(0) issues elsewhere)
    y = [math.exp(b_true * xi) for xi in x]
    b_fit = fit_forced_exponential(x, y)
    assert b_fit == pytest.approx(b_true, rel=1e-9)


def test_log10_inactivation_consistent_with_survival(legacy):
    ri = legacy["reactor_inputs"]
    L, Uavg, U0 = ri["L"], uavg(ri["Q"], ri["a"]), u0(ri["Q"], ri["a"])
    q = mean_residence_time(L, Uavg)
    rows = build_sweep_table(L, U0, q, 1.0, [0.5, 1.0, 2.0])
    for row in rows:
        assert row["log10_inactivation"] == pytest.approx(-math.log10(row["C_C0_rtd"]), rel=1e-12)
