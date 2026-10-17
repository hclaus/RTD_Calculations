import json
from pathlib import Path

import pytest

from app.core.reactor import u0, uavg
from app.core.rtd_engine import integrate_rtd
from app.core.solvers import BisectionError, find_e, find_k

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture
def hydraulics():
    legacy = json.loads((FIXTURES / "spreadsheet_known_values.json").read_text())
    ri = legacy["reactor_inputs"]
    return ri["L"], uavg(ri["Q"], ri["a"]), u0(ri["Q"], ri["a"])


def test_find_k_round_trip(hydraulics):
    L, Uavg, U0 = hydraulics
    k_true, Eavg = 1.2, 1.6
    R = integrate_rtd(L, U0, k_true, Eavg)
    k_est = find_k(R, L, U0, Eavg)
    assert k_est == pytest.approx(k_true, rel=1e-6)


def test_find_e_round_trip(hydraulics):
    L, Uavg, U0 = hydraulics
    k, Eavg_true = 0.87, 1.3
    R = integrate_rtd(L, U0, k, Eavg_true)
    eavg_est = find_e(R, L, U0, k)
    assert eavg_est == pytest.approx(Eavg_true, rel=1e-6)


def test_find_k_matches_legacy_within_expected_tolerance():
    """Against the real measured bioassay scenario from the spreadsheet
    (measured C/C0=0.00202, Eavg=1.6), the back-solved k should be close to
    (but not identical to) the legacy VBA's cached estimate, since the
    underlying integration method changed deliberately."""
    legacy = json.loads((FIXTURES / "spreadsheet_known_values.json").read_text())
    ri = legacy["reactor_inputs"]
    L, Uavg, U0 = ri["L"], uavg(ri["Q"], ri["a"]), u0(ri["Q"], ri["a"])
    R = legacy["headline"]["measured_C_C0"]
    Eavg = legacy["headline"]["Eavg"]

    k_est = find_k(R, L, U0, Eavg)
    legacy_k = legacy["headline"]["estimated_k"]
    assert k_est == pytest.approx(legacy_k, rel=0.05)


def test_find_k_raises_bisection_error_for_unreachable_target(hydraulics):
    L, Uavg, U0 = hydraulics
    with pytest.raises(BisectionError):
        find_k(1.5, L, U0, 1.6)  # survival ratio > 1 is physically unreachable


def test_find_e_raises_bisection_error_for_unreachable_target(hydraulics):
    L, Uavg, U0 = hydraulics
    with pytest.raises(BisectionError):
        find_e(1.5, L, U0, 0.87)


def test_find_k_custom_bracket_changes_behavior(hydraulics):
    L, Uavg, U0 = hydraulics
    k_true, Eavg = 10.0, 1.6
    R = integrate_rtd(L, U0, k_true, Eavg)
    # Default bracket (0.001-50) should bracket k_true=10 fine.
    assert find_k(R, L, U0, Eavg) == pytest.approx(k_true, rel=1e-6)
    # A narrow bracket that excludes k_true must fail instead of silently
    # returning a wrong answer.
    with pytest.raises(BisectionError):
        find_k(R, L, U0, Eavg, k_min=0.001, k_max=5.0)
