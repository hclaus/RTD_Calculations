import json
from pathlib import Path

import numpy as np
import pytest

from app.core.rtd_engine import integrate_rtd, rtd_kernel

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture
def legacy():
    return json.loads((FIXTURES / "spreadsheet_known_values.json").read_text())


def _hydraulics(reactor_inputs):
    from app.core.reactor import u0, uavg

    L, a, Q = reactor_inputs["L"], reactor_inputs["a"], reactor_inputs["Q"]
    return L, uavg(Q, a), u0(Q, a)


def test_integrate_rtd_against_independent_fine_trapezoid(legacy):
    """Cross-check scipy.quad against a very fine (1,000,000-point), correctly
    weighted trapezoidal rule -- i.e. an independent numerical method, not
    the legacy VBA's buggy one -- to validate the adaptive quadrature itself."""
    L, Uavg, U0 = _hydraulics(legacy["reactor_inputs"])
    k, Eavg = legacy["headline"]["k"], legacy["headline"]["Eavg"]

    t_min, t_max = L / U0, 100.0 * L / U0
    t = np.linspace(t_min, t_max, 1_000_001)
    y = np.exp(-k * Eavg * t) * L / U0 / t**2
    reference = float(np.trapezoid(y, t))

    result = integrate_rtd(L, U0, k, Eavg)
    assert result == pytest.approx(reference, rel=1e-5)


def test_integrate_rtd_deviates_from_legacy_by_the_documented_bug_magnitude(legacy):
    """The legacy VBA Integrate_RTD has a known off-by-one that overestimates
    the integral. This test documents that the new (correct) scipy-based
    result differs from the legacy cached value by a small, expected amount --
    so nobody mistakes a ~4% gap here for a porting error."""
    L, Uavg, U0 = _hydraulics(legacy["reactor_inputs"])
    k, Eavg = legacy["headline"]["k"], legacy["headline"]["Eavg"]
    legacy_value = legacy["headline"]["C_C0_reactor"]

    result = integrate_rtd(L, U0, k, Eavg)
    rel_diff = (legacy_value - result) / legacy_value

    assert result < legacy_value  # new engine removes an overestimate
    assert 0.01 < rel_diff < 0.10  # observed ~4% for the headline scenario


def test_rtd_kernel_decreases_with_t(legacy):
    L, Uavg, U0 = _hydraulics(legacy["reactor_inputs"])
    k, Eavg = legacy["headline"]["k"], legacy["headline"]["Eavg"]
    t_min = L / U0
    values = [rtd_kernel(t, L, U0, k, Eavg) for t in (t_min, t_min * 2, t_min * 10, t_min * 100)]
    assert values == sorted(values, reverse=True)


def test_integrate_rtd_decreases_with_k_and_eavg(legacy):
    L, Uavg, U0 = _hydraulics(legacy["reactor_inputs"])
    base = integrate_rtd(L, U0, 1.0, 1.0)
    higher_k = integrate_rtd(L, U0, 2.0, 1.0)
    higher_eavg = integrate_rtd(L, U0, 1.0, 2.0)
    assert higher_k < base
    assert higher_eavg < base


def test_integrate_rtd_sweep_rows_within_expected_tolerance_of_legacy(legacy):
    """All 12 sweep rows should differ from the legacy (buggy) cached values by
    a bounded, non-negative relative amount, consistent with the documented
    off-by-one fix -- this guards against an unrelated regression slipping in."""
    L, Uavg, U0 = _hydraulics(legacy["reactor_inputs"])
    k = legacy["sweep"]["k"]
    for row in legacy["sweep"]["rows"]:
        result = integrate_rtd(L, U0, k, row["Eavg"])
        legacy_value = row["C_C0_rtd"]
        rel_diff = (legacy_value - result) / legacy_value
        assert result < legacy_value
        assert 0.0 < rel_diff < 0.15, f"Eavg={row['Eavg']}: rel_diff={rel_diff}"
