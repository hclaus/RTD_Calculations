"""Solver tests for the CORRECTED build.

The solvers here work against the flux-weighted RTD and, when a theta is
supplied, the molecular-diffusion correction evaluated INSIDE the
root-find. These tests are deliberately not comparable with the legacy
area-weighted fixtures.
"""

import json
from pathlib import Path

import pytest

from app.core.diffusion import correction_factor, damkohler, theta_from
from app.core.reactor import mean_residence_time, u0, uavg
from app.core.rtd_engine import integrate_rtd_corrected, integrate_rtd_flux
from app.core.solvers import BisectionError, find_e, find_k

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture
def ref():
    return json.loads((FIXTURES / "corrected_reference_values.json").read_text())


@pytest.fixture
def geom(ref):
    ri = ref["reactor_inputs"]
    Uavg = uavg(ri["Q"], ri["a"])
    return ri["L"], Uavg, u0(ri["Q"], ri["a"]), mean_residence_time(ri["L"], Uavg), ri["a"]


# --- uncorrected (pure flux-weighted) round trips -------------------------

def test_find_k_round_trip_flux(geom):
    L, _, U0, _, _ = geom
    k_true, Eavg = 1.2, 1.6
    R = integrate_rtd_flux(L, U0, k_true, Eavg)
    assert find_k(R, L, U0, Eavg) == pytest.approx(k_true, rel=1e-6)


def test_find_e_round_trip_flux(geom):
    L, _, U0, _, _ = geom
    k, E_true = 0.87, 1.3
    R = integrate_rtd_flux(L, U0, k, E_true)
    assert find_e(R, L, U0, k) == pytest.approx(E_true, rel=1e-6)


# --- corrected round trips: the correction must be inverted consistently ---

def test_find_k_round_trip_with_correction(ref, geom):
    """If the correction were applied post-hoc instead of inside the
    root-find, these round trips would not close."""
    L, _, U0, _, _ = geom
    for case in ref["find_k_round_trip_corrected"]:
        got = find_k(case["C_C0"], L, U0, case["Eavg"], theta=case["theta"])
        assert got == pytest.approx(case["k_true"], rel=1e-6)


def test_find_e_round_trip_with_correction(ref, geom):
    L, _, U0, _, _ = geom
    for case in ref["find_e_round_trip_corrected"]:
        got = find_e(case["C_C0"], L, U0, case["k"], theta=case["theta"])
        assert got == pytest.approx(case["Eavg_true"], rel=1e-6)


def test_correction_actually_changes_the_back_solved_k(geom):
    """Guards against theta being silently ignored: solving the same
    measurement with and without the correction must differ for a
    gas-phase theta."""
    L, _, U0, tau, a = geom
    theta = theta_from(0.068, tau, a)
    R, Eavg = 0.30, 2.8
    k_plain = find_k(R, L, U0, Eavg)
    k_corr = find_k(R, L, U0, Eavg, theta=theta)
    assert k_corr != pytest.approx(k_plain, rel=1e-3)
    # ignoring diffusion over-predicts survival, so it inflates the
    # inferred rate constant
    assert k_plain > k_corr


def test_aerosol_theta_leaves_solution_unchanged(geom):
    """At aerosol diffusivity the correction is exactly 1, so passing
    theta must be a no-op."""
    L, _, U0, tau, a = geom
    theta = theta_from(2.77e-7, tau, a)
    R, Eavg = 0.05, 1.6
    assert find_k(R, L, U0, Eavg, theta=theta) == pytest.approx(
        find_k(R, L, U0, Eavg), rel=1e-9
    )


# --- bracket behaviour (unchanged semantics) ------------------------------

def test_find_k_raises_for_unreachable_target(geom):
    L, _, U0, _, _ = geom
    with pytest.raises(BisectionError):
        find_k(1.5, L, U0, 1.6)


def test_find_e_raises_for_unreachable_target(geom):
    L, _, U0, _, _ = geom
    with pytest.raises(BisectionError):
        find_e(1.5, L, U0, 0.87)


def test_find_k_custom_bracket_changes_behavior(geom):
    L, _, U0, _, _ = geom
    k_true, Eavg = 10.0, 1.6
    R = integrate_rtd_flux(L, U0, k_true, Eavg)
    assert find_k(R, L, U0, Eavg) == pytest.approx(k_true, rel=1e-6)
    with pytest.raises(BisectionError):
        find_k(R, L, U0, Eavg, k_min=0.001, k_max=5.0)
