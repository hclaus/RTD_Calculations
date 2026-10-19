"""Ozone mode.

The point of most of these tests is to pin down what does NOT happen:
the RTD shape, the RTD weighting, and molecular diffusion all cancel out
of zero-order production. Those cancellations are exact, not approximate,
so they are asserted at machine precision rather than with loose
tolerances -- if any of them ever becomes approximate, something is wrong.
"""

import math

import pytest
from fastapi.testclient import TestClient

from app.core.ozone import (
    generation_rate_si,
    ozone_outlet,
    ozone_well_mixed,
    ozone_zero_order,
)
from app.core.reactor import mean_residence_time, u0, uavg
from app.core.rtd_engine import integrate_rtd_flux
from app.main import app

client = TestClient(app)

L, A, Q_CM3_MIN = 19.7, 1.1, 1200.0
RI = {"L": L, "a": A, "Q": Q_CM3_MIN, "mu": 0.00018, "rho": 0.001204, "D": 0.144}

UAVG = uavg(Q_CM3_MIN, A)
U0 = u0(Q_CM3_MIN, A)
TAU = mean_residence_time(L, UAVG)
G_SI = generation_rate_si(10.5)


def test_generation_rate_unit_conversion():
    """10.5 ppb/hr/(uW/cm^2) -> ppb/s per (mW/cm^2): x1000 for uW->mW, /3600
    for hr->s."""
    assert generation_rate_si(10.5) == pytest.approx(2.9166666666, rel=1e-9)
    assert generation_rate_si(0.0) == 0.0


def test_zero_order_sfm_integral_equals_first_moment():
    """The heart of the mode: because [O3]_batch(t) = G*Eavg*t is LINEAR in
    t, the segregated-flow integral collapses onto the RTD's first moment,
    which is tau. Integrate it numerically and confirm."""
    from scipy.integrate import quad

    eavg = 1.6

    def e_flux(t):
        return TAU**2 / (2.0 * t**3)

    sfm, _ = quad(lambda t: G_SI * eavg * t * e_flux(t), TAU / 2, math.inf)
    assert sfm == pytest.approx(G_SI * eavg * TAU, rel=1e-12)


def test_rtd_weighting_is_irrelevant_at_zero_order():
    """Flux-weighted RTD and the naive well-mixed model give bit-identical
    answers. For FIRST-order kinetics the same two differ by ~7x -- that
    contrast is the whole reason this mode exists."""
    for eavg in (0.1, 1.0, 1.6, 2.8):
        rtd = ozone_outlet(G_SI, eavg, TAU, L, U0)
        wm = ozone_well_mixed(G_SI, eavg, TAU)
        assert rtd == wm

    # the contrast case: first order, same reactor
    first_order_rtd = integrate_rtd_flux(L, U0, 1.0, 1.6)
    first_order_wm = math.exp(-1.0 * 1.6 * TAU)
    assert first_order_rtd / first_order_wm > 5.0


def test_mass_balance_identity():
    """[O3] = G*Eavg*V/Q must equal G*Eavg*tau, since tau IS V/Q. This is
    the diffusion-independence argument in one line."""
    volume = math.pi * A**2 * L
    q_cm3_s = Q_CM3_MIN / 60.0
    eavg = 1.6
    assert ozone_zero_order(G_SI, eavg, TAU) == pytest.approx(
        G_SI * eavg * volume / q_cm3_s, rel=1e-12
    )


def test_total_production_is_independent_of_flow_rate():
    """Raising Q dilutes the per-pass concentration but delivers the same
    total ozone downstream: Q * delta_C is invariant."""
    products = []
    for q_cm3_min in (600.0, 1200.0, 2400.0, 6000.0):
        tau = mean_residence_time(L, uavg(q_cm3_min, A))
        delta_c = ozone_zero_order(G_SI, 1.6, tau)
        products.append((q_cm3_min / 60.0) * delta_c)
    for p in products[1:]:
        assert p == pytest.approx(products[0], rel=1e-12)


def test_known_operating_points():
    """Values quoted in 'Ozone generation.md'."""
    assert ozone_zero_order(G_SI, 1.0, TAU) == pytest.approx(10.921, abs=0.001)
    assert ozone_zero_order(G_SI, 1.6, TAU) == pytest.approx(17.473, abs=0.001)
    assert ozone_zero_order(G_SI, 2.8, TAU) == pytest.approx(30.578, abs=0.001)


def test_output_is_linear_in_fluence_rate():
    """Zero order => exact proportionality, so doubling Eavg doubles [O3]."""
    base = ozone_zero_order(G_SI, 1.0, TAU)
    assert ozone_zero_order(G_SI, 2.0, TAU) == pytest.approx(2 * base, rel=1e-12)
    assert ozone_zero_order(G_SI, 5.0, TAU) == pytest.approx(5 * base, rel=1e-12)


# --------------------------------------------------------- photolytic loss


def test_photolytic_loss_agrees_with_the_survival_engine_at_working_kd():
    """At a realistic kd the stable expm1 form and the naive
    (G/kd)*(1 - survival) form must agree -- this pins the physics, while
    the test below pins the numerics where the naive form breaks down."""
    kd, eavg = 0.0036879, 1.6
    naive = (G_SI / kd) * (1.0 - integrate_rtd_flux(L, U0, kd, eavg))
    assert ozone_outlet(G_SI, eavg, TAU, L, U0, kd=kd) == pytest.approx(naive, rel=1e-8)


def test_photolytic_loss_reduces_output_and_vanishes_as_kd_goes_to_zero():
    """As kd -> 0 the loss term converges back onto G*Eavg*tau -- but only
    LOGARITHMICALLY, because the flux-weighted RTD's divergent second
    moment puts the leading correction in x^2*ln(x) rather than x^2. The
    tolerances below track that real behaviour; demanding linear
    convergence here would be asserting the wrong physics."""
    eavg = 1.6
    zero = ozone_zero_order(G_SI, eavg, TAU)
    assert ozone_outlet(G_SI, eavg, TAU, L, U0, kd=0.0036879) < zero

    for kd, tol in ((1e-4, 2e-3), (1e-6, 1e-5), (1e-8, 1e-7), (1e-10, 1e-8)):
        assert ozone_outlet(G_SI, eavg, TAU, L, U0, kd=kd) == pytest.approx(zero, rel=tol)

    # monotone approach: smaller kd must land strictly closer to the limit
    errors = [
        abs(ozone_outlet(G_SI, eavg, TAU, L, U0, kd=kd) / zero - 1)
        for kd in (1e-4, 1e-6, 1e-8, 1e-10)
    ]
    assert errors == sorted(errors, reverse=True)


def test_small_kd_formulation_beats_the_naive_difference():
    """Below about kd = 1e-9, (G/kd)*(1 - S) stops converging because
    (1 - S) is a catastrophic cancellation, while folding the subtraction
    inside the integral with expm1 keeps going. This is the whole reason
    ozone_outlet does its own quadrature instead of calling
    integrate_rtd_flux."""
    eavg, kd = 1.6, 1e-13
    zero = ozone_zero_order(G_SI, eavg, TAU)
    naive = (G_SI / kd) * (1.0 - integrate_rtd_flux(L, U0, kd, eavg))
    stable = ozone_outlet(G_SI, eavg, TAU, L, U0, kd=kd)
    assert abs(stable / zero - 1) < 1e-10
    assert abs(naive / zero - 1) > abs(stable / zero - 1)


def test_photolytic_loss_separates_rtd_from_well_mixed():
    """Once kd > 0 the kinetics are non-linear again, so the RTD stops
    cancelling and the two models must disagree."""
    kd, eavg = 0.0036879, 2.8
    assert ozone_outlet(G_SI, eavg, TAU, L, U0, kd=kd) != ozone_well_mixed(
        G_SI, eavg, TAU, kd=kd
    )


# --------------------------------------------------------------- endpoints


def test_headline_endpoint_known_value():
    resp = client.post("/api/ozone/headline", json={"reactor": RI, "Eavg": 1.6})
    assert resp.status_code == 200
    data = resp.json()
    assert data["ozone_ppb"] == pytest.approx(17.473, abs=0.001)
    assert data["G_si"] == pytest.approx(2.9166666666, rel=1e-9)
    assert data["Da_loss"] == 0.0
    assert data["correction_applied"] is False


def test_headline_endpoint_three_models_coincide_without_loss():
    resp = client.post("/api/ozone/headline", json={"reactor": RI, "Eavg": 2.8})
    data = resp.json()
    assert data["ozone_ppb"] == data["ozone_ppb_zero_order"] == data["ozone_ppb_well_mixed"]


def test_headline_endpoint_uses_default_generation_rate():
    """G omitted -> falls back to the measured 10.5 from config."""
    resp = client.post("/api/ozone/headline", json={"reactor": RI, "Eavg": 1.6})
    assert resp.json()["G_si"] == pytest.approx(generation_rate_si(10.5), rel=1e-12)


def test_headline_endpoint_with_photolytic_loss():
    resp = client.post(
        "/api/ozone/headline", json={"reactor": RI, "Eavg": 1.6, "kd": 0.0036879}
    )
    data = resp.json()
    assert data["ozone_ppb"] == pytest.approx(16.949, abs=0.001)
    assert data["ozone_ppb"] < data["ozone_ppb_zero_order"]
    assert data["Da_loss"] == pytest.approx(0.02209, abs=1e-5)
    # the correction is reported for transparency but must NOT be applied:
    # ozone depends on (1 - S), which amplifies table error by ~1/Da
    assert data["correction_applied"] is False
    assert data["correction_warning"] is not None


def test_headline_endpoint_rejects_bad_inputs():
    for body in (
        {"reactor": RI, "Eavg": -1.0},
        {"reactor": RI, "Eavg": 1.6, "G": -1.0},
        {"reactor": RI, "Eavg": 1.6, "kd": -0.1},
        {"reactor": {**RI, "a": -1}, "Eavg": 1.6},
    ):
        assert client.post("/api/ozone/headline", json=body).status_code == 422


def test_sweep_endpoint_returns_all_rows_and_is_linear():
    resp = client.post("/api/ozone/sweep", json={"reactor": RI})
    assert resp.status_code == 200
    rows = resp.json()["rows"]
    assert len(rows) == 12
    # exact proportionality through the origin
    slope = rows[0]["ozone_ppb"] / rows[0]["Eavg"]
    for row in rows:
        assert row["ozone_ppb"] / row["Eavg"] == pytest.approx(slope, rel=1e-12)
        assert row["ozone_ppb"] == row["ozone_ppb_well_mixed"] == row["ozone_ppb_zero_order"]


def test_sweep_endpoint_separates_models_once_kd_is_set():
    resp = client.post("/api/ozone/sweep", json={"reactor": RI, "kd": 0.0036879})
    rows = resp.json()["rows"]
    for row in rows:
        assert row["ozone_ppb"] < row["ozone_ppb_zero_order"]
        assert row["ozone_ppb"] != row["ozone_ppb_well_mixed"]


def test_diffusivity_does_not_change_the_result():
    """The headline claim, end to end: vary D over six orders of magnitude
    (bioaerosol -> gas phase) and the outlet concentration must not move."""
    values = []
    for d in (2.77e-7, 1e-3, 0.0432, 0.144, 1.0):
        resp = client.post(
            "/api/ozone/headline", json={"reactor": {**RI, "D": d}, "Eavg": 1.6}
        )
        values.append(resp.json()["ozone_ppb"])
    for v in values[1:]:
        assert v == values[0]
