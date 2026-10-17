import pytest
from fastapi.testclient import TestClient

from app.core.chemical import (
    PLANCK_CONSTANT_J_S,
    SPEED_OF_LIGHT_M_S,
    diffusivity_in_air,
    kp,
    photon_energy_j,
    photon_flux_constant,
)
from app.main import app

client = TestClient(app)


def test_photon_energy_matches_pptx_worked_example():
    # Theoretical Foundation.pptx slide 15: Up(222nm) = 8.948e-19 J/photon.
    up = photon_energy_j(PLANCK_CONSTANT_J_S, SPEED_OF_LIGHT_M_S, 222.0)
    assert up == pytest.approx(8.948e-19, rel=1e-3)


def test_kp_formula_composition():
    up = photon_energy_j(PLANCK_CONSTANT_J_S, SPEED_OF_LIGHT_M_S, 222.0)
    j = photon_flux_constant(up)
    value = kp(sigma=1.32e-17, phi=1.0, J=j)
    # sigma * phi * (1/(1000*up)) -- sanity-check order of magnitude and sign.
    assert value > 0
    assert value == pytest.approx(1.32e-17 * 1.0 * j, rel=1e-12)


def test_compute_kp_endpoint():
    resp = client.post("/api/chemical/kp", json={"sigma": 1.32e-17, "phi": 1.0, "wavelength_nm": 222.0})
    assert resp.status_code == 200
    data = resp.json()
    assert data["photon_energy_j"] == pytest.approx(8.948e-19, rel=1e-3)
    assert data["kp"] == pytest.approx(0.014751760253892977, rel=1e-6)


def test_compute_kp_rejects_non_positive_wavelength():
    resp = client.post("/api/chemical/kp", json={"sigma": 1.32e-17, "phi": 1.0, "wavelength_nm": -5})
    assert resp.status_code == 422


def test_kp_reusable_in_existing_headline_endpoint():
    """kp should plug directly into the generic /headline endpoint as 'k',
    proving the RTD engine is genuinely domain-agnostic -- no new engine
    code was needed for chemical mode."""
    kp_resp = client.post("/api/chemical/kp", json={"sigma": 1.32e-17, "phi": 1.0, "wavelength_nm": 222.0})
    kp_value = kp_resp.json()["kp"]

    ri = {"L": 19.7, "a": 1.1, "Q": 1200, "mu": 0.00018, "rho": 0.001204, "D": 2.77e-7}
    resp = client.post("/api/microbiology/headline", json={"reactor": ri, "k": kp_value, "Eavg": 2.8})
    assert resp.status_code == 200
    assert 0 < resp.json()["C_C0_reactor"] < 1


def test_back_solve_eavg_returns_eavg_for_each_value():
    kp_resp = client.post("/api/chemical/kp", json={"sigma": 1.32e-17, "phi": 1.0, "wavelength_nm": 222.0})
    kp_value = kp_resp.json()["kp"]

    ri = {"L": 19.7, "a": 1.1, "Q": 1200, "mu": 0.00018, "rho": 0.001204, "D": 2.77e-7}
    body = {"reactor": ri, "k": kp_value, "c_c0_values": [0.95, 0.90, 0.85, 0.80]}
    resp = client.post("/api/chemical/back-solve-eavg", json=body)
    assert resp.status_code == 200
    data = resp.json()
    assert len(data["rows"]) == 4
    for row in data["rows"]:
        assert row["Eavg"] is not None
        assert row["error"] is None
        assert row["Eavg"] > 0


def test_back_solve_eavg_default_values():
    ri = {"L": 19.7, "a": 1.1, "Q": 1200, "mu": 0.00018, "rho": 0.001204, "D": 2.77e-7}
    resp = client.post("/api/chemical/back-solve-eavg", json={"reactor": ri, "k": 0.0147})
    assert resp.status_code == 200
    data = resp.json()
    assert [row["C_C0"] for row in data["rows"]] == [0.95, 0.90, 0.85, 0.80]


def test_back_solve_eavg_one_unreachable_value_does_not_fail_whole_batch():
    ri = {"L": 19.7, "a": 1.1, "Q": 1200, "mu": 0.00018, "rho": 0.001204, "D": 2.77e-7}
    body = {"reactor": ri, "k": 0.0147, "c_c0_values": [0.95, 0.01, 0.90]}  # 0.01 unreachable at this kp
    resp = client.post("/api/chemical/back-solve-eavg", json=body)
    assert resp.status_code == 200
    data = resp.json()
    assert len(data["rows"]) == 3
    assert data["rows"][0]["Eavg"] is not None
    assert data["rows"][1]["Eavg"] is None
    assert data["rows"][1]["error"] is not None
    assert data["rows"][2]["Eavg"] is not None


def test_back_solve_eavg_round_trips_with_find_e_directly():
    """Cross-check the API against calling find_e directly, to confirm the
    route isn't silently using the wrong solver."""
    from app.core.reactor import u0, uavg
    from app.core.rtd_engine import integrate_rtd
    from app.core.solvers import find_e

    L, a, Q = 19.7, 1.1, 1200
    Uavg, U0 = uavg(Q, a), u0(Q, a)
    kp_value = 0.0147
    true_eavg = 1.5
    c_c0 = integrate_rtd(L, U0, kp_value, true_eavg)
    expected_eavg = find_e(c_c0, L, U0, kp_value)

    ri = {"L": L, "a": a, "Q": Q, "mu": 0.00018, "rho": 0.001204, "D": 2.77e-7}
    resp = client.post(
        "/api/chemical/back-solve-eavg",
        json={"reactor": ri, "k": kp_value, "c_c0_values": [c_c0]},
    )
    assert resp.json()["rows"][0]["Eavg"] == pytest.approx(expected_eavg, rel=1e-9)


def test_diffusivity_in_air_matches_reference_compounds_within_tolerance():
    """The power-law fit should reproduce the slide 10 reference compounds
    (that it was fit to) within a few percent -- this is the whole point of
    using a fit instead of an unfounded generic formula."""
    reference_pairs = [
        (106.12, 0.073),  # Benzaldehyde
        (120.19, 0.070),  # 1,2,4-Trimethylbenzene
        (104.15, 0.071),  # Styrene
        (319.62, 0.068),  # Carbon Tetrabromide (CBr4)
    ]
    for mw, actual in reference_pairs:
        predicted = diffusivity_in_air(mw)
        assert predicted == pytest.approx(actual, rel=0.03)


def test_diffusivity_in_air_decreases_with_molecular_weight():
    assert diffusivity_in_air(50.0) > diffusivity_in_air(500.0)


def test_compute_diffusivity_endpoint():
    resp = client.post("/api/chemical/diffusivity", json={"molecular_weight": 319.62})
    assert resp.status_code == 200
    assert resp.json()["diffusivity_cm2_s"] == pytest.approx(diffusivity_in_air(319.62), rel=1e-9)


def test_compute_diffusivity_endpoint_default_value():
    """No molecular_weight supplied -> falls back to the CBr4 default."""
    resp = client.post("/api/chemical/diffusivity", json={})
    assert resp.status_code == 200
    assert resp.json()["diffusivity_cm2_s"] == pytest.approx(diffusivity_in_air(319.62), rel=1e-6)


def test_compute_diffusivity_rejects_non_positive_molecular_weight():
    resp = client.post("/api/chemical/diffusivity", json={"molecular_weight": -10})
    assert resp.status_code == 422


def test_diffusivity_feeds_peclet_number_via_existing_headline_endpoint():
    """The computed diffusivity should plug straight into the existing
    shared 'D' reactor input and produce a sensible Peclet number -- no new
    Pe computation path was needed, confirming reuse of reactor.py."""
    diff_resp = client.post("/api/chemical/diffusivity", json={"molecular_weight": 319.62})
    d_value = diff_resp.json()["diffusivity_cm2_s"]

    ri = {"L": 19.7, "a": 1.1, "Q": 1200, "mu": 0.00018, "rho": 0.001204, "D": d_value}
    resp = client.post("/api/microbiology/headline", json={"reactor": ri, "k": 0.0147, "Eavg": 2.8})
    assert resp.status_code == 200
    data = resp.json()
    # Gas-phase D (~0.07 cm^2/s) is ~6 orders of magnitude larger than the
    # bioaerosol default (2.77e-7), so Pe should be correspondingly smaller.
    assert data["derived"]["Pe"] < 1e4
