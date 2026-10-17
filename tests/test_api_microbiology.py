import json
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from app.main import app

FIXTURES = Path(__file__).parent / "fixtures"
client = TestClient(app)


@pytest.fixture
def scipy_ref():
    return json.loads((FIXTURES / "scipy_reference_values.json").read_text())


def test_headline_matches_scipy_reference(scipy_ref):
    ri = scipy_ref["reactor_inputs"]
    body = {
        "reactor": ri,
        "k": scipy_ref["headline"]["k"],
        "Eavg": scipy_ref["headline"]["Eavg"],
        "measured_C_C0": scipy_ref["headline"]["measured_C_C0"],
    }
    resp = client.post("/api/microbiology/headline", json=body)
    assert resp.status_code == 200
    data = resp.json()

    assert data["C_C0_reactor"] == pytest.approx(scipy_ref["headline"]["C_C0_reactor"], rel=1e-6)
    assert data["estimated_k"] == pytest.approx(scipy_ref["headline"]["estimated_k"], rel=1e-6)
    assert data["residual"] == pytest.approx(0.0, abs=1e-9)

    derived = data["derived"]
    expected = scipy_ref["derived"]
    for key in ("Uavg", "U0", "Re", "Pe", "q", "L_over_U0"):
        assert derived[key] == pytest.approx(expected[key], rel=1e-6)
    assert derived["flow_regime"] == expected["flow_regime"]
    assert derived["diffusion_regime"] == expected["diffusion_regime"]


def test_headline_without_measured_c_c0_skips_solve(scipy_ref):
    ri = scipy_ref["reactor_inputs"]
    body = {"reactor": ri, "k": scipy_ref["headline"]["k"], "Eavg": scipy_ref["headline"]["Eavg"]}
    resp = client.post("/api/microbiology/headline", json=body)
    assert resp.status_code == 200
    data = resp.json()
    assert data["estimated_k"] is None
    assert data["residual"] is None


def test_headline_rejects_negative_length(scipy_ref):
    ri = dict(scipy_ref["reactor_inputs"])
    ri["L"] = -5.0
    body = {"reactor": ri, "k": 0.87, "Eavg": 1.6}
    resp = client.post("/api/microbiology/headline", json=body)
    assert resp.status_code == 422


def test_headline_rejects_invalid_solver_bracket(scipy_ref):
    ri = scipy_ref["reactor_inputs"]
    body = {
        "reactor": ri,
        "k": 0.87,
        "Eavg": 1.6,
        "solver_settings": {"k_min": 10.0, "k_max": 5.0, "eavg_min": 0.0001, "eavg_max": 3.0},
    }
    resp = client.post("/api/microbiology/headline", json=body)
    assert resp.status_code == 422


def test_headline_returns_422_on_unreachable_bisection_target(scipy_ref):
    ri = scipy_ref["reactor_inputs"]
    body = {"reactor": ri, "k": 0.87, "Eavg": 1.6, "measured_C_C0": 0.9999999}
    resp = client.post("/api/microbiology/headline", json=body)
    assert resp.status_code == 422
    assert "bracket" in resp.json()["detail"].lower() or "sign" in resp.json()["detail"].lower()


def test_sweep_returns_all_requested_rows_not_truncated(scipy_ref):
    ri = scipy_ref["reactor_inputs"]
    eavg_values = [row["Eavg"] for row in scipy_ref["sweep"]["rows"]]
    body = {"reactor": ri, "k": scipy_ref["sweep"]["k"], "eavg_values": eavg_values}
    resp = client.post("/api/microbiology/sweep", json=body)
    assert resp.status_code == 200
    data = resp.json()
    assert len(data["rows"]) == 12

    for row, expected in zip(data["rows"], scipy_ref["sweep"]["rows"]):
        assert row["C_C0_rtd"] == pytest.approx(expected["C_C0_rtd"], rel=1e-6)
        assert row["Davg"] == pytest.approx(expected["Davg"], rel=1e-6)


def test_back_solve_k_returns_k_for_each_pair(scipy_ref):
    ri = scipy_ref["reactor_inputs"]
    Eavg = scipy_ref["headline"]["Eavg"]
    k_true = scipy_ref["headline"]["k"]

    # Build a measured C/C0 value from a known-good forward prediction so we
    # can assert the back-solved k round-trips to the true value.
    sweep_body = {"reactor": ri, "k": k_true, "eavg_values": [Eavg]}
    sweep_resp = client.post("/api/microbiology/sweep", json=sweep_body)
    c_c0_known = sweep_resp.json()["rows"][0]["C_C0_rtd"]

    body = {
        "reactor": ri,
        "c_c0_values": [c_c0_known, 0.5, 0.05],
        "eavg_values": [Eavg, 1.2, 0.8],  # each C/C0 measured at its own Eavg
    }
    resp = client.post("/api/microbiology/back-solve-k", json=body)
    assert resp.status_code == 200
    data = resp.json()
    assert len(data["rows"]) == 3

    first = data["rows"][0]
    assert first["C_C0"] == pytest.approx(c_c0_known, rel=1e-9)
    assert first["Eavg"] == pytest.approx(Eavg, rel=1e-9)
    assert first["k"] == pytest.approx(k_true, rel=1e-5)
    assert first["error"] is None

    for row, expected_eavg in zip(data["rows"], [Eavg, 1.2, 0.8]):
        assert row["Eavg"] == pytest.approx(expected_eavg, rel=1e-9)
        assert row["k"] is not None
        assert row["error"] is None


def test_back_solve_k_default_values(scipy_ref):
    ri = scipy_ref["reactor_inputs"]
    body = {"reactor": ri}
    resp = client.post("/api/microbiology/back-solve-k", json=body)
    assert resp.status_code == 200
    data = resp.json()
    assert [row["C_C0"] for row in data["rows"]] == [0.5, 0.10, 0.05, 0.01, 0.005]
    assert [row["Eavg"] for row in data["rows"]] == [1.6, 1.6, 1.6, 1.6, 1.6]


def test_back_solve_k_one_unreachable_value_does_not_fail_whole_batch(scipy_ref):
    ri = scipy_ref["reactor_inputs"]
    body = {
        "reactor": ri,
        "c_c0_values": [0.5, 1.5, 0.05],  # 1.5 is physically unreachable (>1)
        "eavg_values": [1.6, 1.6, 1.6],
    }
    resp = client.post("/api/microbiology/back-solve-k", json=body)
    assert resp.status_code == 200
    data = resp.json()
    assert len(data["rows"]) == 3
    assert data["rows"][0]["k"] is not None
    assert data["rows"][1]["k"] is None
    assert data["rows"][1]["error"] is not None
    assert data["rows"][2]["k"] is not None


def test_back_solve_k_rejects_empty_values_list(scipy_ref):
    ri = scipy_ref["reactor_inputs"]
    body = {"reactor": ri, "c_c0_values": [], "eavg_values": []}
    resp = client.post("/api/microbiology/back-solve-k", json=body)
    assert resp.status_code == 422


def test_back_solve_k_rejects_mismatched_list_lengths(scipy_ref):
    ri = scipy_ref["reactor_inputs"]
    body = {"reactor": ri, "c_c0_values": [0.5, 0.1, 0.05], "eavg_values": [1.6, 1.2]}
    resp = client.post("/api/microbiology/back-solve-k", json=body)
    assert resp.status_code == 422
    assert "same length" in resp.json()["detail"][0]["msg"].lower()


def test_get_variables():
    resp = client.get("/api/variables")
    assert resp.status_code == 200
    data = resp.json()
    ids = {v["id"] for v in data}
    assert {"L", "a", "Q", "mu", "rho", "D", "k", "Eavg", "measured_C_C0"} <= ids


def test_get_solver_settings_defaults():
    resp = client.get("/api/solver-settings/defaults")
    assert resp.status_code == 200
    data = resp.json()
    assert data == {"k_min": 0.001, "k_max": 50.0, "eavg_min": 0.0001, "eavg_max": 3.0}
