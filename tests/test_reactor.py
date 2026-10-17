import pytest

from app.core.reactor import (
    diffusion_regime,
    flow_regime,
    mean_residence_time,
    min_residence_time,
    peclet,
    reynolds,
    u0,
    uavg,
)


@pytest.fixture
def legacy_inputs():
    return {"L": 19.7, "a": 1.1, "Q": 1200, "mu": 0.00018, "rho": 0.001204}


def test_uavg_and_u0(legacy_inputs):
    Uavg = uavg(legacy_inputs["Q"], legacy_inputs["a"])
    assert Uavg == pytest.approx(5.261320432789928, rel=1e-9)
    assert u0(legacy_inputs["Q"], legacy_inputs["a"]) == pytest.approx(2 * Uavg, rel=1e-12)


def test_reynolds(legacy_inputs):
    Uavg = uavg(legacy_inputs["Q"], legacy_inputs["a"])
    Re = reynolds(legacy_inputs["rho"], Uavg, legacy_inputs["a"], legacy_inputs["mu"])
    assert Re == pytest.approx(77.42325312429978, rel=1e-9)


def test_mean_and_min_residence_time(legacy_inputs):
    Uavg = uavg(legacy_inputs["Q"], legacy_inputs["a"])
    U0 = u0(legacy_inputs["Q"], legacy_inputs["a"])
    assert mean_residence_time(legacy_inputs["L"], Uavg) == pytest.approx(3.7443072041809953, rel=1e-9)
    assert min_residence_time(legacy_inputs["L"], U0) == pytest.approx(1.8721536020904976, rel=1e-9)


@pytest.mark.parametrize(
    "Re, expected",
    [
        (77.4, "Laminar"),
        (2299.9, "Laminar"),
        (2300.0, "Transitional (uncertain)"),
        (3000.0, "Transitional (uncertain)"),
        (4000.0, "Transitional (uncertain)"),
        (4000.1, "Turbulent"),
        (10000.0, "Turbulent"),
    ],
)
def test_flow_regime_boundaries(Re, expected):
    assert flow_regime(Re) == expected


@pytest.mark.parametrize(
    "Pe, expected_substring",
    [
        (1e5, "Advection-dominated"),
        (1e5 - 1, "Diffusion may be significant"),
        (3.74e8, "Advection-dominated"),
        (10.0, "Diffusion may be significant"),
    ],
)
def test_diffusion_regime_boundaries(Pe, expected_substring):
    assert diffusion_regime(Pe).startswith(expected_substring)


def test_peclet_matches_pptx_worked_example(legacy_inputs):
    # Theoretical Foundation.pptx slide 9 worked example: U=5.26 cm/s, L=19.7 cm,
    # D<=2.77e-7 cm^2/s -> Pe = 3.74e8
    Uavg = uavg(legacy_inputs["Q"], legacy_inputs["a"])
    Pe = peclet(Uavg, legacy_inputs["L"], 2.77e-7)
    assert Pe == pytest.approx(3.74e8, rel=1e-2)
