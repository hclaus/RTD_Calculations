"""Shared reactor hydraulics for a laminar-flow cylindrical UV reactor.

Formulas follow Taylor (1953) laminar (parabolic) pipe flow, as presented in
"UV-C Dose-Response Behavior of Airborne (Aerosolized) Challenge Agents I:
Theoretical Foundation" (Blatchley et al.), slides 9/21/23. Reused as-is by
both the microbiological (Phase 1) and chemical/actinometry (Phase 2) engines.
"""

import math

# Standard pipe-flow Reynolds-number thresholds (textbook values; exact
# transition point is empirical/approximate, not a sharp boundary).
_RE_LAMINAR_MAX = 2300.0
_RE_TURBULENT_MIN = 4000.0

# Peclet-number threshold above which advection dominates diffusion, per
# Bowers and Miller (2025) as cited on slide 9 of Theoretical Foundation.pptx.
_PE_ADVECTION_THRESHOLD = 1e5


def uavg(Q: float, a: float) -> float:
    """Mean (cross-section-averaged) velocity, cm/s. Q in cm^3/min, a in cm."""
    return Q / 60.0 / (math.pi * a**2)


def u0(Q: float, a: float) -> float:
    """Centerline (maximum) velocity for laminar/parabolic flow, cm/s."""
    return 2.0 * uavg(Q, a)


def reynolds(rho: float, Uavg: float, a: float, mu: float) -> float:
    """Reynolds number using the tube diameter (2a) as the length scale."""
    return rho * Uavg * 2.0 * a / mu


def peclet(Uavg: float, L: float, D: float) -> float:
    """Peclet number Pe = U*L/D, comparing advective to diffusive transport."""
    return Uavg * L / D


def mean_residence_time(L: float, Uavg: float) -> float:
    """Mean residence time q = L/Uavg, s."""
    return L / Uavg


def min_residence_time(L: float, U0: float) -> float:
    """Minimum (centerline) residence time L/U0, s -- the lower integration bound."""
    return L / U0


def flow_regime(Re: float) -> str:
    """Classify flow regime from the Reynolds number using standard pipe-flow thresholds."""
    if Re < _RE_LAMINAR_MAX:
        return "Laminar"
    if Re <= _RE_TURBULENT_MIN:
        return "Transitional (uncertain)"
    return "Turbulent"


def diffusion_regime(Pe: float, threshold: float = _PE_ADVECTION_THRESHOLD) -> str:
    """Classify whether advection dominates diffusion from the Peclet number."""
    if Pe >= threshold:
        return "Advection-dominated (diffusion negligible) — RTD model assumptions hold"
    return "Diffusion may be significant — zero-diffusion RTD assumption may not hold"
