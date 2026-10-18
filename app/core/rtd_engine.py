"""Core RTD-weighted microbial survival model.

Physical model (confirmed against Theoretical Foundation.pptx):
  - Taylor (1953) closed-form laminar-flow residence time distribution
    density function: E(t) = L/(U0*t^2) for t > L/U0 (slide 23).
  - First-order (Chick-Watson) batch inactivation kinetics:
    (C/C0)_batch(t) = exp(-k*Eavg*t).
  - Segregated Flow Model integral (slide 13):
    (C/C0)_reactor = integral_0^inf (C/C0)_batch(t) * E(t) dt
    which, since E(t)=0 for t <= L/U0, reduces to integrating over
    [L/U0, infinity). The original VBA (Module1.bas) truncated this at
    100*L/U0 as a practical "infinity"; we keep that same truncation.

Numerics: the original VBA used a fixed 10,000-point trapezoidal rule with
a genuine off-by-one (its loop re-initializes `t` to `tMin` on the first
interior iteration, so it double-weights the t=tMin endpoint and silently
drops the t=tMax-h point). That bug is NOT reproduced here -- per explicit
product decision, this port uses scipy.integrate.quad (adaptive quadrature)
instead, which is both more accurate and faster. Expect small (~4-10% in
tested cases) differences from the legacy spreadsheet's cached values as a
direct, documented consequence -- see tests/fixtures/spreadsheet_known_values.json.
"""

import math

from scipy.integrate import quad


def rtd_kernel(t: float, L: float, U0: float, k: float, Eavg: float) -> float:
    """RTD(t) = exp(-k*Eavg*t) * L/U0 / t**2, for t > 0."""
    return math.exp(-k * Eavg * t) * L / U0 / t**2


def integrate_rtd(L: float, U0: float, k: float, Eavg: float) -> float:
    """Predicted overall survival ratio N/N0 (= C/C0), integrating rtd_kernel
    over [L/U0, 100*L/U0] via adaptive quadrature."""
    t_min = L / U0
    t_max = 100.0 * t_min
    value, _abserr = quad(rtd_kernel, t_min, t_max, args=(L, U0, k, Eavg))
    return value


# ---------------------------------------------------------------------------
# Flux-weighted ("mixed-cup") formulation -- see "RTD calculation review.md".
#
# The functions above use the AREA-weighted laminar RTD, E(t) = L/(U0 t^2),
# as given in Theoretical Foundation.pptx (slides 22-23) and the legacy
# spreadsheet. That distribution has a DIVERGENT mean residence time, which
# cannot describe a reactor of finite volume at finite flow.
#
# The functions below implement the flux-weighted (mixed-cup) laminar RTD,
#     E(t) = tau^2 / (2 t^3),   t >= tau/2,   tau = L/Uavg = 2L/U0
# the standard form in the chemical reaction engineering literature
# (Levenspiel; Fogler). Its mean residence time equals tau = V/Q exactly.
# This is the quantity a bioassay measures, since collected effluent is
# weighted by the flow each streamline delivers, not by cross-sectional area.
#
# These are ADDITIONS: nothing above is modified, so existing app behaviour
# is unchanged pending expert confirmation of the finding.
# ---------------------------------------------------------------------------


def rtd_density_area(t: float, L: float, U0: float) -> float:
    """Area-weighted laminar RTD density (as currently used), 1/s."""
    return L / U0 / t**2


def rtd_density_flux(t: float, L: float, U0: float) -> float:
    """Flux-weighted (mixed-cup) laminar RTD density, 1/s."""
    tau = 2.0 * L / U0
    return tau**2 / (2.0 * t**3)


def rtd_kernel_flux(t: float, L: float, U0: float, k: float, Eavg: float) -> float:
    """First-order survival weighted by the flux-weighted RTD."""
    return math.exp(-k * Eavg * t) * rtd_density_flux(t, L, U0)


def integrate_rtd_flux(L: float, U0: float, k: float, Eavg: float) -> float:
    """Predicted survival ratio C/C0 using the flux-weighted laminar RTD.

    Integrated to infinity rather than a 100*L/U0 cutoff: the t^-3 tail
    converges cleanly, so no truncation is needed (unlike the area-weighted
    kernel, whose heavier t^-2 tail makes truncation unavoidable)."""
    t_min = L / U0
    value, _abserr = quad(rtd_kernel_flux, t_min, math.inf, args=(L, U0, k, Eavg))
    return value


def survival_well_mixed(L: float, Uavg: float, k: float, Eavg: float) -> float:
    """Naive 'well-mixed' survival: every element receives the same average
    dose Davg = Eavg * q, with q = L/Uavg. Ignores the residence time
    distribution entirely -- the comparison baseline already shown in the
    app's sweep table."""
    q = L / Uavg
    return math.exp(-k * Eavg * q)


def integrate_rtd_corrected(
    L: float, U0: float, k: float, Eavg: float, theta: float | None = None
) -> float:
    """Survival ratio C/C0 using the flux-weighted laminar RTD, optionally
    multiplied by the tabulated molecular-diffusion correction.

    theta = D*tau/a^2. Pass None (or 0) to skip the correction entirely and
    get the pure zero-diffusion flux-weighted result -- useful for
    comparison and for the aerosol case where the correction is 1 anyway.

    Note tau is derived here as 2*L/U0 (= L/Uavg), so no separate flow
    argument is needed.
    """
    base = integrate_rtd_flux(L, U0, k, Eavg)
    if not theta:
        return base

    from app.core.diffusion import correction_factor  # local: avoids a cycle

    tau = 2.0 * L / U0
    return base * correction_factor(theta, k * Eavg * tau)


def mean_residence_time_from_rtd(
    L: float, U0: float, weighting: str = "flux", t_max_factor: float | None = None
) -> float:
    """Mean residence time obtained BY INTEGRATING the RTD itself,
    <t> = integral t*E(t) dt -- as opposed to the algebraic V/Q used in
    reactor.mean_residence_time(). Comparing the two is the diagnostic that
    distinguishes a physically consistent RTD from an inconsistent one.

    weighting:     "flux" or "area".
    t_max_factor:  upper limit as a multiple of t_min. None means infinity.
                   The area-weighted first moment DIVERGES, so a finite
                   factor must be supplied for it (use 100.0 to reproduce
                   what the app's own truncation implies).
    """
    t_min = L / U0
    if weighting == "flux":
        density = rtd_density_flux
    elif weighting == "area":
        density = rtd_density_area
    else:
        raise ValueError("weighting must be 'flux' or 'area'")

    upper = math.inf if t_max_factor is None else t_max_factor * t_min
    if weighting == "area" and upper == math.inf:
        raise ValueError(
            "The area-weighted RTD has a divergent first moment; supply "
            "t_max_factor to evaluate it over a finite window."
        )
    value, _abserr = quad(lambda t: t * density(t, L, U0), t_min, upper, limit=400)
    return value
