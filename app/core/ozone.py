"""Ozone generation at 222 nm.

Unlike inactivation, which is FIRST order in the thing being destroyed,
ozone production is ZERO order: atmospheric O2 is in vast excess and a few
tens of ppb of product never depletes it, so

    d[O3]/dt = G * Eavg        ->   [O3]_batch(t) = G * Eavg * t

The batch solution is LINEAR in t. Feeding that through the segregated
flow model therefore picks out only the FIRST MOMENT of the RTD:

    [O3]_out = int G*Eavg*t * E_rtd(t) dt = G*Eavg * <t> = G*Eavg*tau

The RTD shape cancels completely, and so does diffusion -- the latter
follows from mass balance alone (integrate the advection-diffusion
equation over the cross-section: with no wall flux the diffusive term
vanishes and [O3]_out = G*Eavg*V/Q = G*Eavg*tau for any D).

That does NOT mean the ozone field is uniform. Without diffusion the
radial profile is c(r) = G*Eavg*L/u(r), which diverges at the wall -- the
near-wall fluid is genuinely loaded with ozone. The mixed-cup average
survives this because the slow fluid carries more ozone (~1/u) but
delivers less of itself (~u), and the product is constant. See
'Ozone generation.md' for the derivation, simulations and figures.

Back-solving needs no solver here: at kd == 0 the forward model is linear,
so Eavg = [O3] / (G*tau) directly -- unlike the microbiological and
chemical modes, which bracket-and-bisect. Hence no back-solve route.

Two refinements restore the RTD's relevance, both optional:
  * photolytic loss (kd > 0) -- implemented here; off by default because
    O3 photolysis at 222 nm is largely a null cycle and the measured G is
    almost certainly already a net rate.
  * wall destruction -- not implemented; would need a measured deposition
    velocity for the actual tube material.
"""

import math

from scipy.integrate import quad

from app.core.rtd_engine import rtd_density_flux

# 1 hr / 3600 s, and 1000 uW/cm^2 per mW/cm^2.
_PPB_HR_UW_TO_PPB_S_MW = 1000.0 / 3600.0


def generation_rate_si(G_ppb_per_hr_per_uW: float) -> float:
    """Convert the measured generation rate to working units.

    Input:  ppb/hr per (uW/cm^2)   -- how the measurement is reported
    Output: ppb/s  per (mW/cm^2)   -- matches Eavg in mW/cm^2 and tau in s

    The project's measured value, 10.5 ppb/hr/(uW/cm^2), becomes
    2.9167 ppb/s per (mW/cm^2).
    """
    return G_ppb_per_hr_per_uW * _PPB_HR_UW_TO_PPB_S_MW


def ozone_zero_order(G_si: float, Eavg: float, tau: float) -> float:
    """Outlet [O3] in ppb with no loss term: exactly G*Eavg*tau.

    Independent of the RTD weighting and of diffusion. This is the
    headline prediction.
    """
    return G_si * Eavg * tau


def ozone_outlet(
    G_si: float,
    Eavg: float,
    tau: float,
    L: float,
    U0: float,
    kd: float = 0.0,
) -> float:
    """Outlet [O3] in ppb.

    kd == 0 (default): the zero-order result G*Eavg*tau.

    kd  > 0: first-order photolytic loss makes the batch solution
    non-linear again, so the RTD genuinely re-enters:

        [O3]_batch(t) = (G/kd) * (1 - exp(-kd*Eavg*t))
        [O3]_out      = (G/kd) * (1 - int exp(-kd*Eavg*t) E_rtd(t) dt)

    NUMERICS: the bracketed form above is evaluated by folding the
    subtraction INSIDE the integral,

        [O3]_out = (G/kd) * int -expm1(-kd*Eavg*t) * E_flux(t) dt

    rather than calling the survival engine and subtracting from 1. The
    two are algebraically identical and agree for every kd anyone would
    actually enter; they diverge only in the deep limit, where S -> 1 and
    (1 - S) becomes a catastrophic cancellation. Measured against the
    analytic kd -> 0 limit, relative error is:

        kd      naive (1-S)     this form
        1e-8      9.6e-8          9.6e-8
        1e-10     2.5e-7          9.6e-10
        1e-12     1.6e-5          9.6e-12
        1e-14     1.1e-3          9.6e-14

    so the naive form stops converging below about kd = 1e-9 while this
    one keeps going. Cheap insurance, and expm1 also keeps the integrand
    accurate for small arguments.

    CONVERGENCE TO THE ZERO-ORDER LIMIT IS LOGARITHMIC, not linear: the
    flux-weighted RTD has a divergent second moment, so expanding
    1 - 2*E3(x) gives a leading correction in x^2*ln(x) rather than x^2.
    In practice a loss constant as small as kd = 1e-4 still shifts the
    answer by 0.135%, far more than Da/2 would suggest. Physically, the
    long tail means the slowest fluid has had the most time both to
    accumulate ozone AND to lose it, so even a tiny sink bites.

    DELIBERATELY UNCORRECTED FOR DIFFUSION, unlike the microbiological and
    chemical modes. The tabulated factor F multiplies the SURVIVAL S, but
    ozone depends on (1 - S), and at the Da this term operates at
    (kd*Eavg*tau ~ 0.02) that inverts into a ~45x amplification:

        d[O3]/[O3] = -S*dF/(1 - S*F) ~ -45 * dF

    The table is good to ~0.2%, which would become ~9% here. Worse, at
    Da -> 0 the correction must be exactly 1 on physical grounds (no
    conversion, nothing to correct) while the interpolated surface returns
    0.9996 -- so the amplified quantity is pure interpolation residue, not
    signal. The flux-weighted RTD alone is the defensible choice.
    """
    if kd <= 0:
        return ozone_zero_order(G_si, Eavg, tau)

    def integrand(t):
        return -math.expm1(-kd * Eavg * t) * rtd_density_flux(t, L, U0)

    converted, _abserr = quad(integrand, L / U0, math.inf)
    return (G_si / kd) * converted


def ozone_well_mixed(G_si: float, Eavg: float, tau: float, kd: float = 0.0) -> float:
    """Naive well-mixed (single mean dose) baseline, for comparison.

    At kd == 0 this is IDENTICAL to ozone_outlet -- not approximately, but
    algebraically. That identity is the central result and is asserted in
    tests/test_ozone.py. It only diverges once kd > 0 makes the kinetics
    non-linear.
    """
    if kd <= 0:
        return ozone_zero_order(G_si, Eavg, tau)
    return (G_si / kd) * -math.expm1(-kd * Eavg * tau)
