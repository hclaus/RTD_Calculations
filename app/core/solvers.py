"""Inverse (back-solving) routines: given a measured survival ratio, find
the rate constant (find_k) or average fluence rate (find_e) consistent
with it under the RTD model.

This is the CORRECTED build: it uses the flux-weighted laminar RTD (see
"RTD calculation review.md") and optionally applies the tabulated
molecular-diffusion correction (see "Diffusion correction.md").

CRITICAL -- why the correction lives inside the residual
--------------------------------------------------------
The correction factor F is a function of (theta, Da) with

    Da = k * Ebar * tau

so it depends on whichever quantity is being solved for. Applying F once
as a post-hoc multiplier on an already-solved answer would be wrong: the
solver would converge against the uncorrected model and then be nudged,
rather than converging against the corrected model. Here F is evaluated
inside the residual at every iterate, so brentq converges to the root of
the *corrected* equation.

Two properties keep this safe:
  * theta is CONSTANT during the solve (it depends only on D, flow and
    geometry, never on k or Ebar), so only the Da axis of the table is
    traversed and the interpolant stays smooth.
  * Monotonicity is preserved -- raising k lowers both the baseline
    survival and F -- so the residual remains monotone and the existing
    bracket semantics are unchanged.
"""

from scipy.optimize import brentq

from app.core.rtd_engine import integrate_rtd_corrected


class BisectionError(RuntimeError):
    """Raised when the measured target isn't bracketed by the current solver
    bracket (i.e. f(a) and f(b) don't have opposite signs) -- the bracket
    needs widening in the advanced-settings panel."""


def find_k(
    R: float,
    L: float,
    U0: float,
    Eavg: float,
    k_min: float = 0.001,
    k_max: float = 50.0,
    theta: float | None = None,
) -> float:
    """Back-solve k so that the corrected model reproduces the measured R.

    theta: pass the reactor's D*tau/a^2 to include the diffusion
    correction, or None to solve against the pure flux-weighted RTD.
    """
    def f(k: float) -> float:
        return R - integrate_rtd_corrected(L, U0, k, Eavg, theta)

    try:
        return brentq(f, k_min, k_max, xtol=1e-10, maxiter=100)
    except ValueError as exc:
        raise BisectionError(
            f"No sign change for k in [{k_min}, {k_max}]: measured R={R} is not "
            "reachable with this bracket. Widen the k range in advanced settings."
        ) from exc


def find_e(
    R: float,
    L: float,
    U0: float,
    k: float,
    eavg_min: float = 0.0001,
    eavg_max: float = 3.0,
    theta: float | None = None,
) -> float:
    """Back-solve Eavg so that the corrected model reproduces the measured R.

    theta: pass the reactor's D*tau/a^2 to include the diffusion
    correction, or None to solve against the pure flux-weighted RTD.
    """
    def f(Eavg: float) -> float:
        return R - integrate_rtd_corrected(L, U0, k, Eavg, theta)

    try:
        return brentq(f, eavg_min, eavg_max, xtol=1e-10, maxiter=100)
    except ValueError as exc:
        raise BisectionError(
            f"No sign change for Eavg in [{eavg_min}, {eavg_max}]: measured R={R} is "
            "not reachable with this bracket. Widen the Eavg range in advanced settings."
        ) from exc
