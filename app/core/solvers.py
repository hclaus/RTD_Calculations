"""Inverse (back-solving) routines: given a measured bioassay survival ratio,
find the rate constant (find_k) or average fluence rate (find_e) consistent
with it under the RTD model.

Ports Module1.bas's Find_k/Find_E, which used a hand-rolled bisection search
over hardcoded brackets. Per explicit product decision, scipy.optimize.brentq
is used instead (same bracket semantics, more robust/efficient than manual
bisection) -- see rtd_engine.py for the matching decision on Integrate_RTD.

Note: the VBA Find_k took an unused 4th argument (a vestigial, dead
parameter -- it was never read in the function body) and Find_E likewise
took an unused Eavg seed. Those dead parameters are dropped here.
"""

from scipy.optimize import brentq

from app.core.rtd_engine import integrate_rtd


class BisectionError(RuntimeError):
    """Raised when the measured target isn't bracketed by the current solver
    bracket (i.e. f(a) and f(b) don't have opposite signs) -- the bracket
    needs widening in the advanced-settings panel."""


def find_k(R: float, L: float, U0: float, Eavg: float,
           k_min: float = 0.001, k_max: float = 50.0) -> float:
    """Back-solve for k such that integrate_rtd(L, U0, k, Eavg) == R."""
    def f(k: float) -> float:
        return R - integrate_rtd(L, U0, k, Eavg)

    try:
        return brentq(f, k_min, k_max, xtol=1e-10, maxiter=100)
    except ValueError as exc:
        raise BisectionError(
            f"No sign change for k in [{k_min}, {k_max}]: measured R={R} is not "
            "reachable with this bracket. Widen the k range in advanced settings."
        ) from exc


def find_e(R: float, L: float, U0: float, k: float,
           eavg_min: float = 0.0001, eavg_max: float = 3.0) -> float:
    """Back-solve for Eavg such that integrate_rtd(L, U0, k, Eavg) == R."""
    def f(Eavg: float) -> float:
        return R - integrate_rtd(L, U0, k, Eavg)

    try:
        return brentq(f, eavg_min, eavg_max, xtol=1e-10, maxiter=100)
    except ValueError as exc:
        raise BisectionError(
            f"No sign change for Eavg in [{eavg_min}, {eavg_max}]: measured R={R} is "
            "not reachable with this bracket. Widen the Eavg range in advanced settings."
        ) from exc
