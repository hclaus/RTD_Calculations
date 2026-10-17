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
