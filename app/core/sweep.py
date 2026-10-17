"""Builds the Eavg sweep table (RTD-weighted survival vs. naive well-mixed
survival) that drives the comparison chart -- the tool's core pedagogical
point that dose non-uniformity under laminar flow matters.

Reproduces Calculations!N6:S17 (sweep table) and the chart's two exponential
trendlines, with one deliberate fix: the legacy chart only plotted rows 6-10
(Eavg 0.1-1.0) out of 12 due to a hardcoded range; this module has no such
truncation -- callers get the full requested set of Eavg values.
"""

import math

from app.core.rtd_engine import integrate_rtd


def build_sweep_table(L: float, U0: float, q: float, k: float, eavg_values: list[float]) -> list[dict]:
    """For each Eavg, compute RTD-weighted survival, log10 inactivation, the
    naive well-mixed-dose survival, and their relative difference."""
    rows = []
    for Eavg in eavg_values:
        c_c0_rtd = integrate_rtd(L, U0, k, Eavg)
        davg = Eavg * q
        c_c0_well_mixed = math.exp(-k * davg)
        rows.append(
            {
                "Eavg": Eavg,
                "C_C0_rtd": c_c0_rtd,
                "log10_inactivation": -math.log10(c_c0_rtd),
                "Davg": davg,
                "C_C0_well_mixed": c_c0_well_mixed,
                "pct_relative_difference": (c_c0_well_mixed - c_c0_rtd) / c_c0_well_mixed,
            }
        )
    return rows


def fit_forced_exponential(x: list[float], y: list[float]) -> float:
    """OLS fit of y = exp(b*x) forced through (0,1) -- i.e. ln(y) = b*x with
    no intercept term, matching Excel's exponential trendline with
    "Set Intercept" = 1. Returns b."""
    sum_x_lny = sum(xi * math.log(yi) for xi, yi in zip(x, y))
    sum_x2 = sum(xi * xi for xi in x)
    return sum_x_lny / sum_x2
