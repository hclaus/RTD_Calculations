"""Molecular-diffusion correction for the laminar-flow RTD model.

Taylor's laminar RTD assumes no radial transport between streamlines.
That holds for bioaerosols (Pe ~ 1e8) but not for gas-phase chemical
actinometry (Pe ~ 1.5e3). The reactor sits in an intermediate regime
where neither Taylor's zero-diffusion RTD nor the Taylor-Aris dispersion
model applies, so the correction is precomputed by Monte Carlo and
tabulated here.

Because first-order kinetics make conversion depend only on the RTD
(Danckwerts; Zwietering), the correction is a clean multiplicative factor
on the zero-diffusion result:

    (C/C0)_corrected = F(theta, Da) * (C/C0)_zero-diffusion, flux-weighted

with the two dimensionless groups

    theta = D * tau / a^2        (radial mixing during transit)
    Da    = k * Ebar * tau       (reaction progress during transit)
    tau   = L / Uavg = pi a^2 L / Q

Both are computed from the user's inputs; no new physical input is
required beyond the diffusivity D the app already collects.

IMPORTANT: Da contains k and Ebar, so when back-solving either of them
the correction must be evaluated INSIDE the root-find, not applied
afterwards. See solvers.py.

The table is geometry-specific (a = 1.1 cm, L = 19.7 cm); the aspect
ratio a/L is the third dimensionless group and is fixed hardware.
See "Diffusion correction.md".
"""

from __future__ import annotations

import json
import math
from functools import lru_cache
from pathlib import Path

_TABLE_PATH = Path(__file__).parent.parent / "data" / "diffusion_correction_table.json"


class DiffusionTableUnavailable(RuntimeError):
    """Raised when the correction table file is missing or malformed."""


@lru_cache(maxsize=1)
def _load_table() -> dict:
    if not _TABLE_PATH.exists():
        raise DiffusionTableUnavailable(
            f"Diffusion correction table not found at {_TABLE_PATH}. "
            "Generate it with scripts/build_correction_table.py."
        )
    with open(_TABLE_PATH, encoding="utf-8") as f:
        raw = json.load(f)

    theta = raw["theta_grid"]
    da = raw["Da_grid"]
    corr = raw["correction_vs_flux"]
    if len(corr) != len(theta) or any(len(row) != len(da) for row in corr):
        raise DiffusionTableUnavailable("Correction table shape does not match its grids.")

    return {
        "log_theta": [math.log(t) for t in theta],
        "log_da": [math.log(d) for d in da],
        # interpolate log F: keeps F positive and handles the steeper
        # high-Da region better than linear interpolation of F itself
        "log_corr": [[math.log(max(c, 1e-12)) for c in row] for row in corr],
        "theta_min": theta[0],
        "theta_max": theta[-1],
        "da_min": da[0],
        "da_max": da[-1],
        "geometry": raw.get("geometry", {}),
        "monte_carlo": raw.get("monte_carlo", {}),
    }


def table_is_available() -> bool:
    try:
        _load_table()
        return True
    except DiffusionTableUnavailable:
        return False


def table_metadata() -> dict:
    t = _load_table()
    return {
        "theta_range": [t["theta_min"], t["theta_max"]],
        "Da_range": [t["da_min"], t["da_max"]],
        "geometry": t["geometry"],
        "monte_carlo": t["monte_carlo"],
    }


def _bisect_left(seq: list[float], x: float) -> int:
    lo, hi = 0, len(seq) - 1
    while lo < hi - 1:
        mid = (lo + hi) // 2
        if seq[mid] <= x:
            lo = mid
        else:
            hi = mid
    return lo


def theta_from(D: float, tau: float, a: float) -> float:
    """theta = D*tau/a^2 -- radial diffusion during one mean transit."""
    return D * tau / a**2


def damkohler(k: float, Eavg: float, tau: float) -> float:
    """Da = k*Ebar*tau -- reaction progress during one mean transit."""
    return k * Eavg * tau


def correction_factor(theta: float, Da: float) -> float:
    """Interpolated diffusion correction F(theta, Da).

    Bilinear in (log theta, log Da) on log F. Returns 1.0 when the
    correction is physically negligible rather than extrapolating:

      * Da below the tabulated range -> essentially no conversion, so no
        dose-distribution effect can exist (e.g. O2 at 222 nm sits at
        Da ~ 1e-7). Exactly 1.0.
      * theta below the tabulated range -> the aerosol regime, where the
        zero-diffusion RTD is the correct model. Exactly 1.0.

    Outside the upper edges the value is clamped to the table edge; callers
    that care should check `correction_in_range()` and surface a warning.
    """
    if Da <= 0.0 or theta <= 0.0:
        return 1.0

    t = _load_table()
    if Da < t["da_min"] or theta < t["theta_min"]:
        return 1.0

    lt = math.log(min(theta, t["theta_max"]))
    ld = math.log(min(Da, t["da_max"]))

    i = _bisect_left(t["log_theta"], lt)
    j = _bisect_left(t["log_da"], ld)

    t0, t1 = t["log_theta"][i], t["log_theta"][i + 1]
    d0, d1 = t["log_da"][j], t["log_da"][j + 1]
    wt = 0.0 if t1 == t0 else (lt - t0) / (t1 - t0)
    wd = 0.0 if d1 == d0 else (ld - d0) / (d1 - d0)

    c = t["log_corr"]
    top = c[i][j] * (1 - wd) + c[i][j + 1] * wd
    bot = c[i + 1][j] * (1 - wd) + c[i + 1][j + 1] * wd
    return math.exp(top * (1 - wt) + bot * wt)


def correction_in_range(theta: float, Da: float) -> tuple[bool, str | None]:
    """Whether (theta, Da) lies inside the tabulated grid, plus a message
    describing any extrapolation being applied."""
    t = _load_table()
    if Da < t["da_min"]:
        return True, None          # legitimately 1.0, not extrapolation
    if theta < t["theta_min"]:
        return True, None          # aerosol regime, legitimately 1.0
    if theta > t["theta_max"]:
        return False, (
            f"theta={theta:.3g} exceeds the tabulated maximum "
            f"({t['theta_max']:.3g}); correction clamped to the table edge. "
            "Very low flow with a highly diffusive gas -- treat with caution."
        )
    if Da > t["da_max"]:
        return False, (
            f"Da={Da:.3g} exceeds the tabulated maximum ({t['da_max']:.3g}); "
            "correction clamped to the table edge. C/C0 here is below 1e-12, "
            "far past any measurable range."
        )
    return True, None
