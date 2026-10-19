"""Ozone-generation endpoints.

Chemical mode reuses /api/microbiology/headline by passing kp as k, because
its output quantity is the same survival ratio. Ozone cannot do that: the
output is a concentration in ppb, and the kinetics are zero- not
first-order. Hence its own routes -- but still built on the shared
_derive() helper and the same RTD engine.
"""

from fastapi import APIRouter

from app.api.routes_microbiology import _derive
from app.core.diffusion import correction_factor, damkohler
from app.core.ozone import (
    generation_rate_si,
    ozone_outlet,
    ozone_well_mixed,
    ozone_zero_order,
)
from app.core.schemas import (
    OzoneRequest,
    OzoneResponse,
    OzoneSweepRequest,
    OzoneSweepResponse,
    OzoneSweepRow,
)

router = APIRouter()


@router.post("/headline", response_model=OzoneResponse)
def headline(req: OzoneRequest) -> OzoneResponse:
    derived = _derive(req.reactor)
    tau = derived.q
    g_si = generation_rate_si(req.G)

    # At kd == 0 these three are algebraically identical -- that identity is
    # the whole point of the mode, so all three are returned and the UI
    # plots them on top of each other.
    zero_order = ozone_zero_order(g_si, req.Eavg, tau)
    value = ozone_outlet(g_si, req.Eavg, tau, req.reactor.L, derived.U0, kd=req.kd)
    well_mixed = ozone_well_mixed(g_si, req.Eavg, tau, kd=req.kd)

    # Only meaningful once kd > 0 re-introduces first-order kinetics. The
    # factor is reported but NOT applied -- ozone depends on (1 - S), which
    # amplifies any error in the correction by ~1/Da. See core/ozone.py.
    da_loss = damkohler(req.kd, req.Eavg, tau) if req.kd > 0 else 0.0
    if req.kd > 0:
        factor = correction_factor(derived.theta, da_loss)
        warning = (
            "Diffusion correction not applied to the photolytic-loss term: "
            "ozone depends on (1 - survival), which amplifies any error in "
            f"the tabulated factor by ~{1 / da_loss:.0f}x at this Damkohler "
            "number. The flux-weighted RTD is used uncorrected."
        )
    else:
        factor, warning = 1.0, None

    return OzoneResponse(
        derived=derived,
        ozone_ppb=value,
        ozone_ppb_zero_order=zero_order,
        ozone_ppb_well_mixed=well_mixed,
        G_si=g_si,
        Da_loss=da_loss,
        correction_factor=factor,
        correction_applied=False,
        correction_warning=warning,
    )


@router.post("/sweep", response_model=OzoneSweepResponse)
def sweep(req: OzoneSweepRequest) -> OzoneSweepResponse:
    derived = _derive(req.reactor)
    tau = derived.q
    g_si = generation_rate_si(req.G)

    rows = []
    for eavg in req.eavg_values:
        da_loss = damkohler(req.kd, eavg, tau) if req.kd > 0 else 0.0
        rows.append(
            OzoneSweepRow(
                Eavg=eavg,
                ozone_ppb=ozone_outlet(
                    g_si, eavg, tau, req.reactor.L, derived.U0, kd=req.kd
                ),
                ozone_ppb_well_mixed=ozone_well_mixed(g_si, eavg, tau, kd=req.kd),
                ozone_ppb_zero_order=ozone_zero_order(g_si, eavg, tau),
                # Reported per row (Da varies with Eavg, as in sweep.py) but
                # informational only -- not applied. See core/ozone.py.
                correction_factor=(
                    correction_factor(derived.theta, da_loss) if req.kd > 0 else 1.0
                ),
            )
        )
    return OzoneSweepResponse(rows=rows)
