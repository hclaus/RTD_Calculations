import math

from fastapi import APIRouter, HTTPException

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
from app.core.rtd_engine import integrate_rtd
from app.core.schemas import (
    BackSolveKRequest,
    BackSolveKResponse,
    BackSolveKRow,
    HeadlineRequest,
    HeadlineResponse,
    ReactorDerived,
    SweepRequest,
    SweepResponse,
)
from app.core.solvers import BisectionError, find_k
from app.core.sweep import build_sweep_table, fit_forced_exponential

router = APIRouter()


def _derive(reactor) -> ReactorDerived:
    Uavg = uavg(reactor.Q, reactor.a)
    U0 = u0(reactor.Q, reactor.a)
    Re = reynolds(reactor.rho, Uavg, reactor.a, reactor.mu)
    Pe = peclet(Uavg, reactor.L, reactor.D)
    return ReactorDerived(
        Uavg=Uavg,
        U0=U0,
        Re=Re,
        Pe=Pe,
        q=mean_residence_time(reactor.L, Uavg),
        L_over_U0=min_residence_time(reactor.L, U0),
        flow_regime=flow_regime(Re),
        diffusion_regime=diffusion_regime(Pe),
    )


@router.post("/headline", response_model=HeadlineResponse)
def headline(req: HeadlineRequest) -> HeadlineResponse:
    derived = _derive(req.reactor)
    c_c0 = integrate_rtd(req.reactor.L, derived.U0, req.k, req.Eavg)

    estimated_k = None
    residual = None
    if req.measured_C_C0 is not None:
        try:
            estimated_k = find_k(
                req.measured_C_C0,
                req.reactor.L,
                derived.U0,
                req.Eavg,
                k_min=req.solver_settings.k_min,
                k_max=req.solver_settings.k_max,
            )
        except BisectionError as exc:
            raise HTTPException(status_code=422, detail=str(exc)) from exc
        residual = integrate_rtd(req.reactor.L, derived.U0, estimated_k, req.Eavg) - req.measured_C_C0

    return HeadlineResponse(
        derived=derived,
        C_C0_reactor=c_c0,
        estimated_k=estimated_k,
        residual=residual,
    )


@router.post("/sweep", response_model=SweepResponse)
def sweep(req: SweepRequest) -> SweepResponse:
    derived = _derive(req.reactor)
    rows = build_sweep_table(req.reactor.L, derived.U0, derived.q, req.k, req.eavg_values)

    x = [r["Eavg"] for r in rows]
    rtd_b = fit_forced_exponential(x, [r["C_C0_rtd"] for r in rows])
    wm_b = fit_forced_exponential(x, [r["C_C0_well_mixed"] for r in rows])

    x_min, x_max = min(x), max(x)
    x_curve = [x_min + i * (x_max - x_min) / 49 for i in range(50)]

    return SweepResponse(
        rows=rows,
        rtd_trendline={"b": rtd_b, "x": x_curve, "y": [math.exp(rtd_b * xi) for xi in x_curve]},
        well_mixed_trendline={"b": wm_b, "x": x_curve, "y": [math.exp(wm_b * xi) for xi in x_curve]},
    )


@router.post("/back-solve-k", response_model=BackSolveKResponse)
def back_solve_k(req: BackSolveKRequest) -> BackSolveKResponse:
    """Batch back-solve: for each (measured C/C0, Eavg) pair -- each C/C0 was
    measured at its own fluence rate -- find the rate constant k consistent
    with it. A bracket failure on one pair doesn't fail the whole batch --
    that row just reports its error."""
    derived = _derive(req.reactor)
    rows = []
    for c_c0, eavg in zip(req.c_c0_values, req.eavg_values):
        try:
            k = find_k(
                c_c0,
                req.reactor.L,
                derived.U0,
                eavg,
                k_min=req.solver_settings.k_min,
                k_max=req.solver_settings.k_max,
            )
            rows.append(BackSolveKRow(C_C0=c_c0, Eavg=eavg, k=k))
        except BisectionError as exc:
            rows.append(BackSolveKRow(C_C0=c_c0, Eavg=eavg, error=str(exc)))
    return BackSolveKResponse(rows=rows)
