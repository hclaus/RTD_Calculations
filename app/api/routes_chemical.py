from fastapi import APIRouter

from app.api.routes_microbiology import _derive
from app.core.chemical import PLANCK_CONSTANT_J_S, SPEED_OF_LIGHT_M_S
from app.core.chemical import diffusivity_in_air
from app.core.chemical import kp as compute_kp_value
from app.core.chemical import photon_energy_j, photon_flux_constant
from app.core.schemas import (
    BackSolveEavgRequest,
    BackSolveEavgResponse,
    BackSolveEavgRow,
    ChemicalDiffusivityRequest,
    ChemicalDiffusivityResponse,
    ChemicalKpRequest,
    ChemicalKpResponse,
)
from app.core.solvers import BisectionError, find_e

router = APIRouter()


@router.post("/kp", response_model=ChemicalKpResponse)
def compute_kp(req: ChemicalKpRequest) -> ChemicalKpResponse:
    """kp = sigma * phi * J, where J is the actinic photon-flux conversion
    constant derived from the photon energy at the given wavelength."""
    up = photon_energy_j(PLANCK_CONSTANT_J_S, SPEED_OF_LIGHT_M_S, req.wavelength_nm)
    j = photon_flux_constant(up)
    return ChemicalKpResponse(photon_energy_j=up, kp=compute_kp_value(req.sigma, req.phi, j))


@router.post("/diffusivity", response_model=ChemicalDiffusivityResponse)
def compute_diffusivity(req: ChemicalDiffusivityRequest) -> ChemicalDiffusivityResponse:
    """Estimate the actinometer compound's diffusivity in air from its
    molecular weight (Graham's-law-style scaling from a CBr4 reference
    point) -- feeds the Peclet number / diffusion-regime verdict, the same
    way the shared diffusivity input does in microbiological mode."""
    return ChemicalDiffusivityResponse(
        diffusivity_cm2_s=diffusivity_in_air(req.molecular_weight)
    )


@router.post("/back-solve-eavg", response_model=BackSolveEavgResponse)
def back_solve_eavg(req: BackSolveEavgRequest) -> BackSolveEavgResponse:
    """Batch back-solve: given a fixed (known) rate constant and a measured
    C/C0 for each actinometer sample, find the fluence rate (Eavg) each
    sample must have seen. A bracket failure on one value doesn't fail the
    whole batch -- that row just reports its error."""
    derived = _derive(req.reactor)
    rows = []
    for c_c0 in req.c_c0_values:
        try:
            eavg = find_e(
                c_c0,
                req.reactor.L,
                derived.U0,
                req.k,
                eavg_min=req.solver_settings.eavg_min,
                eavg_max=req.solver_settings.eavg_max,
                # inside the root-find: Da = k*Eavg*tau contains the unknown
                theta=derived.theta,
            )
            rows.append(BackSolveEavgRow(C_C0=c_c0, Eavg=eavg))
        except BisectionError as exc:
            rows.append(BackSolveEavgRow(C_C0=c_c0, error=str(exc)))
    return BackSolveEavgResponse(rows=rows)
