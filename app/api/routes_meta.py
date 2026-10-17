from fastapi import APIRouter

from app.core.schemas import SolverSettings
from app.metadata.loader import load_variables

router = APIRouter()


@router.get("/variables")
def get_variables():
    return [v.model_dump() for v in load_variables()]


@router.get("/solver-settings/defaults", response_model=SolverSettings)
def get_solver_settings_defaults() -> SolverSettings:
    return SolverSettings()
