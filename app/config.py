"""App-wide defaults. Legacy solver-bracket values and the demo sweep list
come straight from the original spreadsheet (Module1.bas / Calculations!N6:N17)
so the new app's defaults match what users of the old tool are used to."""

from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    # Legacy Find_k / Find_E bisection bracket defaults (Module1.bas).
    k_min_default: float = 0.001
    k_max_default: float = 50.0
    eavg_min_default: float = 0.0001
    eavg_max_default: float = 3.0

    # Legacy demo Eavg sweep values (Calculations!N6:N17).
    default_sweep_eavg_values: list[float] = [
        0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75,
    ]

    # Default C/C0 + paired Eavg values for the batch "back-solve k from
    # measured C/C0" table -- each C/C0 was measured at its own Eavg, so the
    # two lists are paired by index (same default Eavg for all rows here,
    # since there's no legacy multi-condition default to draw from).
    default_back_solve_c_c0_values: list[float] = [0.5, 0.10, 0.05, 0.01, 0.005]
    default_back_solve_eavg_values: list[float] = [1.6, 1.6, 1.6, 1.6, 1.6]

    # Chemical (Phase 2) defaults: spreadsheet's chemical!C3/C9/C8/C13 cells
    # (wavelength matches the 222nm far-UVC context throughout the project).
    default_wavelength_nm: float = 222.0
    default_sigma: float = 1.32e-17
    default_phi: float = 1.0
    default_eprime: float = 2.8
    default_molecular_weight: float = 319.62  # CBr4, g/mol -- see app/core/chemical.py
    # CORRECTED build: the flux-weighted RTD predicts higher survival than
    # the area-weighted one, so Eavg <= 3 mW/cm^2 now bottoms out around
    # C/C0 ~ 0.85. The legacy defaults (0.85, 0.80) are unreachable within
    # the solver bracket here, so these are raised accordingly.
    default_back_solve_eavg_c_c0_values: list[float] = [0.97, 0.94, 0.91, 0.88]

    # Standard pipe-flow Reynolds-number regime thresholds.
    re_laminar_max: float = 2300.0
    re_turbulent_min: float = 4000.0

    # Peclet-number advection-dominance threshold (Bowers & Miller, 2025).
    pe_advection_threshold: float = 1e5

    model_config = {"env_prefix": "RTD_"}


settings = Settings()
