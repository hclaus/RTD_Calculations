from pydantic import BaseModel, Field, model_validator

from app.config import settings


class ReactorInputs(BaseModel):
    L: float = Field(..., gt=0, description="Irradiated zone length, cm")
    a: float = Field(..., gt=0, description="Reactor radius, cm")
    Q: float = Field(..., gt=0, description="Volumetric flow rate, cm^3/min")
    mu: float = Field(..., gt=0, description="Fluid viscosity, g/(cm*s)")
    rho: float = Field(..., gt=0, description="Fluid density, g/cm^3")
    D: float = Field(..., gt=0, description="Challenge-agent diffusivity in air, cm^2/s")


class ReactorDerived(BaseModel):
    Uavg: float
    U0: float
    Re: float
    Pe: float
    q: float
    L_over_U0: float
    flow_regime: str
    diffusion_regime: str
    # theta = D*tau/a^2 -- radial mixing during one transit. Depends only on
    # the reactor and the compound's diffusivity, never on k or Eavg, so it
    # is constant while back-solving.
    theta: float


class SolverSettings(BaseModel):
    k_min: float = Field(default_factory=lambda: settings.k_min_default, gt=0)
    k_max: float = Field(default_factory=lambda: settings.k_max_default, gt=0)
    eavg_min: float = Field(default_factory=lambda: settings.eavg_min_default, gt=0)
    eavg_max: float = Field(default_factory=lambda: settings.eavg_max_default, gt=0)

    @model_validator(mode="after")
    def check_order(self):
        if self.k_min >= self.k_max:
            raise ValueError("k_min must be < k_max")
        if self.eavg_min >= self.eavg_max:
            raise ValueError("eavg_min must be < eavg_max")
        return self


class HeadlineRequest(BaseModel):
    reactor: ReactorInputs
    k: float = Field(..., gt=0, description="Inactivation rate constant, cm^2/mJ")
    Eavg: float = Field(..., gt=0, description="Average fluence rate, mW/cm^2")
    measured_C_C0: float | None = Field(None, gt=0, lt=1)
    solver_settings: SolverSettings = Field(default_factory=SolverSettings)


class HeadlineResponse(BaseModel):
    derived: ReactorDerived
    C_C0_reactor: float
    estimated_k: float | None = None
    residual: float | None = None
    # Diffusion correction diagnostics (this build only)
    Da: float
    correction_factor: float
    C_C0_uncorrected: float
    correction_warning: str | None = None


class SweepRequest(BaseModel):
    reactor: ReactorInputs
    k: float = Field(..., gt=0)
    eavg_values: list[float] = Field(
        default_factory=lambda: list(settings.default_sweep_eavg_values)
    )


class SweepRow(BaseModel):
    Eavg: float
    C_C0_rtd: float
    log10_inactivation: float
    Davg: float
    C_C0_well_mixed: float
    pct_relative_difference: float
    correction_factor: float = 1.0


class TrendlineFit(BaseModel):
    b: float
    x: list[float]
    y: list[float]


class SweepResponse(BaseModel):
    rows: list[SweepRow]
    rtd_trendline: TrendlineFit
    well_mixed_trendline: TrendlineFit


class BackSolveKRequest(BaseModel):
    reactor: ReactorInputs
    c_c0_values: list[float] = Field(
        default_factory=lambda: list(settings.default_back_solve_c_c0_values),
        min_length=1,
    )
    eavg_values: list[float] = Field(
        default_factory=lambda: list(settings.default_back_solve_eavg_values),
        min_length=1,
        description="Average fluence rate paired 1:1 with c_c0_values (each C/C0 was measured at its own Eavg)",
    )
    solver_settings: SolverSettings = Field(default_factory=SolverSettings)

    @model_validator(mode="after")
    def check_matching_lengths(self):
        if len(self.c_c0_values) != len(self.eavg_values):
            raise ValueError(
                f"c_c0_values ({len(self.c_c0_values)} values) and eavg_values "
                f"({len(self.eavg_values)} values) must have the same length -- "
                "each C/C0 measurement needs its own paired Eavg."
            )
        return self


class BackSolveKRow(BaseModel):
    C_C0: float
    Eavg: float
    k: float | None = None
    error: str | None = None


class BackSolveKResponse(BaseModel):
    rows: list[BackSolveKRow]


class ChemicalKpRequest(BaseModel):
    sigma: float = Field(..., gt=0, description="Molar absorption cross-section, cm^2/molec")
    phi: float = Field(..., gt=0, description="Quantum yield, molec/photon")
    wavelength_nm: float = Field(..., gt=0, description="Wavelength, nm")


class ChemicalKpResponse(BaseModel):
    photon_energy_j: float
    kp: float


class ChemicalDiffusivityRequest(BaseModel):
    molecular_weight: float = Field(
        default_factory=lambda: settings.default_molecular_weight,
        gt=0,
        description="Actinometer compound molecular weight, g/mol",
    )


class ChemicalDiffusivityResponse(BaseModel):
    diffusivity_cm2_s: float


class BackSolveEavgRequest(BaseModel):
    """Mirror image of BackSolveKRequest: here the rate constant (kp) is the
    known/fixed quantity -- derived from chemistry, not back-solved -- and
    each measured C/C0 yields its own estimated Eavg (the reactor fluence
    rate an actinometer sample must have seen to produce that conversion)."""

    reactor: ReactorInputs
    k: float = Field(..., gt=0, description="Fixed rate constant (e.g. kp), cm^2/mJ")
    c_c0_values: list[float] = Field(
        default_factory=lambda: list(settings.default_back_solve_eavg_c_c0_values),
        min_length=1,
    )
    solver_settings: SolverSettings = Field(default_factory=SolverSettings)


class BackSolveEavgRow(BaseModel):
    C_C0: float
    Eavg: float | None = None
    error: str | None = None


class BackSolveEavgResponse(BaseModel):
    rows: list[BackSolveEavgRow]


class OzoneRequest(BaseModel):
    """Ozone generation is ZERO order, so unlike the other two modes there
    is no rate constant to supply or back-solve -- only the measured,
    fluence-normalised generation rate G."""

    reactor: ReactorInputs
    Eavg: float = Field(..., gt=0, description="Average fluence rate, mW/cm^2")
    G: float = Field(
        default_factory=lambda: settings.default_ozone_generation_rate,
        gt=0,
        description="Measured ozone generation rate, ppb/hr per (uW/cm^2)",
    )
    kd: float = Field(
        default_factory=lambda: settings.default_ozone_kd,
        ge=0,
        description="Optional photolytic ozone-loss constant, cm^2/mJ. 0 disables it.",
    )


class OzoneResponse(BaseModel):
    derived: ReactorDerived
    ozone_ppb: float
    # G*Eavg*tau -- the pure zero-order value, shown so the effect of any
    # loss term is visible rather than silently folded in.
    ozone_ppb_zero_order: float
    ozone_ppb_well_mixed: float
    G_si: float = Field(..., description="G converted to ppb/s per (mW/cm^2)")
    # Loss diagnostics. Meaningful only when kd > 0; at kd == 0 the RTD
    # drops out of the problem entirely.
    Da_loss: float
    # What the diffusion table WOULD return at this (theta, Da). Reported
    # for transparency but NOT applied -- see ozone.ozone_outlet for why
    # (the correction multiplies S, ozone depends on 1-S, ~45x amplification).
    correction_factor: float
    correction_applied: bool = False
    correction_warning: str | None = None


class OzoneSweepRequest(BaseModel):
    reactor: ReactorInputs
    G: float = Field(
        default_factory=lambda: settings.default_ozone_generation_rate, gt=0
    )
    kd: float = Field(default_factory=lambda: settings.default_ozone_kd, ge=0)
    eavg_values: list[float] = Field(
        default_factory=lambda: list(settings.default_sweep_eavg_values)
    )


class OzoneSweepRow(BaseModel):
    Eavg: float
    ozone_ppb: float
    ozone_ppb_well_mixed: float
    ozone_ppb_zero_order: float
    correction_factor: float = 1.0


class OzoneSweepResponse(BaseModel):
    """No TrendlineFit here, unlike SweepResponse: at kd == 0 the Eavg ->
    [O3] relation is exactly linear through the origin, so the forced-
    EXPONENTIAL fit used for survival curves would be meaningless."""

    rows: list[OzoneSweepRow]
