"""Phase 2: chemical/actinometry calculations.

Derivation confirmed directly against the live 'chemical' sheet and against
slide 15 of Theoretical Foundation.pptx (homogeneous gas-phase photolysis
kinetics): kp = sigma * phi * J, where J is a photon-flux conversion
constant derived from the photon energy at a given wavelength.

kp plugs into the exact same rtd_engine.integrate_rtd / solvers.find_k /
solvers.find_e functions used by Phase 1, in place of the microbial rate
constant k -- the reactor geometry and RTD engine are already fully generic
(see rtd_engine.py, reactor.py), so no new calculation engine is needed,
only this kp derivation and the API wiring in routes_chemical.py.

Note on sigma (molar absorption cross-section): three different literature
values for CBr4 were found across sources during analysis (slide 10's
table: 1.28e-17 cm^2/molec; slide 15's worked example: 2.85e-18 cm^2/molec;
the spreadsheet's chemical sheet: 1.32e-17 cm^2/molec). The spreadsheet's
own value (1.32e-17) is used as the UI default for parity with the legacy
tool; it remains a user-editable input, not a hardcoded constant, so this
is not a silent pick -- just the most faithful starting point.
"""

import math

# Universal physical constants (spreadsheet's chemical!C3/C4 cells).
PLANCK_CONSTANT_J_S = 6.626e-34
SPEED_OF_LIGHT_M_S = 2.998e8


def photon_energy_j(h: float, c: float, wavelength_nm: float) -> float:
    """Energy per photon, J, from Planck's constant h [J*s], speed of light
    c [m/s], and wavelength [nm]."""
    return h * c * 1e9 / wavelength_nm


def photon_flux_constant(photon_energy_j_: float) -> float:
    """J = 1/(1000*Up): converts a fluence rate in mW/cm^2 to an actinic
    photon flux in photons/(cm^2*s), matching the spreadsheet's chemical!C6."""
    return 1.0 / (1000.0 * photon_energy_j_)


def kp(sigma: float, phi: float, J: float) -> float:
    """Photochemical pseudo-first-order rate constant, cm^2/mJ.
    sigma: molar absorption cross-section, cm^2/molec.
    phi: quantum yield, molec/photon.
    J: photon-flux conversion constant from photon_flux_constant().
    """
    return sigma * phi * J


# Reference data for diffusivity_in_air(): the 4 gas-phase compounds from
# the "Properties of Challenge Agents" table on slide 10 of
# Theoretical Foundation.pptx (molecular weight, g/mol -> diffusivity in
# air, cm^2/s). Used to fit a simple power-law D = a * MW^b below.
_REFERENCE_MW_D_PAIRS = [
    (106.12, 0.073),  # Benzaldehyde
    (120.19, 0.070),  # 1,2,4-Trimethylbenzene
    (104.15, 0.071),  # Styrene
    (319.62, 0.068),  # Carbon Tetrabromide (CBr4)
]


def _fit_diffusivity_power_law(pairs: list[tuple[float, float]]) -> tuple[float, float]:
    """Least-squares fit of D = a * MW^b (linear in log-log space) to the
    reference (MW, D) pairs. Returns (a, b).

    Note on why this (and not Graham's law / pure 1/sqrt(MW) scaling) is
    used: for a trace compound diffusing through air, the Fuller-Schettler-
    Giddings correlation's mass term uses a *reduced* molar mass
    2/(1/M_compound + 1/M_air) -- since air's molar mass (~29 g/mol) is
    much smaller than these compounds' (104-320 g/mol), that reduced mass
    is only weakly sensitive to the compound's own weight, which is exactly
    why the slide 10 data shows D varying only ~0.068-0.073 cm^2/s despite
    MW spanning 104-320 g/mol -- a far flatter relationship than 1/sqrt(MW)
    predicts (which overshoots by >60% at the light end of this range).
    A direct power-law fit to the real reference points tracks this
    observed flatness far more accurately (<=2% error on all 4 points)
    than assuming Graham's law would.
    """
    n = len(pairs)
    xs = [math.log(mw) for mw, _ in pairs]
    ys = [math.log(d) for _, d in pairs]
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    b = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys)) / sum((x - mean_x) ** 2 for x in xs)
    log_a = mean_y - b * mean_x
    return math.exp(log_a), b


_DIFFUSIVITY_FIT_A, _DIFFUSIVITY_FIT_B = _fit_diffusivity_power_law(_REFERENCE_MW_D_PAIRS)


def diffusivity_in_air(molecular_weight_g_mol: float) -> float:
    """Estimate a gas-phase actinometer compound's diffusivity in air from
    its molecular weight alone, via the power-law fit to the slide 10
    reference compounds (see _fit_diffusivity_power_law). This needs only
    the one input the user actually has on hand (molecular weight), unlike
    the full Fuller-Schettler-Giddings correlation (which needs molecular
    diffusion volumes). Treat this as an estimate suitable for computing
    the Peclet number's advection-vs-diffusion verdict, not a replacement
    for a literature diffusivity value when one is available.
    """
    return _DIFFUSIVITY_FIT_A * molecular_weight_g_mol**_DIFFUSIVITY_FIT_B
