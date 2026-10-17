# RTD Reactor Calculator

A web app reimplementing the calculations from `RTD for Laminar Flow Reactor
calculations V2.xlsm` -- a laminar-flow UV reactor validation tool based on
the Residence Time Distribution (RTD) methodology described in *UV-C
Dose-Response Behavior of Airborne (Aerosolized) Challenge Agents I:
Theoretical Foundation* (Blatchley, Jongewaard, Claus, Hernandez, Shorno &
Bowers). A "Microbiological / Chemical" toggle on the calculator page
switches between the two applications of the same underlying RTD engine.

## Running it

**Easiest (Windows):** double-click `RTD_cal.bat`. It builds and starts the
container with Docker Compose and opens the app in your browser at
`http://localhost:8000`. Use `RTD_cal_stop.bat` to stop it. Requires Docker
Desktop to be running.

**Manually:**
```
docker compose up --build
```
then open `http://localhost:8000`.

**Without Docker (local Python dev):**
```
python -m venv .venv
.venv\Scripts\pip install -e ".[dev]"
.venv\Scripts\uvicorn app.main:app --reload
```

## Running the tests

```
.venv\Scripts\pytest
```

## What's in here

- `app/core/` -- the calculation engine (reactor hydraulics, RTD integration,
  bisection solvers, sweep-table builder), pure functions with no FastAPI
  dependency.
- `app/api/` -- the FastAPI routers: `microbiology` (headline/sweep/batch
  back-solve-k calculations), `chemical` (kp derivation from
  wavelength/sigma/phi, and batch back-solve-Eavg), `meta` (variable
  metadata / solver defaults). Chemical mode calls the *same*
  `/api/microbiology/headline` and `/sweep` endpoints as microbiological
  mode -- it just computes `kp` client-visibly via `/api/chemical/kp` and
  passes it in as `k`, since the underlying RTD engine never knew or cared
  which domain "k" came from.
- `app/metadata/variables.yaml` -- the single source of truth for every
  variable's label/units/tooltip/description. Both the calculator's input
  tooltips and the Help page glossary render from this one file via
  `GET /api/variables`, so they can't drift apart.
- `static/` -- plain HTML/CSS/JS frontend (no build step), charting via a
  vendored `plotly.min.js` (works offline/air-gapped).
- `tests/fixtures/` -- two regression fixture sets, see below.
- `RTD for Laminar Flow Reactor calculations V2.xlsm`, `Theoretical Foundation.pptx`
  -- the original spreadsheet and the theoretical-background slide deck these
  calculations are grounded in. Kept as reference; excluded from the Docker
  image via `.dockerignore`.

## Deliberate deviations from the legacy spreadsheet

These were found while porting the VBA and confirmed by direct inspection
of the live workbook -- documented here so nobody re-discovers them as
"bugs" later, and so results can be sensibly compared against the old tool.

1. **Numerics: scipy instead of hand-rolled VBA.** `Integrate_RTD`'s
   original fixed 10,000-point trapezoidal loop has a genuine off-by-one
   (it double-weights `t=L/U0` and silently drops `t=tMax-h`), overestimating
   the integral by roughly 4-10% in tested scenarios. `Find_k`/`Find_E`'s
   hand-rolled bisection has the same character as a plain root-find. Per
   explicit product decision, this port uses `scipy.integrate.quad`
   (adaptive quadrature) and `scipy.optimize.brentq` instead of reproducing
   the legacy numerics bug-for-bug. **Expect headline/sweep results here to
   differ from the legacy spreadsheet's cached values by that same small,
   documented amount** -- see `app/core/rtd_engine.py`'s docstring and
   `tests/fixtures/spreadsheet_known_values.json` vs.
   `tests/fixtures/scipy_reference_values.json`.
2. **`k` unit label fixed.** The spreadsheet labels the microbial rate
   constant's units as "cm²/mW"; dimensionally it must be cm²/mJ (it
   multiplies a dose, not a fluence rate). Fixed in this app's UI/metadata;
   no calculation changes.
3. **Sweep chart plots the full range.** The legacy chart was hardcoded to
   only plot 5 of the 12 sweep rows (Eavg 0.1-1.0 of 0.1-2.75). This app
   plots all requested points.
4. **One shared `k` for headline and sweep.** The legacy sheet used two
   disconnected `k` values (`O4`=1.83 for the sweep table, `B16`=0.87 for
   the headline calc) that looked like leftover demo data rather than an
   intentional feature. This app uses one `k` input for both by default; the
   sweep API still accepts its own `k` override for deliberate what-if
   sweeps.
5. **New: Peclet number + flow/diffusion regime verdicts.** Not present in
   the original spreadsheet. Added per explicit request: a Reynolds-number
   based laminar/transitional/turbulent verdict, and a Peclet-number based
   advection-vs-diffusion verdict (the latter needs a new diffusivity input,
   `D`, not in the original sheet -- defaulted to the ≥1 µm bioaerosol value
   worked out in the theoretical-foundation slide deck).
6. **Stateless, no persistence/auth.** By design for Phase 1 -- see the
   implementation plan for the fuller rationale.

## Regression fixtures

- `tests/fixtures/spreadsheet_known_values.json` -- exact cached cell
  values read directly from the live `.xlsm` (documentation/comparison only
  -- **not** a pass/fail target, since it reflects the legacy integration
  bug described above).
- `tests/fixtures/scipy_reference_values.json` -- computed with this repo's
  own scipy-based engine; this is the actual regression target for
  `test_rtd_engine.py`, `test_sweep.py`, and `test_api_microbiology.py`.

## Phase 2: chemical actinometry

Built on top of the exact same RTD engine as Phase 1, per the spreadsheet's
`chemical` sheet: `kp = σ·Φ·J`, where `J = 1/(1000·Up)` and
`Up = h·c/λ` (photon energy at the chosen wavelength). No new calculation
engine was needed -- `kp` just plugs into `integrate_rtd`/`find_k`/`find_e`
in place of the microbiological `k`. The "Calculation type" selector on the
calculator page swaps the left-column inputs (and a few headings/table
columns) between the two domains; the results panel, chart, and most of
the table machinery are shared, not duplicated.

One open item inherited from the source material: three different
literature values for CBr4's absorption cross-section (σ) were found across
sources (the pptx's slide 10 table: 1.28e-17 cm²/molec; slide 15's worked
example: 2.85e-18 cm²/molec; the spreadsheet's `chemical` sheet: 1.32e-17
cm²/molec). The app defaults to the spreadsheet's own value (1.32e-17) for
parity with the legacy tool, but it's a plain editable input, not a
hardcoded constant -- replace it with the literature value appropriate to
your actual compound and wavelength.
