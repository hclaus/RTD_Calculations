// Wires the input fields to POST /api/microbiology/headline, /sweep, and
// /back-solve-k, and renders the results panel and both tables/chart from
// the responses. Editing any input auto-recalculates everything (debounced
// 400ms), like a spreadsheet; the two "Calculate" buttons in the middle
// column remain as an explicit, non-debounced trigger.
(function () {
  // Maps each input element id to its default value (legacy spreadsheet
  // defaults for L/a/Q/mu/rho/k/Eavg/measured_C_C0; pptx slide 9 for D).
  // Chemical defaults (wavelength/sigma/phi/Eprime) mirror the legacy
  // spreadsheet's chemical!C3/C9/C8/C13 cells -- see app/config.py.
  const DEFAULTS = {
    "in-L": 19.7, "in-a": 1.1, "in-Q": 1200, "in-mu": 0.00018, "in-rho": 0.001204,
    "in-D": 2.77e-7, "in-k": 0.87, "in-Eavg": 1.6, "in-measured": 0.00202,
    "in-wavelength": 222.0, "in-sigma": 1.32e-17, "in-phi": 1.0, "in-eprime": 2.8,
    "in-molecular-weight": 319.62,
  };
  // Bioaerosol default -- restored when switching back to microbiological
  // mode, since chemical mode overwrites #in-D with a gas-phase estimate
  // (~0.07 cm^2/s, ~6 orders of magnitude larger than this) derived from
  // molecular weight.
  const MICROBIO_DEFAULT_D = 2.77e-7;
  const DEFAULT_SWEEP = [0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75];
  const DEFAULT_BACKSOLVE_CC0_MICROBIO = [0.5, 0.10, 0.05, 0.01, 0.005];
  const DEFAULT_BACKSOLVE_EAVG = [1.6, 1.6, 1.6, 1.6, 1.6];
  const DEFAULT_BACKSOLVE_CC0_CHEMICAL = [0.95, 0.90, 0.85, 0.80];

  function $(id) { return document.getElementById(id); }

  function getMode() { return $("in-calc-type").value; }

  // Fires fn 400ms after the last call -- used so edits feel "live" (like a
  // spreadsheet recalculating as you type) without hammering the backend on
  // every keystroke.
  function debounce(fn, waitMs) {
    let timer;
    return (...args) => {
      clearTimeout(timer);
      timer = setTimeout(() => fn(...args), waitMs);
    };
  }

  function fillDefaults() {
    for (const [elementId, value] of Object.entries(DEFAULTS)) {
      const el = $(elementId);
      if (el) el.value = value;
    }
    $("in-sweep-values").value = DEFAULT_SWEEP.join(", ");
    $("in-backsolve-values").value = DEFAULT_BACKSOLVE_CC0_MICROBIO.join(", ");
    $("in-backsolve-eavg-values").value = DEFAULT_BACKSOLVE_EAVG.join(", ");
  }

  // Swaps which fieldset/fields/labels are visible between Microbiological
  // and Chemical mode. Both modes reuse the same results panel, both tables,
  // the chart, and the sweep-values / C-C0-values input fields -- only the
  // mode-specific bits (rate-constant inputs, the paired-Eavg field, a few
  // headings/columns) actually differ.
  async function applyModeVisibility() {
    const mode = getMode();
    const isChemical = mode === "chemical";

    // Keep the shared diffusivity input physically appropriate for the
    // active domain: gas-phase estimate (from MW) in chemical mode,
    // bioaerosol default when back in microbiological mode.
    if (isChemical) {
      try {
        await updateDiffusivityFromMolecularWeight();
      } catch (err) {
        showError(err.message || String(err));
      }
    } else {
      $("in-D").value = MICROBIO_DEFAULT_D;
    }

    $("microbio-fieldset").classList.toggle("hidden", isChemical);
    $("chemical-fieldset").classList.toggle("hidden", !isChemical);

    // Chemical mode uses one shared kp for the whole back-solve batch, so
    // the per-row paired-Eavg field (microbio-only) is hidden.
    $("backsolve-eavg-group").classList.toggle("hidden", isChemical);
    $("backsolve-summary").textContent = isChemical
      ? "Calculate E' from measured C/C0"
      : "Calculate k from measured C/C0";
    $("backsolve-section-title").textContent = isChemical
      ? "Calculated E' (mW/cm²) from C/C0 values"
      : "Calculated k (cm²/mJ) from C/C0 and Eavg values";
    $("backsolve-table-head").innerHTML = isChemical
      ? "<tr><th>C/C0</th><th>E'</th></tr>"
      : "<tr><th>C/C0</th><th>Eavg</th><th>k</th></tr>";
    $("in-backsolve-values").value = (isChemical ? DEFAULT_BACKSOLVE_CC0_CHEMICAL : DEFAULT_BACKSOLVE_CC0_MICROBIO).join(", ");

    $("sweep-values-summary").textContent = isChemical ? "Sweep E' values" : "Sweep Eavg values";

    // The single-measurement headline back-solve (tied to microbio's
    // "Measured C/C0" field) doesn't apply in chemical mode -- kp is
    // derived from known chemistry, not something to back-solve per-point;
    // that back-solve lives entirely in the batch table above instead.
    $("result-item-est-k").classList.toggle("hidden", isChemical);
    $("result-item-residual").classList.toggle("hidden", isChemical);

    // Clear stale content from the other mode so nothing mismatched lingers.
    document.querySelector("#backsolve-table tbody").innerHTML = "";
    document.querySelector("#sweep-table tbody").innerHTML = "";
    $("sweep-heading-detail").textContent = "";
  }

  // Computes kp = sigma*phi*J from the Chemical fieldset's current inputs.
  async function getChemicalKp() {
    const body = {
      sigma: parseFloat($("in-sigma").value),
      phi: parseFloat($("in-phi").value),
      wavelength_nm: parseFloat($("in-wavelength").value),
    };
    const resp = await fetch("/api/chemical/kp", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body),
    });
    if (!resp.ok) {
      const err = await resp.json();
      throw new Error(err.detail || "Failed to compute kp from wavelength/sigma/phi.");
    }
    return resp.json();
  }

  // Estimates diffusivity in air from the Chemical fieldset's molecular
  // weight and writes it into the shared #in-D field, so the existing
  // Peclet-number computation (reactor.py, unchanged) picks it up. Only
  // called when molecular weight itself changes or chemical mode is
  // entered -- not on every chemical-field edit -- so a manual override of
  // #in-D isn't silently clobbered by e.g. editing wavelength.
  async function updateDiffusivityFromMolecularWeight() {
    const molecular_weight = parseFloat($("in-molecular-weight").value);
    const resp = await fetch("/api/chemical/diffusivity", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ molecular_weight }),
    });
    if (!resp.ok) {
      const err = await resp.json();
      throw new Error(err.detail || "Failed to estimate diffusivity from molecular weight.");
    }
    const { diffusivity_cm2_s } = await resp.json();
    $("in-D").value = diffusivity_cm2_s;
  }

  async function fillSolverDefaults() {
    try {
      const resp = await fetch("/api/solver-settings/defaults");
      const defaults = await resp.json();
      $("in-k-min").value = defaults.k_min;
      $("in-k-max").value = defaults.k_max;
      $("in-eavg-min").value = defaults.eavg_min;
      $("in-eavg-max").value = defaults.eavg_max;
    } catch (err) {
      console.error("Failed to load solver-settings defaults:", err);
    }
  }

  function readReactorInputs() {
    return {
      L: parseFloat($("in-L").value),
      a: parseFloat($("in-a").value),
      Q: parseFloat($("in-Q").value),
      mu: parseFloat($("in-mu").value),
      rho: parseFloat($("in-rho").value),
      D: parseFloat($("in-D").value),
    };
  }

  function readSolverSettings() {
    return {
      k_min: parseFloat($("in-k-min").value),
      k_max: parseFloat($("in-k-max").value),
      eavg_min: parseFloat($("in-eavg-min").value),
      eavg_max: parseFloat($("in-eavg-max").value),
    };
  }

  // Accepts commas AND/OR any run of whitespace as the separator, so values
  // typed or pasted space-, tab-, or newline-separated (e.g. a row or column
  // copied straight out of Excel) parse correctly even without going through
  // the paste-normalizer below (e.g. if the browser's paste event is bypassed).
  function readCommaList(elementId) {
    return $(elementId).value
      .split(/[,\s]+/)
      .map((s) => parseFloat(s.trim()))
      .filter((v) => !isNaN(v));
  }

  // On paste, rewrites whatever was pasted (Excel row -> tab-separated,
  // Excel column -> newline-separated, or plain space-separated) into the
  // field's displayed "1, 2, 3" comma format, so the textarea always looks
  // consistent regardless of where the numbers came from.
  function attachPasteNormalizer(elementId) {
    const el = $(elementId);
    el.addEventListener("paste", (event) => {
      const clipboard = event.clipboardData || window.clipboardData;
      if (!clipboard) return;
      const text = clipboard.getData("text");
      const values = text
        .split(/[,\s]+/)
        .map((s) => s.trim())
        .filter((s) => s !== "" && !isNaN(parseFloat(s)));
      if (values.length === 0) return; // not number-like; let the default paste happen
      event.preventDefault();
      el.value = values.join(", ");
      // Programmatic value assignment doesn't fire 'input' on its own --
      // dispatch one so the auto-recalculate listener (attached below) picks
      // up pasted values too, not just typed ones.
      el.dispatchEvent(new Event("input", { bubbles: true }));
    });
  }

  function showError(message) {
    const box = $("error-box");
    box.textContent = message;
    box.classList.remove("hidden");
  }

  function clearError() {
    $("error-box").classList.add("hidden");
  }

  function regimeBadgeClass(text) {
    if (text.startsWith("Laminar")) return "laminar";
    if (text.startsWith("Transitional")) return "transitional";
    if (text.startsWith("Turbulent")) return "turbulent";
    if (text.startsWith("Advection")) return "advection";
    return "diffusion";
  }

  function renderHeadline(data) {
    $("res-c-c0").textContent = data.C_C0_reactor.toExponential(4);
    $("res-est-k").textContent = data.estimated_k != null ? data.estimated_k.toFixed(4) : "—";
    $("res-residual").textContent = data.residual != null ? data.residual.toExponential(2) : "—";
    $("res-q").textContent = data.derived.q.toFixed(4) + " s";
    $("res-uavg").textContent = data.derived.Uavg.toFixed(3) + " cm/s";
    $("res-u0").textContent = data.derived.U0.toFixed(3) + " cm/s";

    const flowBadge = regimeBadgeClass(data.derived.flow_regime);
    $("res-flow-regime").innerHTML =
      `<span class="badge ${flowBadge}">${data.derived.flow_regime}</span> (Re=${data.derived.Re.toFixed(1)})`;

    const diffBadge = regimeBadgeClass(data.derived.diffusion_regime);
    $("res-diffusion-regime").innerHTML =
      `<span class="badge ${diffBadge}">${data.derived.diffusion_regime.split(" — ")[0]}</span> (Pe=${data.derived.Pe.toExponential(2)})`;
  }

  function renderSweepTable(sweepResponse) {
    const tbody = document.querySelector("#sweep-table tbody");
    tbody.innerHTML = "";
    for (const row of sweepResponse.rows) {
      const tr = document.createElement("tr");
      tr.innerHTML = `
        <td>${row.Eavg}</td>
        <td>${row.Davg.toFixed(3)}</td>
        <td>${row.C_C0_rtd.toExponential(3)}</td>
        <td>${row.log10_inactivation.toFixed(3)}</td>
        <td>${row.C_C0_well_mixed.toExponential(3)}</td>
        <td>${(row.pct_relative_difference * 100).toFixed(1)}%</td>
      `;
      tbody.appendChild(tr);
    }
  }

  function renderBackSolveTable(backSolveResponse, mode) {
    const tbody = document.querySelector("#backsolve-table tbody");
    tbody.innerHTML = "";
    for (const row of backSolveResponse.rows) {
      const tr = document.createElement("tr");
      if (mode === "chemical") {
        const eavgCell = row.Eavg != null
          ? row.Eavg.toFixed(4)
          : `<span class="badge turbulent" title="${row.error || ""}">unreachable</span>`;
        tr.innerHTML = `<td>${row.C_C0}</td><td>${eavgCell}</td>`;
      } else {
        const kCell = row.k != null
          ? row.k.toFixed(4)
          : `<span class="badge turbulent" title="${row.error || ""}">unreachable</span>`;
        tr.innerHTML = `<td>${row.C_C0}</td><td>${row.Eavg}</td><td>${kCell}</td>`;
      }
      tbody.appendChild(tr);
    }
  }

  async function calculateHeadlineAndSweep(autoOpen = true) {
    clearError();
    const mode = getMode();
    const reactor = readReactorInputs();
    const solver_settings = readSolverSettings();

    try {
      let k, Eavg, measured_C_C0, headingPrefix;
      if (mode === "chemical") {
        const { kp } = await getChemicalKp();
        k = kp;
        Eavg = parseFloat($("in-eprime").value);
        measured_C_C0 = null; // chemical mode's single-point back-solve isn't used -- see batch table below
        headingPrefix = `(kp=${kp.toExponential(4)} cm²/mJ)`;
      } else {
        k = parseFloat($("in-k").value);
        Eavg = parseFloat($("in-Eavg").value);
        const measuredRaw = $("in-measured").value.trim();
        measured_C_C0 = measuredRaw === "" ? null : parseFloat(measuredRaw);
        headingPrefix = `(k=${k} cm²/mJ)`;
      }

      const headlineResp = await fetch("/api/microbiology/headline", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ reactor, k, Eavg, measured_C_C0, solver_settings }),
      });
      if (!headlineResp.ok) {
        const err = await headlineResp.json();
        throw new Error(err.detail || "Headline calculation failed.");
      }
      renderHeadline(await headlineResp.json());

      const sweepResp = await fetch("/api/microbiology/sweep", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ reactor, k, eavg_values: readCommaList("in-sweep-values") }),
      });
      if (!sweepResp.ok) {
        const err = await sweepResp.json();
        throw new Error(err.detail || "Sweep calculation failed.");
      }
      const sweepData = await sweepResp.json();
      renderSweepTable(sweepData);
      renderSweepChart(sweepData);
      $("sweep-heading-detail").textContent = headingPrefix;
      if (autoOpen) $("sweep-table-details").open = true;
    } catch (err) {
      showError(err.message || String(err));
    }
  }

  async function calculateBackSolveK(autoOpen = true) {
    clearError();
    const mode = getMode();
    const reactor = readReactorInputs();
    const solver_settings = readSolverSettings();
    const c_c0_values = readCommaList("in-backsolve-values");

    try {
      if (mode === "chemical") {
        const { kp } = await getChemicalKp();
        const resp = await fetch("/api/chemical/back-solve-eavg", {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ reactor, k: kp, c_c0_values, solver_settings }),
        });
        if (!resp.ok) {
          const err = await resp.json();
          throw new Error(err.detail || "Back-solve calculation failed.");
        }
        renderBackSolveTable(await resp.json(), mode);
      } else {
        const eavg_values = readCommaList("in-backsolve-eavg-values");
        if (c_c0_values.length !== eavg_values.length) {
          showError(
            `C/C0 values (${c_c0_values.length}) and their paired Eavg values (${eavg_values.length}) ` +
            "must have the same count -- each C/C0 needs its own Eavg."
          );
          return;
        }
        const resp = await fetch("/api/microbiology/back-solve-k", {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ reactor, c_c0_values, eavg_values, solver_settings }),
        });
        if (!resp.ok) {
          const err = await resp.json();
          throw new Error(err.detail || "Back-solve calculation failed.");
        }
        renderBackSolveTable(await resp.json(), mode);
      }
      if (autoOpen) $("backsolve-table-details").open = true;
    } catch (err) {
      showError(err.message || String(err));
    }
  }

  // Auto-recalc is scoped exactly like the two "Calculate" buttons: each
  // input only triggers the table(s) whose calculation actually reads it --
  // e.g. editing the back-solve C/C0 list never touches the sweep table, and
  // vice versa. Shared reactor/chemistry inputs feed both, since both
  // calculations read them.
  const HEADLINE_SWEEP_ONLY_IDS = ["in-k", "in-Eavg", "in-measured", "in-eprime", "in-sweep-values"];
  const BACKSOLVE_ONLY_IDS = ["in-backsolve-values", "in-backsolve-eavg-values"];
  const SHARED_IDS = [
    "in-L", "in-a", "in-Q", "in-mu", "in-rho", "in-D",
    "in-wavelength", "in-sigma", "in-phi",
    "in-k-min", "in-k-max", "in-eavg-min", "in-eavg-max",
  ];

  document.addEventListener("DOMContentLoaded", async () => {
    fillDefaults();
    await fillSolverDefaults();
    await applyModeVisibility();
    attachPasteNormalizer("in-sweep-values");
    attachPasteNormalizer("in-backsolve-values");
    attachPasteNormalizer("in-backsolve-eavg-values");
    $("sweep-calc-btn").addEventListener("click", () => calculateHeadlineAndSweep(true));
    $("backsolve-calc-btn").addEventListener("click", () => calculateBackSolveK(true));
    $("in-calc-type").addEventListener("change", async () => {
      await applyModeVisibility();
      calculateHeadlineAndSweep(true);
      calculateBackSolveK(true);
    });

    // Molecular weight only affects diffusivity (not kp), so it gets its
    // own debounced handler: update #in-D from MW, then silently
    // recalculate both tables (D is shared, so both read it).
    const debouncedMolecularWeight = debounce(async () => {
      try {
        await updateDiffusivityFromMolecularWeight();
      } catch (err) {
        showError(err.message || String(err));
        return;
      }
      calculateHeadlineAndSweep(false);
      calculateBackSolveK(false);
    }, 400);
    $("in-molecular-weight").addEventListener("input", debouncedMolecularWeight);

    // Live edits recalculate silently (no forced table expansion) -- only an
    // explicit Calculate click or a mode switch is a deliberate "show me
    // the table" action.
    const debouncedHeadlineSweep = debounce(() => calculateHeadlineAndSweep(false), 400);
    const debouncedBackSolve = debounce(() => calculateBackSolveK(false), 400);

    for (const id of HEADLINE_SWEEP_ONLY_IDS) {
      const el = $(id);
      if (el) el.addEventListener("input", debouncedHeadlineSweep);
    }
    for (const id of BACKSOLVE_ONLY_IDS) {
      const el = $(id);
      if (el) el.addEventListener("input", debouncedBackSolve);
    }
    for (const id of SHARED_IDS) {
      const el = $(id);
      if (!el) continue;
      el.addEventListener("input", debouncedHeadlineSweep);
      el.addEventListener("input", debouncedBackSolve);
    }

    // Auto-run both once with the defaults so the page isn't empty on load,
    // but don't force the tables open just because of this automatic run --
    // they should start collapsed and only pop open from an actual user
    // action (a Calculate click or editing a value).
    calculateHeadlineAndSweep(false);
    calculateBackSolveK(false);
  });
})();
