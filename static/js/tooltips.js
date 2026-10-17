// Fetches /api/variables once and attaches each variable's tooltip text (via
// a native `title` attribute) to every matching `[data-var-id]` element. The
// info icon is itself a link to /help#<anchor> -- hover shows the tooltip,
// click jumps to the full explanation -- so no separate "Learn more" line is
// needed. Applies to both input field rows (label + input) and result-panel
// items (a bare `.result-label` span), whichever the element contains.
(function () {
  async function applyTooltips() {
    let variables;
    try {
      const resp = await fetch("/api/variables");
      variables = await resp.json();
    } catch (err) {
      console.error("Failed to load variable metadata for tooltips:", err);
      return;
    }

    const byId = {};
    for (const v of variables) byId[v.id] = v;

    document.querySelectorAll("[data-var-id]").forEach((el) => {
      const meta = byId[el.dataset.varId];
      if (!meta) return;

      // Input rows: a wrapping `.field-row` with a <label>. Result items:
      // the `.result-label` span itself carries `data-var-id`.
      const label = el.matches("label, .result-label") ? el : el.querySelector("label, .result-label");
      if (!label) return;

      if (label.classList.contains("result-label")) {
        label.title = meta.tooltip;
        return;
      }

      if (!label.querySelector(".info-icon")) {
        const icon = document.createElement("a");
        icon.className = "info-icon";
        icon.textContent = "i";
        icon.title = meta.tooltip;
        icon.href = "/help#" + meta.methodology_section_anchor;
        icon.target = "_blank";
        label.appendChild(icon);
      }
    });

    window.__variableMetadata = byId;
  }

  document.addEventListener("DOMContentLoaded", applyTooltips);
})();
