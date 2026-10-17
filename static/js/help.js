// Renders the variable glossary table on the Help page from GET /api/variables
// -- the same metadata source that drives the calculator's input tooltips, so
// the two can never drift apart.
(function () {
  const CATEGORY_LABELS = {
    shared: "Shared reactor & flow inputs",
    microbiological: "Microbiological inputs",
    chemical: "Chemical (actinometry) inputs",
    solver_settings: "Solver settings",
    derived: "Derived / computed results",
  };
  const CATEGORY_ORDER = ["shared", "microbiological", "chemical", "solver_settings", "derived"];

  async function renderGlossary() {
    let variables;
    try {
      const resp = await fetch("/api/variables");
      variables = await resp.json();
    } catch (err) {
      console.error("Failed to load variable metadata for the glossary:", err);
      return;
    }

    const tbody = document.querySelector("#glossary-table tbody");
    tbody.innerHTML = "";

    for (const category of CATEGORY_ORDER) {
      const rows = variables.filter((v) => v.category === category);
      if (rows.length === 0) continue;

      const heading = document.createElement("tr");
      heading.className = "category-heading";
      heading.innerHTML = `<td colspan="5">${CATEGORY_LABELS[category] || category}</td>`;
      tbody.appendChild(heading);

      for (const v of rows) {
        const tr = document.createElement("tr");
        tr.innerHTML = `
          <td>${v.symbol}</td>
          <td>${v.label}</td>
          <td>${v.units || "—"}</td>
          <td>${v.description}</td>
          <td><a href="#${v.methodology_section_anchor}">#${v.methodology_section_anchor}</a></td>
        `;
        tbody.appendChild(tr);
      }
    }
  }

  document.addEventListener("DOMContentLoaded", renderGlossary);
})();
