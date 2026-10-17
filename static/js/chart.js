// Renders the RTD vs. well-mixed survival chart (Eavg on X, log-scale C/C0 on
// Y), including both forced-intercept exponential trendlines. Plots the FULL
// swept range, unlike the legacy spreadsheet's hardcoded first-5-of-12-rows bug.
function renderSweepChart(sweepResponse) {
  const rows = sweepResponse.rows;
  const rtdTrace = {
    x: rows.map((r) => r.Eavg),
    y: rows.map((r) => r.C_C0_rtd),
    mode: "markers",
    type: "scatter",
    name: "C/C0 (RTD)",
    marker: { color: "#2563eb", size: 7 },
  };
  const wellMixedTrace = {
    x: rows.map((r) => r.Eavg),
    y: rows.map((r) => r.C_C0_well_mixed),
    mode: "markers",
    type: "scatter",
    name: "C/C0 (well-mixed)",
    marker: { color: "#b45309", size: 7 },
  };
  const rtdTrend = {
    x: sweepResponse.rtd_trendline.x,
    y: sweepResponse.rtd_trendline.y,
    mode: "lines",
    type: "scatter",
    name: `RTD trend (b=${sweepResponse.rtd_trendline.b.toFixed(3)})`,
    line: { color: "#2563eb", dash: "dot", width: 1.5 },
  };
  const wellMixedTrend = {
    x: sweepResponse.well_mixed_trendline.x,
    y: sweepResponse.well_mixed_trendline.y,
    mode: "lines",
    type: "scatter",
    name: `Well-mixed trend (b=${sweepResponse.well_mixed_trendline.b.toFixed(3)})`,
    line: { color: "#b45309", dash: "dot", width: 1.5 },
  };

  const layout = {
    margin: { t: 10, r: 10, b: 60, l: 70 },
    xaxis: { title: { text: "Eavg (mW/cm²)" } },
    yaxis: { title: { text: "C/C0" }, type: "log" },
    legend: { orientation: "h", y: -0.25 },
    autosize: true,
  };

  Plotly.newPlot(
    "chart-container",
    [rtdTrace, wellMixedTrace, rtdTrend, wellMixedTrend],
    layout,
    { responsive: true, displaylogo: false }
  );
}
