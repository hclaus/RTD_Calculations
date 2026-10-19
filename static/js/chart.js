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
    // generous bottom margin so the x-axis title clears the horizontal
    // legend in the narrow right-hand column
    margin: { t: 10, r: 10, b: 95, l: 70 },
    xaxis: { title: { text: "Eavg (mW/cm²)", standoff: 12 } },
    yaxis: { title: { text: "C/C0" }, type: "log" },
    legend: { orientation: "h", y: -0.38, yanchor: "top" },
    autosize: true,
  };

  Plotly.newPlot(
    "chart-container",
    [rtdTrace, wellMixedTrace, rtdTrend, wellMixedTrend],
    layout,
    { responsive: true, displaylogo: false }
  );
}

// Ozone output vs. fluence rate. Deliberately a LINEAR y-axis with no
// exponential trendline: zero-order kinetics make [O3] = G*Ebar*tau, a
// straight line through the origin, so the forced-exponential fit used for
// survival curves would be meaningless here.
//
// The three traces are the point of the chart. At kd = 0 they lie exactly
// on top of one another -- that visible coincidence IS the result (the RTD
// cancels) -- and they only separate once a photolytic-loss constant makes
// the kinetics non-linear again. Line widths are staggered so the overlap
// reads as three coincident curves rather than one.
function renderOzoneChart(sweepResponse) {
  const rows = sweepResponse.rows;
  const x = rows.map((r) => r.Eavg);

  const traces = [
    {
      x, y: rows.map((r) => r.ozone_ppb_zero_order),
      mode: "lines", type: "scatter", name: "zero-order G·Ē·τ",
      line: { color: "#94a3b8", width: 9 },
    },
    {
      x, y: rows.map((r) => r.ozone_ppb_well_mixed),
      mode: "lines", type: "scatter", name: "well-mixed",
      line: { color: "#b45309", width: 4.5 },
    },
    {
      x, y: rows.map((r) => r.ozone_ppb),
      mode: "lines+markers", type: "scatter", name: "flux-weighted RTD",
      line: { color: "#2563eb", width: 1.8 },
      marker: { color: "#2563eb", size: 6 },
    },
  ];

  const layout = {
    margin: { t: 10, r: 10, b: 95, l: 70 },
    xaxis: { title: { text: "Eavg (mW/cm²)", standoff: 12 } },
    yaxis: { title: { text: "Outlet [O₃] (ppb)" }, rangemode: "tozero" },
    legend: { orientation: "h", y: -0.38, yanchor: "top" },
    autosize: true,
  };

  Plotly.newPlot("chart-container", traces, layout, {
    responsive: true,
    displaylogo: false,
  });
}
