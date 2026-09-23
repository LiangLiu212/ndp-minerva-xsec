# scripts/comparison_plot.py

The single shared builder for two-panel comparison figures with a **mandatory
lower ratio panel**.

## Purpose

Every data-vs-MC or model-vs-MC comparison plot built with these helpers carries a lower
ratio panel beneath the main panel (shared x-axis, unity reference line +
optional band). This module writes that layout + the ratio math once, so it is
not re-implemented per script. It is the enforcement point for the
ratio-panel coverage convention (documented here; a host project may also
restate it in its own conventions file).

All colors, markers, the unity line/band, and the panel geometry come from
`scripts/plot_style.py` (→ `docs/style_preferences.yaml`); no literal colors
live here. Model-overlay colors are caller-supplied and must themselves be
sourced from `plot_style` (the caller owns the categorical assignment when
there are many models).

## Builders

Both return `(fig, ax_main, ax_ratio, info)`. The caller sets the x-label on
`ax_ratio`, any title/xlim not passed in, and saves the figure.

### `data_mc_figure(centres, widths, *, data, mc_stack, mc_total, ...)`

Events comparison. Main panel: each `mc_stack` series drawn in order via
`plot_style.draw_hist` (honoring `hist_style`), then data points. Lower panel:
`Data / MC` = `data / mc_total` as points with the unity line + ±band.

- `mc_stack`: list of `{"heights", "role", "alpha", "label"}` — `role` names a
  `plot_style` color (e.g. `"mc_signal"`, `"mc_bkg"`).
- `mc_total`: explicit ratio denominator (a stacked draw is not its own sum).
- `info`: `{"ratio", "mask"}`.

### `model_data_figure(edges, *, data, models, ...)`

Cross-section comparison. Main panel: data points + one step curve per model.
Lower panel: one `Model / Data` step curve per model with the unity line
(+ optional band).

Step curves are drawn with `ax.stairs(values, edges)`. **Not**
`ax.step(edges[:-1], values, where="post")` — that ends the curve at the last
bin's *left* edge, so the final bin is never drawn on either panel, silently
and on every figure built from this shared builder.

- `models`: list of `{"name", "values", "color", "linestyle", "linewidth",
  "alpha"}`. `color` MUST be sourced from `plot_style` by the caller.
- `scale`: multiply data + models for display (e.g. `1/1e-39`).
- `data_err`: the **published** uncertainty on the measured points. Because
  `data` is a measured quantity (a cross-section, not counts), there is no
  Poisson default — omit it and no y-error bars are drawn.
- `info`: `{"ratios": {name: array}, "mask"}`.

## Style knobs

From the `ratio_panel` block of `docs/style_preferences.yaml` (via
`plot_style.ratio_panel()`), with a per-call `panel_override=` dict:

| Key | Role |
|---|---|
| `main_height` / `panel_height` | gridspec `height_ratios` |
| `data_mc_label` / `model_data_label` | ratio-panel y-labels |
| `data_mc_ylim` / `model_data_ylim` | ratio-panel y-limits |
| `model_linewidth` | step-curve width in the ratio panel |

The unity line reuses `reference_line.ratio_linestyle/linewidth`; the band
reuses `ratio_band.lo/hi/alpha` and the `reference` palette color.

## Overlay helpers

### `draw_syst_band(ax_main, ax_ratio, edges, cv, syst_err, *, ratio_denom=None, label="Syst. unc.", color=None, band_override=None)`

Overlay a ±1σ systematic uncertainty band on existing main + ratio panels.
Accepts a per-bin σ array or a full covariance matrix (extracts √diag
automatically). Pass `ax_ratio=None` for a figure with no ratio panel;
`ratio_denom` bins that are ≤ 0 get a zero-width band rather than a division
by zero.

Color defaults to `plot_style.color_for("syst_band")`. Alphas come from the
`syst_band` block via **`plot_style.syst_band()`** — a section-scoped accessor,
not `plot_style.style("alpha")`. `style()` searches several sections first-hit-
wins, and for the key `alpha` it reaches `ratio_band`'s value (0.1) before this
band's (0.25), which is the wrong number. `band_override={...}` merges per call.

`uncertainty_summary_figure` panel (a) is a caller.

### `uncertainty_summary_figure(artifact_path, *, title_prefix="", figsize=(14, 10))`

Standalone 2×2 summary of a portable uncertainty artifact (`.npz` from the
Phase 2 serializer). Panels: (a) CV with total error band, (b) fractional
uncertainty per bin, (c) error budget bar chart (top-10 bands by diagonal
variance share), (d) total correlation matrix. Returns `(fig, axes_2x2)`.

## Consumers

- `scripts/diagnostic_plots.py::plot_data_mc` (first consumer; refactored onto
  `data_mc_figure`, pixel-parity verified).
- New plotting scripts drawing a data/MC or model/MC comparison **must** use
  these builders (or replicate the ratio panel) — see the coverage convention.
  Grandfathered scripts migrate when next touched.

## Dependencies

- Python 3.10+, numpy, matplotlib
- `scripts/plot_style.py`
