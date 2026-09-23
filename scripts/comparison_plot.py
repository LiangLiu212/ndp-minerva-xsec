"""Two-panel comparison figures with a mandatory lower ratio panel.

Every data-vs-MC or model-vs-MC comparison plot built here carries a lower
ratio panel beneath the main panel (shared x-axis, unity reference line +
band). This module is the single shared builder for that convention, so the
layout/ratio math is written once instead of copied per script. See the
ratio-panel coverage convention documented in `comparison_plot.md`.

All colors, markers, the unity line, the band, and the panel geometry come
from `scripts/plot_style.py` (which reads docs/style_preferences.yaml) — no
literal colors here. Model-overlay colors are caller-supplied and must
themselves be sourced from plot_style (the caller owns the categorical
assignment when there are many models).

Builders, all returning ``(fig, ax_main, ax_ratio, info)``:

  - ``data_mc_figure``   — events: stacked MC (via plot_style.draw_hist) + data
    points; ratio panel = Data / MC (points).
  - ``model_data_figure`` — cross-section: data points + N model step curves;
    ratio panel = each Model / Data (step curves).

Overlay helpers (draw on existing axes):

  - ``draw_syst_band``  — add a systematic uncertainty band (±1σ) to main +
    ratio panels from a per-bin uncertainty array or a covariance matrix.
  - ``uncertainty_summary_figure`` — standalone 2×2 figure summarizing a
    portable uncertainty artifact (.npz from the Phase 2 serializer).

The caller sets the x-label on ``ax_ratio``, any title/xlim, and saves. ``info``
carries the computed ratio + validity mask for the caller's summary JSON.
"""

import numpy as np

try:
    from . import plot_style as ps
except ImportError:  # pragma: no cover - direct-script execution path
    import plot_style as ps


def _make_panels(figsize, panel_override):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rp = ps.ratio_panel(override=panel_override)
    fig, (ax_main, ax_ratio) = plt.subplots(
        2, 1, figsize=figsize,
        gridspec_kw={"height_ratios": [rp["main_height"], rp["panel_height"]]},
        sharex=True,
    )
    return fig, ax_main, ax_ratio, rp


def _unity_line_and_band(ax_ratio, x_lo, x_hi, draw_band):
    """Dashed unity line at 1.0 plus the optional ±band, both reference-colored."""
    ax_ratio.axhline(1, color=ps.color_for("reference"),
                     ls=ps.style("ratio_linestyle"), lw=ps.style("ratio_linewidth"))
    if draw_band:
        ax_ratio.fill_between([x_lo, x_hi], ps.style("lo"), ps.style("hi"),
                              color=ps.color_for("reference"), alpha=ps.style("alpha"))


def data_mc_figure(centres, widths, *, data, mc_stack, mc_total, data_err=None,
                   width_scale=0.9, figsize=(8, 7), data_label="Data",
                   ylabel="Events", title=None, legend_fontsize=9,
                   draw_band=True, panel_override=None):
    """Events comparison: stacked MC + data points, with a Data/MC ratio panel.

    centres, widths : bin centres / widths (arrays).
    data            : data counts per bin.
    mc_stack        : list of dicts drawn in order on the main panel, each
                      ``{"heights", "role", "alpha", "label"}`` — role names a
                      plot_style color (e.g. "mc_signal", "mc_bkg").
    mc_total        : the ratio denominator (e.g. signal+background), kept
                      explicit because a stacked draw is not its own sum.
    data_err        : data uncertainty (default sqrt(data)).

    Returns (fig, ax_main, ax_ratio, info). The caller sets the x-label on
    ax_ratio and saves the figure.
    """
    centres = np.asarray(centres, dtype=float)
    widths = np.asarray(widths, dtype=float)
    data = np.asarray(data, dtype=float)
    mc_total = np.asarray(mc_total, dtype=float)
    if data_err is None:
        data_err = np.sqrt(data)

    fig, ax_main, ax_ratio, rp = _make_panels(figsize, panel_override)

    for series in mc_stack:
        ps.draw_hist(ax_main, centres, widths, series["heights"],
                     color=ps.color_for(series["role"]),
                     alpha=series.get("alpha", 1.0),
                     label=series.get("label"), width_scale=width_scale)
    ax_main.errorbar(centres, data, yerr=data_err,
                     marker=ps.marker_for("data"), linestyle="none",
                     color=ps.color_for("data"),
                     markersize=ps.style("data_markersize"), label=data_label)
    ax_main.set_ylabel(ylabel)
    ax_main.legend(fontsize=legend_fontsize)
    if title:
        ax_main.set_title(title)
    ax_main.grid(alpha=0.3)

    mask = mc_total > 0
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(mask, data / mc_total, np.nan)
        ratio_err = np.where(mask, data_err / mc_total, np.nan)
    valid = mask & ~np.isnan(ratio)
    ax_ratio.errorbar(centres[valid], ratio[valid], yerr=ratio_err[valid],
                      marker=ps.marker_for("ratio"), linestyle="none",
                      color=ps.color_for("ratio"),
                      markersize=ps.style("ratio_markersize"))
    x_lo = centres[0] - widths[0] / 2
    x_hi = centres[-1] + widths[-1] / 2
    _unity_line_and_band(ax_ratio, x_lo, x_hi, draw_band)
    ax_ratio.set_ylabel(rp["data_mc_label"])
    ax_ratio.set_ylim(*rp["data_mc_ylim"])
    ax_ratio.grid(alpha=0.3)

    return fig, ax_main, ax_ratio, {"ratio": ratio, "mask": mask}


def model_data_figure(edges, *, data, models, data_err=None, data_xerr=None,
                      scale=1.0, figsize=(8, 8), ylabel=None, title=None,
                      data_label="Data", legend_fontsize=9, xlim=None,
                      draw_band=False, panel_override=None):
    """Cross-section comparison: data points + N model step curves, with a
    Model/Data ratio panel (one step curve per model).

    edges  : bin edges (len = nbins + 1); step curves are drawn with
             ``ax.stairs(values, edges)`` so the final bin is drawn.
    data   : measured values per bin (pre-scale).
    models : list of dicts, each
             ``{"name", "values", "color", "linestyle", "linewidth", "alpha"}``.
             `color` MUST be sourced from plot_style by the caller.
    scale  : multiply data + model values for display (e.g. 1/1e-39).

    Returns (fig, ax_main, ax_ratio, info) with info["ratios"][name] per model.
    """
    edges = np.asarray(edges, dtype=float)
    centres = (edges[:-1] + edges[1:]) / 2
    widths = np.diff(edges)
    data = np.asarray(data, dtype=float)
    if data_xerr is None:
        data_xerr = widths / 2
    # data here is a measured quantity (e.g. a cross-section), NOT counts, so
    # there is no sensible Poisson default — draw no y-error bars unless the
    # caller passes the published uncertainty.
    yerr = None if data_err is None else np.asarray(data_err) * scale

    fig, ax_main, ax_ratio, rp = _make_panels(figsize, panel_override)

    ax_main.errorbar(centres, data * scale, yerr=yerr,
                     xerr=data_xerr, marker=ps.marker_for("data"),
                     linestyle="none", color=ps.color_for("data"),
                     markersize=ps.style("data_markersize"),
                     label=data_label, zorder=20)
    for m in models:
        # stairs(values, edges) -- NOT step(edges[:-1], values, where="post"),
        # which ends the curve at the LAST BIN'S LEFT EDGE and so never draws
        # the final bin.
        ax_main.stairs(np.asarray(m["values"], dtype=float) * scale, edges,
                       baseline=None, color=m["color"],
                       lw=m.get("linewidth", 1.5),
                       ls=m.get("linestyle", "-"), alpha=m.get("alpha", 1.0),
                       label=m["name"])
    if ylabel:
        ax_main.set_ylabel(ylabel)
    if title:
        ax_main.set_title(title)
    ax_main.legend(fontsize=legend_fontsize)
    ax_main.set_ylim(bottom=0)
    ax_main.grid(alpha=0.3)
    if xlim:
        ax_main.set_xlim(xlim)

    mask = data > 0
    ratios = {}
    for m in models:
        vals = np.asarray(m["values"], dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(mask, vals / data, np.nan)
        ratios[m["name"]] = r
        ax_ratio.stairs(r, edges, baseline=None, color=m["color"],
                        lw=ps.ratio_panel(override=panel_override)["model_linewidth"],
                        ls=m.get("linestyle", "-"), alpha=m.get("alpha", 1.0))
    x_lo = centres[0] - widths[0] / 2
    x_hi = centres[-1] + widths[-1] / 2
    _unity_line_and_band(ax_ratio, x_lo, x_hi, draw_band)
    ax_ratio.set_ylabel(rp["model_data_label"])
    ax_ratio.set_ylim(*rp["model_data_ylim"])
    ax_ratio.grid(alpha=0.3)

    return fig, ax_main, ax_ratio, {"ratios": ratios, "mask": mask}


def draw_syst_band(ax_main, ax_ratio, edges, cv, syst_err, *,
                   ratio_denom=None, label="Syst. unc.", color=None,
                   band_override=None):
    """Overlay a systematic error band (±1σ) on main + ratio panels.

    edges     : bin edges (len = nbins + 1).
    cv        : central value per bin (what the band is drawn around).
    syst_err  : per-bin systematic uncertainty (√diag of a covariance, or
                a pre-computed σ array). Accepts a full covariance matrix
                (2D array) — extracts √diag automatically.
    ratio_denom : denominator for the ratio panel band. Default = cv.
                  Bins where it is <= 0 get a zero-width band rather than a
                  divide-by-zero.
    label     : legend label for the band.
    color     : band fill color (default: ps.color_for("syst_band")).
    band_override : per-call dict merged over the `syst_band` prefs block
                (e.g. ``{"alpha": 0.4}``).

    ``ax_ratio=None`` draws the main-panel band only, for a figure with no
    ratio panel.
    """
    edges = np.asarray(edges, dtype=float)
    cv = np.asarray(cv, dtype=float)
    syst_err = np.asarray(syst_err, dtype=float)
    if syst_err.ndim == 2:
        syst_err = np.sqrt(np.clip(np.diag(syst_err), 0, None))
    if ratio_denom is None:
        ratio_denom = cv
    ratio_denom = np.asarray(ratio_denom, dtype=float)

    if color is None:
        color = ps.color_for("syst_band")
    # Section-scoped lookup: ps.style("alpha") would return `ratio_band`'s
    # alpha (first hit wins across sections), not this band's.
    band = ps.syst_band(override=band_override)

    centres = (edges[:-1] + edges[1:]) / 2
    widths = np.diff(edges)

    ax_main.bar(centres, 2 * syst_err, bottom=cv - syst_err,
                width=widths, color=color, alpha=band["alpha"],
                edgecolor=band["edgecolor"], label=label, zorder=5)

    if ax_ratio is None:
        return

    mask = ratio_denom > 0
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = np.where(mask, syst_err / ratio_denom, 0)
    ax_ratio.bar(centres, 2 * frac, bottom=1 - frac,
                 width=widths, color=color, alpha=band["ratio_alpha"],
                 edgecolor=band["edgecolor"], zorder=5)


def uncertainty_summary_figure(artifact_path, *, title_prefix="", figsize=(14, 10)):
    """Standalone 2×2 summary of a portable uncertainty artifact (.npz).

    Panels: (a) CV with total error band, (b) fractional uncertainty,
    (c) error budget bar chart, (d) total correlation matrix.

    Returns (fig, axes_2x2).
    """
    import json as _json
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ps.apply_rcparams()

    z = np.load(artifact_path, allow_pickle=False)
    meta = _json.loads(str(z["meta_json"]))
    cv = z["cv"]
    edges = z["edges"]
    nbins = len(cv)
    centres = (edges[:-1] + edges[1:]) / 2
    widths = np.diff(edges)

    d_total = np.diag(z["cov_total"])
    sd_total = np.sqrt(np.clip(d_total, 0, None))
    frac_unc = z["frac_unc_total"]
    corr = z["correlation_total"]

    fig, axes = plt.subplots(2, 2, figsize=figsize)

    # (a) CV + total error band
    ax = axes[0, 0]
    ax.errorbar(centres, cv, xerr=widths / 2,
                marker=ps.marker_for("data"), linestyle="none",
                color=ps.color_for("data"), markersize=ps.style("data_markersize"),
                label="CV", zorder=10)
    draw_syst_band(ax, None, edges, cv, sd_total, label="Total unc.")
    ax.set_ylabel("CV")
    ax.set_title(f"{title_prefix}CV + total uncertainty band")
    ax.legend(fontsize=8)
    ax.set_xlim(edges[0], edges[-1])
    ax.grid(alpha=0.3)

    # (b) Fractional uncertainty
    ax = axes[0, 1]
    ax.stairs(frac_unc, edges, baseline=None, color=ps.color_for("data"), lw=1.5)
    ax.set_ylabel(r"$\sigma / \mathrm{CV}$")
    ax.set_title(f"{title_prefix}Fractional total uncertainty")
    ax.set_xlim(edges[0], edges[-1])
    ax.set_ylim(bottom=0)
    ax.grid(alpha=0.3)

    # (c) Error budget (top 10 bands)
    ax = axes[1, 0]
    band_names = meta.get("band_names", [])
    d_syst = np.diag(z.get("cov_syst_total", z["cov_total"]))
    syst_sum = d_syst.sum()
    budget = {}
    for name in band_names:
        safe = name.replace(" ", "_")
        key = f"cov_{safe}"
        if key in z.files:
            budget[name] = np.diag(z[key]).sum() / syst_sum * 100 if syst_sum > 0 else 0
    top = sorted(budget.items(), key=lambda x: x[1], reverse=True)[:10]
    if top:
        names_top, fracs_top = zip(*top)
        y_pos = np.arange(len(names_top))
        ax.barh(y_pos, fracs_top, color=ps.color_for("mc_signal"), alpha=0.8)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(names_top, fontsize=7)
        ax.invert_yaxis()
        ax.set_xlabel("% of total syst variance")
    ax.set_title(f"{title_prefix}Error budget (top 10)")
    ax.grid(alpha=0.3, axis="x")

    # (d) Total correlation matrix
    ax = axes[1, 1]
    im = ax.imshow(corr, origin="lower", cmap="RdBu_r", vmin=-1, vmax=1,
                   aspect="equal")
    ax.set_title(f"{title_prefix}Correlation matrix")
    ax.set_xlabel("bin")
    ax.set_ylabel("bin")
    fig.colorbar(im, ax=ax, label="correlation", shrink=0.8)

    fig.tight_layout()
    return fig, axes
