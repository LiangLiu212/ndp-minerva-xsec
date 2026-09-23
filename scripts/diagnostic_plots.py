"""Standardized diagnostic plots from a DiagnosticPayload JSON.

Consumes the experiment-agnostic payload defined in
scripts/diagnostic_schema.py and produces four core plot types plus a
machine-readable summary:

  1. Efficiency vs kinematic variable (per stage)
  2. Signal purity vs kinematic variable (per stage)
  3. Data/MC ratio with stacked MC + ratio panel (per stage)
  4. Cumulative cutflow bar chart (across all stages)

No ROOT dependency — works entirely from the JSON payload and numpy +
matplotlib.

Usage:
    python scripts/diagnostic_plots.py payload.json [--output-dir DIR]
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]

# Import siblings without mutating sys.path: relative when imported as part of
# the package, plain when run directly as a script (where sys.path[0] is this
# directory). A sys.path.insert here would make the module unusable from a
# project laid out differently.
try:
    from .diagnostic_schema import DiagnosticPayload, load_payload
    from . import plot_style as ps
except ImportError:  # pragma: no cover - direct-script execution path
    from diagnostic_schema import DiagnosticPayload, load_payload
    import plot_style as ps

# All colors and markers come from plot_style (which reads the active style
# preferences YAML) — never as literals here. A one-off color for a single
# figure goes through the `override=` argument on a plot_style accessor; a
# persistent change goes through `/ndp-plot-prefs`.


def _efficiency(numerator, truth_counts):
    """Per-bin efficiency = numerator / truth_counts."""
    num = np.asarray(numerator, dtype=float)
    truth = np.asarray(truth_counts, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        eff = np.where(truth > 0, num / truth, np.nan)
    return eff


def _eff_numerator(stage):
    """The efficiency numerator for a stage.

    Schema v1.2+ carries `eff_numerator` (truth-binned: truth-signal-in-
    phase-space AND reco-selected) — the correct, <=1 numerator. Legacy
    1.0/1.1 payloads have no such field, so fall back to the reco-binned
    `mc_signal`; that is the original (buggy) behavior, but the eff>=1
    sanity flag still fires on it, surfacing that the payload predates the
    fix rather than silently reporting a wrong-but-plausible number.
    """
    if stage.eff_numerator is not None:
        return stage.eff_numerator
    return stage.mc_signal


def _efficiency_error(eff, truth_counts):
    """Binomial standard error: sqrt(ε(1-ε)/N)."""
    truth = np.asarray(truth_counts, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        err = np.where(
            truth > 0,
            np.sqrt(np.clip(eff, 0, 1) * np.clip(1 - eff, 0, 1) / truth),
            np.nan,
        )
    return err


def _purity(mc_signal, mc_background):
    """Per-bin signal purity = signal / (signal + background)."""
    sig = np.asarray(mc_signal, dtype=float)
    bkg = np.asarray(mc_background, dtype=float)
    total = sig + bkg
    with np.errstate(divide="ignore", invalid="ignore"):
        pur = np.where(total > 0, sig / total, np.nan)
    return pur


def _pot_scale(payload):
    """MC->data normalization factor from payload metadata POT.

    Experiment-agnostic: pot_data/pot_mc are filled per-analysis by the
    emit step (e.g. summed from each file's POT bookkeeping), so no
    experiment-specific constant lives in this tool.
    """
    return payload.metadata.pot_data / payload.metadata.pot_mc


def _cut_variable_ratio(payload, cutvar):
    """Per-bin Data / (POT-scaled MC total) for one cut variable.

    Computed independently of plotting so the summary carries the ratio
    even when the shared PDF was skipped (emit_shared dedup). NaN where
    the MC total is zero — matches the ratio panel's gap convention.
    """
    pot_scale = _pot_scale(payload)
    mc_total = (np.asarray(cutvar.mc_signal, dtype=float)
                + np.asarray(cutvar.mc_background, dtype=float)) * pot_scale
    data = np.asarray(cutvar.data, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(mc_total > 0, data / mc_total, np.nan)


def _pearson_chi2(data, mc_total):
    """Pearson χ²/ndf between data and MC (POT-scaled)."""
    d = np.asarray(data, dtype=float)
    m = np.asarray(mc_total, dtype=float)
    mask = m > 0
    ndf = int(mask.sum()) - 1
    if ndf <= 0:
        return float("nan"), 0
    chi2 = float(np.sum((d[mask] - m[mask]) ** 2 / m[mask]))
    return chi2, ndf


def plot_efficiency(payload, stage, centres, widths, output_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    eff = _efficiency(_eff_numerator(stage), payload.truth.counts)
    err = _efficiency_error(eff, payload.truth.counts)

    fig, ax = plt.subplots(figsize=(8, 5))
    valid = ~np.isnan(eff)
    ax.errorbar(centres[valid], eff[valid], yerr=err[valid],
                marker=ps.marker_for("efficiency"), linestyle="none",
                color=ps.color_for("efficiency"),
                markersize=ps.style("efficiency_markersize"),
                capsize=ps.style("efficiency_capsize"))
    ax.axhline(1, color=ps.color_for("reference"),
               ls=ps.style("guide_linestyle"), lw=ps.style("guide_linewidth"))
    ax.set_xlabel(f"{payload.metadata.variable_label} [{payload.metadata.units}]")
    ax.set_ylabel("Efficiency (reco selected / truth)")
    ax.set_title(f"Efficiency — {payload.metadata.channel}, "
                 f"stage: {stage.name}")
    ax.set_ylim(0, min(1.3, np.nanmax(eff) * 1.2) if np.any(valid) else 1.3)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    name = f"efficiency_{stage.name}_{payload.metadata.variable}.pdf"
    fig.savefig(str(output_dir / name), bbox_inches="tight")
    plt.close()
    return eff


def plot_purity(payload, stage, centres, widths, output_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pur = _purity(stage.mc_signal, stage.mc_background)

    fig, ax = plt.subplots(figsize=(8, 5))
    valid = ~np.isnan(pur)
    ps.draw_hist(ax, centres[valid], widths[valid], pur[valid],
                 color=ps.color_for("purity"), alpha=0.7, width_scale=0.85)
    ax.axhline(1, color=ps.color_for("reference"),
               ls=ps.style("guide_linestyle"), lw=ps.style("guide_linewidth"))
    ax.set_xlabel(f"{payload.metadata.variable_label} [{payload.metadata.units}]")
    ax.set_ylabel("Signal purity")
    ax.set_title(f"Signal purity — {payload.metadata.channel}, "
                 f"stage: {stage.name}")
    y_lo = max(0, np.nanmin(pur) - 0.05) if np.any(valid) else 0
    ax.set_ylim(y_lo, 1.02)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    name = f"purity_{stage.name}_{payload.metadata.variable}.pdf"
    fig.savefig(str(output_dir / name), bbox_inches="tight")
    plt.close()
    return pur


def plot_data_mc(payload, stage, centres, widths, output_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    try:
        from .comparison_plot import data_mc_figure
    except ImportError:  # pragma: no cover - direct-script execution
        from comparison_plot import data_mc_figure

    pot_scale = _pot_scale(payload)
    mc_sig_scaled = np.asarray(stage.mc_signal, dtype=float) * pot_scale
    mc_bkg_scaled = np.asarray(stage.mc_background, dtype=float) * pot_scale
    mc_total = mc_sig_scaled + mc_bkg_scaled
    data = np.asarray(stage.data, dtype=float)

    # The mandatory lower ratio panel comes from the shared comparison-plot
    # builder (the shared ratio-panel convention), so the layout/ratio math lives in one
    # place. Stacked MC is drawn total then background-on-top; mc_total is the
    # explicit ratio denominator. The first series is the signal+background
    # total drawn in the mc_signal color, so its label says "MC total" — in
    # step (outline) mode the orange curve traces the total, not the signal
    # alone (the signal is the orange-minus-green gap).
    fig, ax1, ax2, info = data_mc_figure(
        centres, widths, data=data,
        mc_stack=[
            {"heights": mc_sig_scaled + mc_bkg_scaled, "role": "mc_signal",
             "alpha": 0.6, "label": "MC total (signal+bkg, POT-scaled)"},
            {"heights": mc_bkg_scaled, "role": "mc_bkg",
             "alpha": 0.7, "label": "MC bkg (POT-scaled)"},
        ],
        mc_total=mc_total,
        title=f"Data vs MC — {payload.metadata.channel}, "
              f"stage: {stage.name}",
    )
    ax2.set_xlabel(f"{payload.metadata.variable_label} "
                   f"[{payload.metadata.units}]")

    plt.tight_layout()
    name = f"data_mc_{stage.name}_{payload.metadata.variable}.pdf"
    fig.savefig(str(output_dir / name), bbox_inches="tight")
    plt.close()

    ratio = info["ratio"]
    chi2, ndf = _pearson_chi2(data, mc_total)
    return ratio, chi2, ndf


def plot_cut_variable_distribution(payload, cutvar, output_dir):
    """Data vs POT-scaled MC distribution for one cut variable.

    Observable-independent (vertex Z, apothem, BDT score, boolean
    pass/fail, ...), so written once to the shared output dir. MC is
    scaled to data via POT (pot_data/pot_mc); data shown as points with
    Poisson error bars.

    This is a data-vs-MC comparison, so it carries the mandatory lower
    Data/MC ratio panel from the shared comparison-plot builder (the
    §8 convention) — same builder the per-observable data/MC plots use, so
    the layout/ratio math lives in one place. Returns the per-bin Data/MC
    ratio (NaN where MC total is zero) for the caller's summary JSON.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    try:
        from .comparison_plot import data_mc_figure
    except ImportError:  # pragma: no cover - direct-script execution
        from comparison_plot import data_mc_figure

    centres = np.array(cutvar.centres)
    widths = np.array(cutvar.widths)

    pot_scale = _pot_scale(payload)
    mc_sig = np.asarray(cutvar.mc_signal, dtype=float) * pot_scale
    mc_bkg = np.asarray(cutvar.mc_background, dtype=float) * pot_scale
    mc_total = mc_sig + mc_bkg
    data = np.asarray(cutvar.data, dtype=float)

    # Stacked MC drawn total then background-on-top; mc_total is the explicit
    # ratio denominator (a stacked draw is not its own sum). The first series
    # is the signal+background total drawn in the mc_signal color, so its
    # label says "MC total" — in step (outline) mode the orange curve traces
    # the total, not the signal alone (the signal is the orange-minus-green gap).
    fig, ax1, ax2, info = data_mc_figure(
        centres, widths, data=data,
        mc_stack=[
            {"heights": mc_total, "role": "mc_signal",
             "alpha": 0.6, "label": "MC total (signal+bkg, POT-scaled)"},
            {"heights": mc_bkg, "role": "mc_bkg",
             "alpha": 0.7, "label": "MC bkg (POT-scaled)"},
        ],
        mc_total=mc_total,
        ylabel="Events (MC POT-scaled to data)",
        title=f"Cut variable: {cutvar.label} — "
              f"{payload.metadata.channel} (stage: {cutvar.stage})",
    )
    units = f" [{cutvar.units}]" if cutvar.units else ""
    ax2.set_xlabel(f"{cutvar.label}{units}")

    plt.tight_layout()
    name = f"dist_{cutvar.name}.pdf"
    fig.savefig(str(output_dir / name), bbox_inches="tight")
    plt.close()

    return info["ratio"]


def cutflow_counts(payload):
    """Per-stage (stage_names, n_data, n_mc) with MC POT-scaled to data.

    MC counts use the same pot_data/pot_mc factor every other plot in
    this module applies (_pot_scale). Experiment-agnostic: works for any
    analysis whose payload fills pot_data/pot_mc, so the data/MC bars are
    directly comparable. Factored out of plot_cutflow so the scaling is
    unit-testable without introspecting a rendered PDF.
    """
    pot_scale = _pot_scale(payload)
    stage_names = [s.name for s in payload.stages]
    n_data = [s.n_data for s in payload.stages]
    n_mc = [(s.n_mc_signal + s.n_mc_background) * pot_scale
            for s in payload.stages]
    return stage_names, n_data, n_mc


def plot_cutflow(payload, output_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    stage_names, n_data, n_mc = cutflow_counts(payload)

    x = np.arange(len(stage_names))
    bar_w = 0.35

    fig, ax = plt.subplots(figsize=(max(8, len(stage_names) * 1.5), 5))
    # Grouped bars at x±offset (not a kinematic histogram), so this keeps
    # ax.bar with role-sourced colors rather than the hist_style helper.
    ax.bar(x - bar_w / 2, n_data, bar_w, label="Data",
           color=ps.color_for("data"), alpha=0.7)
    ax.bar(x + bar_w / 2, n_mc, bar_w, label="MC (sig+bkg, POT-scaled)",
           color=ps.color_for("mc_signal"), alpha=0.7)
    ax.set_xlabel("Selection stage")
    ax.set_ylabel("Events (MC POT-scaled to data)")
    ax.set_title(f"Cutflow — {payload.metadata.channel} "
                 f"({payload.metadata.experiment})")
    ax.set_xticks(x)
    ax.set_xticklabels(stage_names, rotation=45, ha="right", fontsize=9)
    ax.legend()
    ax.grid(alpha=0.3, axis="y")

    plt.tight_layout()
    # Cutflow is observable-independent (integrated counts per stage),
    # so the filename carries no variable suffix — written once to the
    # shared dir, shared across pT/pll/... runs.
    fig.savefig(str(output_dir / "cutflow.pdf"), bbox_inches="tight")
    plt.close()


def _should_write(path, emit_shared):
    """Decide whether to (re)write a shared artifact.

    emit_shared: True  -> always write
                 False -> never write
                 "auto"-> write only if absent (idempotent "create once")
    """
    if emit_shared is True:
        return True
    if emit_shared is False:
        return False
    return not Path(path).exists()


# Supported expected_ranges keys. Declaring an analysis range up front turns
# a "looks plausible but is wrong" benchmark into an immediate fail. Keys are
# explicit (no magic) so a typo in the run config is caught rather than
# silently checking nothing.
_EXPECTED_RANGE_KEYS = (
    "efficiency_per_bin",           # every bin, every stage
    "integrated_efficiency",        # every stage
    "final_integrated_efficiency",  # last stage only
    "mean_purity",                  # every stage
    "final_mean_purity",            # last stage only
)


def check_expected_ranges(summary, expected_ranges):
    """Return violation strings for analysis-declared expected ranges.

    `expected_ranges` maps a benchmark key (one of _EXPECTED_RANGE_KEYS) to
    an inclusive (lo, hi) range. An unknown key raises ValueError (typo
    protection). None benchmark values (undefined for a stage) are skipped.
    """
    if not expected_ranges:
        return []
    unknown = set(expected_ranges) - set(_EXPECTED_RANGE_KEYS)
    if unknown:
        raise ValueError(
            f"unknown expected_ranges key(s) {sorted(unknown)}; "
            f"supported: {list(_EXPECTED_RANGE_KEYS)}")

    violations = []
    stages = summary["stages"]

    def _chk(value, lo, hi, label):
        if value is not None and (value < lo or value > hi):
            violations.append(
                f"{label} = {value:.4f} outside declared range [{lo}, {hi}]")

    for key, (lo, hi) in expected_ranges.items():
        if key == "efficiency_per_bin":
            for s in stages:
                for i, v in enumerate(s["efficiency_per_bin"]):
                    _chk(v, lo, hi,
                         f"efficiency bin {i} stage '{s['name']}'")
        elif key == "integrated_efficiency":
            for s in stages:
                _chk(s["integrated_efficiency"], lo, hi,
                     f"integrated_efficiency stage '{s['name']}'")
        elif key == "final_integrated_efficiency":
            if stages:
                _chk(stages[-1]["integrated_efficiency"], lo, hi,
                     f"final_integrated_efficiency stage '{stages[-1]['name']}'")
        elif key == "mean_purity":
            for s in stages:
                _chk(s["mean_purity"], lo, hi,
                     f"mean_purity stage '{s['name']}'")
        elif key == "final_mean_purity":
            if stages:
                _chk(stages[-1]["mean_purity"], lo, hi,
                     f"final_mean_purity stage '{stages[-1]['name']}'")
    return violations


def run_diagnostics(payload_path, output_dir, shared_dir=None,
                    emit_shared="auto", expected_ranges=None):
    """Produce all diagnostic plots + summary JSON from a payload file.

    Per-observable artifacts (efficiency, purity, data/MC) go to
    `output_dir`. Observable-independent artifacts (cutflow, cut-variable
    distributions) go to `shared_dir` (default: `output_dir`) and are
    written once per `emit_shared` policy ("auto" => only if absent), so
    running for pT then pll into the same shared_dir does not duplicate
    them.

    `expected_ranges` (optional): analysis-declared sanity ranges checked
    against the emitted benchmarks. See `check_expected_ranges` for the
    supported keys. Any violation — or any violation of a universal
    physical bound (efficiency / purity outside [0, 1]) — sets
    `summary["sanity_status"] = "fail"`. A `fail` run is not a quotable
    result: callers (CLI `main`, producing event loops) exit non-zero on
    it. The artifacts are still written, stamped `fail`, for debugging.

    Returns the summary dict (also written to output_dir/diagnostic_summary.json).
    """
    payload = load_payload(payload_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    shared_dir = output_dir if shared_dir is None else Path(shared_dir)
    shared_dir.mkdir(parents=True, exist_ok=True)

    centres = np.array(payload.binning.centres)
    widths = np.array(payload.binning.widths)

    summary = {
        "payload": str(payload_path),
        "metadata": {
            "experiment": payload.metadata.experiment,
            "channel": payload.metadata.channel,
            "variable": payload.metadata.variable,
        },
        "stages": [],
        "sanity_flags": [],
        # Subset of sanity_flags that violate a universal physical bound
        # (efficiency / purity outside [0, 1]) — i.e. impossibilities, not
        # merely suspicious values. These force sanity_status = "fail".
        "hard_violations": [],
        # Violations of analysis-declared expected_ranges (also force fail).
        "expected_range_violations": [],
    }

    for stage in payload.stages:
        eff = plot_efficiency(payload, stage, centres, widths, output_dir)
        pur = plot_purity(payload, stage, centres, widths, output_dir)
        ratio, chi2, ndf = plot_data_mc(payload, stage, centres, widths, output_dir)

        eff_valid = eff[~np.isnan(eff)]
        pur_valid = pur[~np.isnan(pur)]
        ratio_valid = ratio[~np.isnan(ratio)]

        stage_summary = {
            "name": stage.name,
            "n_data": stage.n_data,
            "n_mc_signal": stage.n_mc_signal,
            "n_mc_background": stage.n_mc_background,
            "efficiency_per_bin": [
                float(v) if not math.isnan(v) else None for v in eff
            ],
            "purity_per_bin": [
                float(v) if not math.isnan(v) else None for v in pur
            ],
            "data_mc_ratio_per_bin": [
                float(v) if not math.isnan(v) else None for v in ratio
            ],
            "integrated_efficiency": (
                float(sum(_eff_numerator(stage)) / sum(payload.truth.counts))
                if sum(payload.truth.counts) > 0 else None
            ),
            "mean_purity": (
                float(np.nanmean(pur_valid)) if len(pur_valid) > 0 else None
            ),
            "chi2": float(chi2) if not math.isnan(chi2) else None,
            "ndf": ndf,
        }
        summary["stages"].append(stage_summary)

        # --- Hard physical-bound violations (impossibilities -> fail) ------
        # Efficiency and purity are probabilities and cannot reach/exceed 1
        # or go below 0. A violation is an upstream bug, not a fluctuation,
        # so it goes into both sanity_flags (human-readable) and
        # hard_violations (drives sanity_status). Flag eff >= 1.0 (not just
        # > 1.0) so the boundary case is caught too.
        def _hard(msg):
            summary["sanity_flags"].append(msg)
            summary["hard_violations"].append(msg)

        for i, e in enumerate(eff):
            if not math.isnan(e) and e >= 1.0:
                _hard(f"efficiency >= 1.0 in bin {i} of stage "
                      f"'{stage.name}': {e:.4f} (unphysical — efficiency "
                      f"numerator/denominator mismatch upstream)")
        int_eff = stage_summary["integrated_efficiency"]
        if int_eff is not None and int_eff >= 1.0:
            _hard(f"integrated efficiency >= 1.0 in stage '{stage.name}': "
                  f"{int_eff:.4f} (unphysical — efficiency "
                  f"numerator/denominator mismatch upstream)")
        for i, p in enumerate(pur):
            if not math.isnan(p) and (p < 0 or p > 1):
                _hard(f"purity outside [0, 1] in bin {i} of stage "
                      f"'{stage.name}': {p:.4f} (unphysical)")

        # --- Soft anomaly (suspicious, not impossible -> advisory only) ----
        # A large data/MC ratio is not physically forbidden, so it warns but
        # does NOT fail the run.
        for i, r in enumerate(ratio):
            if not math.isnan(r) and (r < 0.1 or r > 10):
                summary["sanity_flags"].append(
                    f"data/MC ratio outside [0.1, 10] in bin {i} of "
                    f"stage '{stage.name}': {r:.4f}"
                )

    # Observable-independent artifacts -> shared_dir, written once.
    summary["shared_artifacts"] = []
    if _should_write(shared_dir / "cutflow.pdf", emit_shared):
        plot_cutflow(payload, shared_dir)
        summary["shared_artifacts"].append("cutflow.pdf")
    for cutvar in payload.cut_variables:
        name = f"dist_{cutvar.name}.pdf"
        if _should_write(shared_dir / name, emit_shared):
            plot_cut_variable_distribution(payload, cutvar, shared_dir)
            summary["shared_artifacts"].append(name)

    summary["cut_variables"] = [
        {"name": cv.name, "stage": cv.stage,
         "n_data": float(sum(cv.data)),
         "n_mc": float((sum(cv.mc_signal) + sum(cv.mc_background))
                       * _pot_scale(payload)),
         "data_mc_ratio_per_bin": [
             float(v) if not math.isnan(v) else None
             for v in _cut_variable_ratio(payload, cv)
         ]}
        for cv in payload.cut_variables
    ]

    # Cutflow cross-check (schema v1.3+): a divergence between the Python
    # companion pass's integrated stage totals and the external C++ `mycuts`
    # table means the two selection implementations disagree — a correctness
    # bug, not a fluctuation — so it escalates to a hard violation (drives
    # sanity_status = "fail"), the same severity as an unphysical efficiency.
    cc = payload.cutflow_crosscheck
    if cc is not None:
        summary["cutflow_crosscheck"] = {
            "source": cc.source, "status": cc.status,
            "rel_tol": cc.rel_tol, "comparisons": cc.comparisons,
        }
        if cc.status == "fail":
            for comp in cc.comparisons:
                if not comp.get("match", True):
                    msg = (f"cutflow cross-check mismatch at stage "
                           f"'{comp.get('stage')}' quantity "
                           f"'{comp.get('quantity')}': python="
                           f"{comp.get('python')} vs external="
                           f"{comp.get('external')} (rel_diff="
                           f"{comp.get('rel_diff')}) — Python and C++ "
                           f"selections diverge, a correctness bug")
                    summary["sanity_flags"].append(msg)
                    summary["hard_violations"].append(msg)

    # Analysis-declared range checks, then the overall verdict. Any hard
    # physical-bound violation or declared-range violation makes the run a
    # non-result: sanity_status = "fail". The artifacts are still written
    # (stamped fail) for debugging, but callers exit non-zero on it.
    summary["expected_range_violations"] = check_expected_ranges(
        summary, expected_ranges)
    summary["sanity_status"] = (
        "fail" if (summary["hard_violations"]
                   or summary["expected_range_violations"]) else "pass")

    summary_path = output_dir / "diagnostic_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))

    return summary


def main():
    parser = argparse.ArgumentParser(
        description="Produce diagnostic plots from a DiagnosticPayload JSON."
    )
    parser.add_argument("payload", help="Path to the DiagnosticPayload JSON file")
    parser.add_argument("--output-dir", default=".",
                        help="Directory for per-observable plots + summary "
                             "(default: cwd)")
    parser.add_argument("--shared-dir", default=None,
                        help="Directory for observable-independent artifacts "
                             "(cutflow, cut-variable distributions). "
                             "Default: same as --output-dir.")
    parser.add_argument("--emit-shared", default="auto",
                        choices=["auto", "always", "never"],
                        help="Whether to (re)write shared artifacts: "
                             "auto=only if absent (default), always, never.")
    args = parser.parse_args()

    emit_shared = {"auto": "auto", "always": True, "never": False}[
        args.emit_shared]
    summary = run_diagnostics(args.payload, args.output_dir,
                              shared_dir=args.shared_dir,
                              emit_shared=emit_shared)

    n_flags = len(summary["sanity_flags"])
    n_stages = len(summary["stages"])
    status = summary["sanity_status"]
    print(f"[diagnostic_plots] {n_stages} stage(s), "
          f"{n_flags} sanity flag(s), status={status.upper()}")
    if n_flags > 0:
        for flag in summary["sanity_flags"]:
            print(f"  WARNING: {flag}")
    for v in summary["expected_range_violations"]:
        print(f"  EXPECTED-RANGE VIOLATION: {v}")
    print(f"[diagnostic_plots] outputs in {args.output_dir}")
    if status == "fail":
        nh = len(summary["hard_violations"])
        nr = len(summary["expected_range_violations"])
        print(f"[diagnostic_plots] FAIL: {nh} hard physical-bound "
              f"violation(s), {nr} expected-range violation(s). "
              f"This run is not a quotable result.")
        sys.exit(1)


if __name__ == "__main__":
    main()
