#!/usr/bin/env python3
"""Side-by-side figures and report: the official GENIE MC and the GiBUU prediction against the same
reconstructed data, on the released muon + leading-proton grids.

    python report/make_genie_gibuu_comparison.py \
        --selection runs/<date>_selection_minerva_me_ccqelike_1mu1p_FHC \
        --overlay   runs/<date>_gibuu_..._all_3 \
        --signal    runs/<date>_signal_minerva_me_ccqelike_1mu1p_FHC \
        --out report/GENIE_vs_GiBUU_FHC.md [--tag cmp]

Run it with PYTHONPATH set to the repo root: the process-composition table loads the channel and the
cached GiBUU sample through `ndp`.

Both panels of a figure are built the same way: the identical data points, the identical MC background
stacked by category, and on top of it the signal predicted by that panel's generator — GENIE's own
selected signal (POT-scaled) on the left, GiBUU folded through the MC response on the right. Only the
signal block differs, so the panels are directly comparable by eye.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
DATA_COLOR = "#08519c"
BKG_CATS = ("bkg 1 pi+-", "bkg 1 pi0", "bkg multi-pi", "bkg other (no pion)")
BKG_COLORS = {"bkg 1 pi+-": "#ff7f0e", "bkg 1 pi0": "#ffbb78", "bkg multi-pi": "#d62728", "bkg other (no pion)": "#9e9e9e"}
SIG_COLOR = {"GENIE": "#1f77b4", "GiBUU": "#2ca25f"}
LABELS = {"muon_p": "muon p [GeV/c]", "muon_theta": "muon θ [deg]", "muon_pt": "muon p_T [GeV/c]",
          "proton_p": "leading proton p [GeV/c]", "proton_theta": "leading proton θ [deg]", "proton_pt": "leading proton p_T [GeV/c]",
          "dpt": "δp_T [GeV/c]", "dpt_fine": "δp_T (fine) [GeV/c]", "dptx": "δp_Tx [GeV/c]", "dpty": "δp_Ty [GeV/c]",
          "alpha": "δα_T [deg]", "phi": "φ_T [deg]", "pl": "δp_L [GeV/c]", "pn": "p_n [GeV/c]"}
ORDER = ["muon_p", "muon_theta", "muon_pt", "proton_p", "proton_theta", "proton_pt",
         "dpt", "dpt_fine", "dptx", "dpty", "alpha", "phi", "pl", "pn"]


def panel(ax, axr, edges, data, bkg_by_cat, signal, color, title, xlabel, ymax):
    """One generator's panel: background stacked by category, its signal on top, data points, ratio."""
    widths = np.diff(edges); centres = 0.5 * (edges[:-1] + edges[1:])
    bottom = np.zeros(len(centres))
    for c in BKG_CATS:
        v = np.asarray(bkg_by_cat.get(c, np.zeros(len(centres))), float)
        if v.sum() > 0:
            ax.bar(edges[:-1], v, width=widths, bottom=bottom, align="edge", color=BKG_COLORS[c], label=c, lw=0)
            bottom += v
    ax.bar(edges[:-1], signal, width=widths, bottom=bottom, align="edge", color=color, alpha=0.9, lw=0,
           label=f"{title.split()[0]} signal (1μ, leading p)")
    pred = bottom + signal
    ax.errorbar(centres, data, xerr=widths / 2, yerr=np.sqrt(np.maximum(data, 0)), fmt="o", ms=4,
                color=DATA_COLOR, capsize=2, label="data", zorder=5)
    ax.set_ylim(0, ymax); ax.set_ylabel("selected events / bin"); ax.legend(fontsize=7, ncol=2, loc="upper right")
    ratio = float(data.sum() / pred.sum()) if pred.sum() else float("nan")
    ax.set_title(f"{title}   data/pred = {ratio:.3f}", fontsize=10)
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where(pred > 0, data / pred, np.nan)
        err = np.where(pred > 0, np.sqrt(np.maximum(data, 0)) / pred, np.nan)
    axr.errorbar(centres, r, xerr=widths / 2, yerr=err, fmt="o", ms=4, color=DATA_COLOR, capsize=2)
    axr.axhline(1, color="k", lw=0.8, ls="--"); axr.set_ylim(0.4, 1.6)
    axr.set_ylabel("data / prediction"); axr.set_xlabel(xlabel)
    return ratio, r


def make_figure(name, genie, gibuu, out: Path) -> tuple:
    edges = np.asarray(genie["edges"], float)
    data = np.asarray(genie["data"], float)
    cats = {c: np.asarray(v, float) for c, v in genie["mc_by_category_scaled"].items() if c in BKG_CATS}
    genie_sig = np.asarray(genie["mc_scaled"], float) - sum(cats.values())
    gibuu_pred = np.asarray(gibuu["pred"], float); gibuu_bkg = np.asarray(gibuu["bkg"], float)
    gibuu_sig = gibuu_pred - gibuu_bkg
    ymax = 1.45 * max(float(np.max(np.asarray(genie["mc_scaled"], float))), float(np.max(gibuu_pred)), float(np.max(data)))
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 6.8), gridspec_kw={"height_ratios": [3, 1]}, sharex="col")
    lab = LABELS.get(name, name)
    r_g, rb_g = panel(axes[0, 0], axes[1, 0], edges, data, cats, genie_sig, SIG_COLOR["GENIE"],
                      "GENIE official MC (POT-scaled)", lab, ymax)
    r_b, rb_b = panel(axes[0, 1], axes[1, 1], edges, data, cats, gibuu_sig, SIG_COLOR["GiBUU"],
                      "GiBUU folded through the MC response", lab, ymax)
    fig.suptitle(f"MINERvA ME FHC CC μ + leading proton, {lab}: the same data and the same MC background, two signal models", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96)); out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=120); plt.close(fig)
    return (r_g, r_b, rb_g, rb_b, float(genie_sig.sum()), float(gibuu_sig.sum()), float(sum(cats.values()).sum()), float(data.sum()))


def composition(signal_run: Path, overlay_run: Path):
    """Process shares of the signal, both sides read from run directories (never typed in).

    GENIE: `composition_int_type_fiducial` of the `ndp signal` run (share of the fiducial truth signal).
    GiBUU: the cached sample named in the overlay run's provenance, restricted to the channel's signal,
    weighted, i.e. the share of the signal cross section.
    """
    from ndp.events import TruthTable, INT_TYPES
    from ndp.channels.registry import load_channel

    sc = json.loads((overlay_run / "scorecard.json").read_text())
    g = {k: v["fraction"] for k, v in
         json.loads((signal_run / "summary.json").read_text())["composition_int_type_fiducial"].items()}
    ch = load_channel(sc["channel"]["name"])
    t = TruthTable.load(Path(sc["provenance"]["source"]) / "truth.npz")
    w, m, it = t["weight"], ch.is_signal(t), np.asarray(t["int_type"])
    tot = float(w[m].sum())
    b = {INT_TYPES.get(int(c), str(c)): float(w[m & (it == c)].sum()) / tot for c in sorted(set(it[m].tolist()))}
    return g, b


def shape_stats(ratio_per_bin, norm):
    r = np.asarray(ratio_per_bin, float)
    r = r[np.isfinite(r)] / norm
    return float(np.std(r)), float(np.max(np.abs(r - 1)))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selection", required=True); ap.add_argument("--overlay", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--tag", default="cmp")
    ap.add_argument("--signal", help="the `ndp signal` run for the GENIE process composition; the section is skipped without it")
    a = ap.parse_args()
    sel = Path(a.selection).resolve(); ov = Path(a.overlay).resolve(); out = Path(a.out).resolve()
    figs = out.parent / "figs"; figs.mkdir(parents=True, exist_ok=True)
    gsum = json.loads((sel / "summary.json").read_text())
    osc = json.loads((ov / "scorecard.json").read_text())
    ots = json.loads((ov / "truth_summary.json").read_text())
    names = [n for n in ORDER if n in gsum["measurements"] and n in osc["measurements"]]
    rows, made = [], []
    for n in names:
        g = gsum["measurements"][n]
        b = osc["measurements"][n]["projections"]["full"]
        if "data" not in g:
            continue
        f = figs / f"{a.tag}_{n}.png"
        r_g, r_b, rb_g, rb_b, sig_g, sig_b, bkg, dat = make_figure(n, g, b, f)
        sg_rms, sg_max = shape_stats(rb_g, r_g); sb_rms, sb_max = shape_stats(rb_b, r_b)
        rows.append({"grid": n, "data": dat, "bkg": bkg, "sig_genie": sig_g, "sig_gibuu": sig_b,
                     "r_genie": r_g, "r_gibuu": r_b, "shape_genie": (sg_rms, sg_max), "shape_gibuu": (sb_rms, sb_max),
                     "bin1_genie": float(rb_g[0]), "bin1_gibuu": float(rb_b[0]),
                     "bin1_hi": float(np.asarray(g["edges"], float)[1])})
        made.append(f.name)

    L = [f"# GENIE and GiBUU side by side against the MINERvA muon + leading-proton data", "",
         f"**The sample.** This is not an exclusive one-muon-one-proton final state. The signal asks for a forward muon "
         f"(2–20 GeV/c, θ < 17°) and **at least one** proton in 0.5–1.1 GeV/c below 70°, vetoing mesons, baryons heavier "
         f"than the neutron and photons above 10 MeV; an event with a second proton in the window is kept. The muon and the "
         f"**highest-momentum** proton define the transverse-imbalance variables, a choice the paper makes explicitly because it "
         f"“is not changed based on the number of protons in the final state since secondary protons may not be reconstructed” "
         f"(arXiv:2503.15047). So every signal block below is one muon plus the leading proton, and the `1mu1p` in the channel "
         f"and run paths is only the identifier those manifests were created with, not a statement of proton multiplicity.", "",
         f"**What is plotted.** Every figure has two panels built identically: the same reconstructed data points, the same "
         f"official-MC background stacked by category, and on top of it the signal predicted by one generator — the official "
         f"GENIE MC scaled by POT on the left, GiBUU folded through the detector response learned from that same MC on the right. "
         f"Only the signal block differs between the panels, so the comparison is visual and direct. Both panels share the vertical scale.", "",
         f"**Sources.** GENIE: `{sel.name}`. GiBUU: `{ov.name}` "
         f"(σ_CC = {ots['meta'].get('sigma_flux_avg_per_nucleon_cm2', float('nan')):.4e} cm²/nucleon, "
         f"{ots['summary']['n_events']:,} events, {ots['meta'].get('flux_fraction_covered', float('nan')):.3f} of the flux). "
         f"Rendered by `report/make_genie_gibuu_comparison.py`; the physics discussion is in "
         f"`report/GiBUU_reco_comparison_FHC.md`.", "", "---", "",
         "## Summary", "",
         "| grid | data | background | GENIE signal | GiBUU signal | data/GENIE | data/GiBUU | GENIE shape rms | GiBUU shape rms |",
         "|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        L.append(f"| {LABELS.get(r['grid'], r['grid'])} | {r['data']:.0f} | {r['bkg']:.0f} | {r['sig_genie']:.0f} | {r['sig_gibuu']:.0f} | "
                 f"{r['r_genie']:.3f} | {r['r_gibuu']:.3f} | {r['shape_genie'][0]:.3f} | {r['shape_gibuu'][0]:.3f} |")
    mg = float(np.mean([r["r_genie"] for r in rows])); mb = float(np.mean([r["r_gibuu"] for r in rows]))
    L += ["", f"Averaged over the {len(rows)} released grids: data/GENIE = {mg:.3f}, data/GiBUU = {mb:.3f}. "
          "The shape columns are the root-mean-square of the per-bin ratio after each generator's own normalisation is divided out, "
          "so they measure shape only.", "",
          "## Figures", ""]
    for r, f in zip(rows, made):
        L += [f"### {LABELS.get(r['grid'], r['grid'])}", "",
              f"![{r['grid']}](figs/{f})", "",
              f"*data/GENIE {r['r_genie']:.3f}, shape rms {r['shape_genie'][0]:.3f} (worst bin {r['shape_genie'][1]:.3f}); "
              f"data/GiBUU {r['r_gibuu']:.3f}, shape rms {r['shape_gibuu'][0]:.3f} (worst bin {r['shape_gibuu'][1]:.3f}).*", ""]
    L += ["## Reading these plots", "",
          "- The background is **the same in both panels**: it comes from the official MC in both cases, because no generator "
          "predicts the non-signal part of this selection on its own. Only the coloured signal block on top is the model under test.",
          "- Both generators sit above the data. That offset is shared with the MC that supplies the background and the response, "
          "which is itself untuned, so it should not be read as a measured discrepancy of either generator.",
          "- The ratio panels are where the models separate: look at the first bins of p_n, δp_T and δp_L for the Fermi-motion peak, "
          "and at the low end of the muon angle and muon p_T for the missing low-Q² suppression.",
          "- The axes use the paper's released binning, so several grids end in one or two very wide bins (p_n to 6 GeV/c, "
          "δp_Ty to 3 GeV/c). Those bins hold few events and squeeze the interesting region; read them together with the numbers "
          "in the summary table.",
          "- The two panels of a figure share the vertical scale, so the height difference between the coloured signal blocks is "
          "the difference between the models at a glance.", ""]

    by = {r["grid"]: r for r in rows}
    if a.signal:
        gen_comp, buu_comp = composition(Path(a.signal).resolve(), ov)
        pretty = [("QE", "quasi-elastic"), ("RES", "resonance with the pion absorbed"), ("MEC", "2p2h"), ("DIS", "deep inelastic")]
        L += ["## Where the two models actually differ", "",
              "| | GENIE | GiBUU |", "|---|---|---|",
              f"| signal, averaged over the grids | {np.mean([r['sig_genie'] for r in rows]):,.0f} events | "
              f"{np.mean([r['sig_gibuu'] for r in rows]):,.0f} events |",
              f"| data / prediction | {mg:.3f} | {mb:.3f} |"]
        for key, label in pretty:
            if key in gen_comp or key in buu_comp:
                L.append(f"| {label} | {100 * gen_comp.get(key, float('nan')):.1f} % | {100 * buu_comp.get(key, float('nan')):.1f} % |")
        L += ["", f"The GENIE column is the share of the **fiducial truth signal** in `{Path(a.signal).name}`; the GiBUU column is "
              f"the share of the **signal cross section** in the cached sample behind `{ov.name}`. GiBUU builds a similar total out "
              f"of a visibly different mixture, and that is what the ratio panels show.", ""]

    al, dt, pn = by.get("alpha"), by.get("dpt"), by.get("pn")
    if al and dt and pn:
        L += [f"The clearest consequences: GiBUU describes δα_T better than GENIE (shape rms {al['shape_gibuu'][0]:.3f} against "
              f"{al['shape_genie'][0]:.3f}) and is the milder of the two in the low-Q² region, while GENIE describes δp_T and p_n "
              f"better ({dt['shape_genie'][0]:.3f} and {pn['shape_genie'][0]:.3f} against {dt['shape_gibuu'][0]:.3f} and "
              f"{pn['shape_gibuu'][0]:.3f}). In the first p_n bin, below {pn['bin1_hi']:.2g} GeV/c, the data stand "
              f"{100 * (pn['bin1_gibuu'] - 1):.0f} % above GiBUU and {100 * (1 - pn['bin1_genie']):.0f} % below GENIE, a genuine "
              f"reversal at the Fermi-motion peak.", ""]
    out.write_text("\n".join(L) + "\n")
    print(f"wrote {out} with {len(made)} side-by-side figures")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
