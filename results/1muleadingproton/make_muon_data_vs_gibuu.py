#!/usr/bin/env python
"""results/1muleadingproton, section 7: the full reconstructed space — data against GiBUU smeared with a VBLL surrogate
plus the official MC background, and the same data against the official MC, side by side.

    python results/1muleadingproton/make_muon_data_vs_gibuu.py [fhc6_het|fhc4_het|x60_het] [--replot]   (default pixi environment; ~6 min)

`--replot` skips the two passes over the playlists and redraws the figures from the histograms already recorded in
muon_data_vs_gibuu_<tag>.json.

Reconstructed muon momentum, angle and cos(theta) of every candidate passing the channel selection
`minerva_ccqelike_1mu1p_v0`, on the fine bins of the sections above:
  * data: the 12 ME FHC playlists (1.057e21 POT), selected candidates in the tuple's reconstructed kinematics;
  * official MC: the 12 ME FHC StandardMC playlists scaled by POT_data / POT_mc, selected candidates split into signal
    (truth signal inside the fiducial volume) and background by category (1 pi+-, 1 pi0, multi-pi, other), in the
    reconstructed kinematics — the reco_selected histogram of section 4 is the signal part;
  * GiBUU + VBLL: the smeared, cuts-applied signal of section 5 (muon_vbll_<tag>.json, or section 3's
    muon_vbll_smeared.json for the ported x60_het), i.e. GiBUU truth signal x selection efficiency x VBLL migration.
Prediction A (official MC) = MC signal + MC background; prediction B (GiBUU + VBLL) = GiBUU smeared signal + the same MC
background. Figures: an overlay of data over prediction B stacked by component with the data / prediction ratio, and
the two predictions against the same data side by side. Goodness of fit: Poisson -2lnL over the bins with a non-zero
prediction (ndp.compare.folded.poisson_gof), statistical only.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PLATFORM = HERE.parents[1]
sys.path.insert(0, str(PLATFORM))
import matplotlib                                                   # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                     # noqa: E402
from ndp.channels import load_channel                               # noqa: E402
from ndp.channels.selections import select                          # noqa: E402
from ndp.compare.folded import poisson_gof                          # noqa: E402
from ndp.config import load_site_config                             # noqa: E402
from ndp.diagnostics import mc_categories                           # noqa: E402
from ndp.events import TruthTable                                   # noqa: E402
from ndp.products import iter_reco_chunks, playlist_dir, total_pot  # noqa: E402
from ndp.surrogate.chunked import BKG_CATEGORIES                    # noqa: E402
from scripts import plot_style as ps                                # noqa: E402

CHANNEL = "minerva_me_ccqelike_1mu1p"
ARGS = [a for a in sys.argv[1:] if not a.startswith("--")]
MODEL = ARGS[0] if ARGS else "fhc6_het"
REPLOT = "--replot" in sys.argv
TAG = MODEL.replace("_het", "")
T0 = time.time()


def log(msg):
    print(f"[{TAG} {time.time() - T0:7.1f}s] {msg}", flush=True)


cfg = load_site_config(); ch = load_channel(CHANNEL); figs = HERE / "figs"
pot_data, pot_mc = total_pot(cfg, ch, "data"), total_pot(cfg, ch, "mc"); scale = pot_data / pot_mc
cos17 = float(np.cos(np.radians(17.0)))
specs = {"muon_p": (np.linspace(2, 20, 37), "reconstructed muon momentum [GeV/c]"),
         "muon_theta": (np.linspace(0, 17, 35), "reconstructed muon angle to the beam [deg]"),
         "muon_costheta": (np.linspace(cos17, 1.0, 45), "reconstructed muon cos(theta) to the beam")}
# the GiBUU + VBLL signal of the sections above, same fine bins
if MODEL == "x60_het":
    src = HERE / "muon_vbll_smeared.json"; j = json.load(open(src)); gib = {n: np.array(j[n]["vbll_smeared_cuts_applied"], float) for n in specs}
    gib_edges = {n: np.array(j[n]["edges"], float) for n in specs}; section = "section 3"
else:
    src = HERE / f"muon_vbll_{TAG}.json"; j = json.load(open(src))["gibuu"]["variables"]
    gib = {n: np.array(j[n][f"{TAG}_smeared_cuts_applied"], float) for n in specs}; gib_edges = {n: np.array(j[n]["edges"], float) for n in specs}
    section = "section 5" if MODEL == "fhc6_het" else "section 6"
for n, (e, _) in specs.items():
    assert np.allclose(e, gib_edges[n]), f"{n}: fine bins differ from {src.name}"


def reco_vars(r, keep):
    p = np.asarray(r["reco_p"], float)[keep]; th = np.degrees(np.asarray(r["reco_theta"], float)[keep])
    return {"muon_p": p, "muon_theta": th, "muon_costheta": np.cos(np.radians(th))}


# ---- data: selected candidates, reconstructed kinematics -----------------------------------------------------------
CATS = ("signal",) + tuple(BKG_CATEGORIES)
D = {n: np.zeros(len(e) - 1) for n, (e, _) in specs.items()}; n_data = 0; chunks = []
M = {n: {c: np.zeros(len(e) - 1) for c in CATS} for n, (e, _) in specs.items()}; n_mc = {c: 0 for c in CATS}
if REPLOT:   # raw MC counts = the recorded POT-scaled histograms / scale
    prev = json.load(open(HERE / f"muon_data_vs_gibuu_{TAG}.json")); n_data = prev["n_data_selected"]; chunks = prev["inputs"]; n_mc = prev["n_mc_selected_raw"]
    for n in specs:
        v = prev["variables"][n]; D[n] = np.array(v["data"], float); M[n]["signal"] = np.array(v["mc_signal"], float) / scale
        for c in BKG_CATEGORIES:
            M[n][c] = np.array(v["mc_background"][c], float) / scale
    log(f"replot from muon_data_vs_gibuu_{TAG}.json")
for label, r, pot, _ in ([] if REPLOT else iter_reco_chunks(cfg, ch, "data")):
    keep = select(ch, r); v = reco_vars(r, keep)
    for n, (e, _) in specs.items():
        D[n] += np.histogram(v[n], bins=e)[0]
    n_data += int(keep.sum()); chunks.append({"kind": "data", "input": label, "pot": pot, "n_selected": int(keep.sum())})
    log(f"data {label}: {keep.sum()} selected of {len(keep)} candidates"); del r
# ---- official MC: selected candidates split into signal and background categories, reconstructed kinematics ---------
for label, r, pot, _ in ([] if REPLOT else iter_reco_chunks(cfg, ch, "mc")):
    beam, pl = label.split("/", 1)
    rt = TruthTable.load(playlist_dir(cfg, ch, beam, pl) / f"reco_{pl}_truthcols_skim.npz")
    keep = select(ch, r); cats = np.asarray(mc_categories(ch, rt))[keep]
    sig = np.asarray(rt["signal_minerva_ccqelike_1mu1p"], bool)[keep] & ch.in_phase_space(rt)[keep]
    n_dis = int((sig != np.char.startswith(cats.astype(str), "signal")).sum())
    if n_dis:
        log(f"WARNING {label}: the skim's signal flag and mc_categories disagree on {n_dis} selected candidates; categories win")
    sig = np.char.startswith(cats.astype(str), "signal")
    v = reco_vars(r, keep)
    for c in CATS:
        m = sig if c == "signal" else (cats == c)
        for n, (e, _) in specs.items():
            M[n][c] += np.histogram(v[n][m], bins=e)[0]
        n_mc[c] += int(m.sum())
    chunks.append({"kind": "mc", "input": label, "pot": pot, "n_selected": int(keep.sum()), "n_signal": int(sig.sum()), "n_background": int((~sig).sum())})
    log(f"mc {label}: {keep.sum()} selected, {sig.sum()} signal, {(~sig).sum()} background"); del r, rt

# ---- figures + numbers -------------------------------------------------------------------------------------------------
oi = ps.okabe_ito()
c_data, c_mc_sig, c_gib = ps.color_for("data"), ps.color_for("mc_signal"), oi["vermilion"]
c_bkg = {"bkg 1 pi+-": oi["sky"], "bkg 1 pi0": oi["green"], "bkg multi-pi": oi["purple"], "bkg other (no pion)": "#9e9e9e"}
short = {"bkg 1 pi+-": "background: 1 pi+-", "bkg 1 pi0": "background: 1 pi0", "bkg multi-pi": "background: multi-pi", "bkg other (no pion)": "background: other (no pion)"}
out = {"model": MODEL, "gibuu_signal_source": src.name, "pot_data": pot_data, "pot_mc": pot_mc, "scale_to_data_pot": scale, "selection": ch.selection["name"],
       "n_data_selected": n_data, "n_mc_selected_raw": {c: n_mc[c] for c in CATS}, "n_mc_selected_at_data_pot": {c: n_mc[c] * scale for c in CATS},
       "mc_purity": n_mc["signal"] / sum(n_mc.values()), "inputs": chunks, "variables": {}}
log(f"data selected {n_data}; MC selected {sum(n_mc.values())} raw ({n_mc['signal']} signal, {sum(n_mc.values()) - n_mc['signal']} background), "
    f"purity {out['mc_purity']:.3f}; POT scale {scale:.4f}")


def stack(ax, e, bkg_parts, sig, sig_color, sig_label):
    """Background categories stacked from the bottom, the signal on top; returns the total."""
    w = np.diff(e); bottom = np.zeros(len(e) - 1)
    for c in BKG_CATEGORIES:
        ax.bar(e[:-1], bkg_parts[c], width=w, bottom=bottom, align="edge", color=c_bkg[c], lw=0, label=short[c]); bottom += bkg_parts[c]
    ax.bar(e[:-1], sig, width=w, bottom=bottom, align="edge", color=sig_color, alpha=0.75, lw=0, label=sig_label)
    tot = bottom + sig; ax.stairs(tot, e, color="black", lw=0.9, baseline=None)
    return tot


def data_points(ax, e, d, label):
    c = 0.5 * (e[:-1] + e[1:])
    ax.errorbar(c, d, xerr=np.diff(e) / 2, yerr=np.sqrt(d), fmt=ps.marker_for("data"), ms=ps.style("data_markersize"), color=c_data, capsize=2, label=label, zorder=5)


def ratio_panel(ax, e, d, pred, var_pred, band_label):
    c = 0.5 * (e[:-1] + e[1:])
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where(pred > 0, d / pred, np.nan); rerr = np.where(pred > 0, np.sqrt(d) / pred, np.nan); band = np.where(pred > 0, np.sqrt(var_pred) / pred, 0.0)
    ax.fill_between(e, 1 - np.append(band, band[-1]), 1 + np.append(band, band[-1]), step="post", color=ps.color_for("reference"), alpha=0.3, lw=0, label=band_label)
    ax.errorbar(c, r, xerr=np.diff(e) / 2, yerr=rerr, fmt=ps.marker_for("ratio"), ms=ps.style("ratio_markersize"), color=ps.color_for("ratio"), capsize=2)
    ax.axhline(1.0, color=ps.color_for("reference"), ls=ps.style("ratio_linestyle"), lw=ps.style("ratio_linewidth"))
    ax.set_ylim(0.4, 1.6); ax.grid(alpha=0.3); ax.legend(fontsize=7, loc="upper left")
    return r


for n, (e, xlabel) in specs.items():
    d = D[n]; sig_mc = M[n]["signal"] * scale; bkg = {c: M[n][c] * scale for c in BKG_CATEGORIES}; bkg_tot = sum(bkg.values())
    var_bkg = sum(M[n][c] for c in BKG_CATEGORIES) * scale ** 2           # MC-stat variance of the POT-scaled background
    var_mc_full = var_bkg + M[n]["signal"] * scale ** 2
    sig_gib = gib[n]; pred_mc = sig_mc + bkg_tot; pred_gib = sig_gib + bkg_tot
    gof_mc, gof_gib = poisson_gof(d, pred_mc, var_mc_full), poisson_gof(d, pred_gib, var_bkg)
    lab_gib = f"GiBUU signal x efficiency, VBLL {MODEL} smeared ({section})"
    lab_mc = "official MC signal (GENIE, MINERvA reconstruction)"
    lab_data = f"data, FHC 1A-1P ({n_data} selected)"
    ymax = 1.5 * max(float(np.max(pred_mc)), float(np.max(pred_gib)), float(np.max(d)))

    # (a) overlay: data over GiBUU + VBLL signal on the MC background
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(7.5, 7.2), sharex=True, gridspec_kw={"height_ratios": [3, 1.2], "hspace": 0.06})
    stack(ax, e, bkg, sig_gib, c_gib, lab_gib); data_points(ax, e, d, lab_data)
    ax.set_ylabel("selected events at 1.057e21 POT / bin"); ax.set_ylim(0, ymax); ax.grid(alpha=0.3); ax.legend(fontsize=7.5, loc="upper right")
    ax.set_title(f"Full reconstructed space, 1mu1p selection: data vs GiBUU + VBLL {MODEL} with the official MC background\n"
                 f"data {d.sum():.0f}, prediction {pred_gib.sum():.0f} (signal {sig_gib.sum():.0f} + background {bkg_tot.sum():.0f}), "
                 f"data/pred {d.sum() / pred_gib.sum():.3f}, -2lnL/ndf {gof_gib['minus2lnL']:.0f}/{gof_gib['ndf']}", fontsize=9)
    r_gib = ratio_panel(axr, e, d, pred_gib, var_bkg, "MC stat of the background")
    axr.set_ylabel("data / prediction"); axr.set_xlabel(xlabel)
    fig.savefig(figs / f"{n}_data_vs_gibuu_{TAG}.png", dpi=150, bbox_inches="tight"); plt.close(fig)

    # (b) side by side: official MC | GiBUU + VBLL, same data, same y range
    fig, axes = plt.subplots(2, 2, figsize=(14.5, 7.2), sharex=True, gridspec_kw={"height_ratios": [3, 1.2], "hspace": 0.06, "wspace": 0.16})
    stack(axes[0, 0], e, bkg, sig_mc, c_mc_sig, lab_mc); data_points(axes[0, 0], e, d, lab_data)
    axes[0, 0].set_title(f"official MC: GENIE signal + background\nprediction {pred_mc.sum():.0f}, data/pred {d.sum() / pred_mc.sum():.3f}, "
                         f"-2lnL/ndf {gof_mc['minus2lnL']:.0f}/{gof_mc['ndf']}", fontsize=9)
    stack(axes[0, 1], e, bkg, sig_gib, c_gib, lab_gib); data_points(axes[0, 1], e, d, lab_data)
    axes[0, 1].set_title(f"GiBUU + VBLL {MODEL} signal + official MC background\nprediction {pred_gib.sum():.0f}, data/pred {d.sum() / pred_gib.sum():.3f}, "
                         f"-2lnL/ndf {gof_gib['minus2lnL']:.0f}/{gof_gib['ndf']}", fontsize=9)
    for a in axes[0]:
        a.set_ylim(0, ymax); a.grid(alpha=0.3); a.legend(fontsize=7, loc="upper right")
    axes[0, 0].set_ylabel("selected events at 1.057e21 POT / bin")
    r_mc = ratio_panel(axes[1, 0], e, d, pred_mc, var_mc_full, "MC stat (signal + background)")
    ratio_panel(axes[1, 1], e, d, pred_gib, var_bkg, "MC stat of the background")
    axes[1, 0].set_ylabel("data / prediction")
    for a in axes[1]:
        a.set_xlabel(xlabel)
    fig.suptitle(f"Full reconstructed space, 1mu1p selection, data FHC 1A-1P ({d.sum():.0f} events): {xlabel}", fontsize=11, y=0.98)
    fig.savefig(figs / f"{n}_data_vs_mc_and_gibuu_{TAG}.png", dpi=150, bbox_inches="tight"); plt.close(fig)

    c = 0.5 * (e[:-1] + e[1:])
    def median(h):
        cum = np.cumsum(h) / h.sum(); return float(c[np.searchsorted(cum, 0.5)])
    out["variables"][n] = {"edges": e.tolist(), "data": d.tolist(), "mc_signal": sig_mc.tolist(), "mc_background": {k: v.tolist() for k, v in bkg.items()},
                           "mc_background_total": bkg_tot.tolist(), "gibuu_vbll_signal": sig_gib.tolist(), "pred_official_mc": pred_mc.tolist(), "pred_gibuu_vbll": pred_gib.tolist(),
                           "ratio_data_over_official_mc": r_mc.tolist(), "ratio_data_over_gibuu_vbll": r_gib.tolist(),
                           "totals": {"data": float(d.sum()), "mc_signal": float(sig_mc.sum()), "mc_background": float(bkg_tot.sum()),
                                      "mc_background_by_category": {k: float(v.sum()) for k, v in bkg.items()}, "gibuu_vbll_signal": float(sig_gib.sum()),
                                      "pred_official_mc": float(pred_mc.sum()), "pred_gibuu_vbll": float(pred_gib.sum()),
                                      "data_over_official_mc": float(d.sum() / pred_mc.sum()), "data_over_gibuu_vbll": float(d.sum() / pred_gib.sum()),
                                      "background_fraction_of_official_mc": float(bkg_tot.sum() / pred_mc.sum())},
                           "gof_official_mc": gof_mc, "gof_gibuu_vbll": gof_gib,
                           "median": {"data": median(d), "pred_official_mc": median(pred_mc), "pred_gibuu_vbll": median(pred_gib)},
                           "peak_bin": {"data": [float(e[d.argmax()]), float(e[d.argmax() + 1])], "pred_official_mc": [float(e[pred_mc.argmax()]), float(e[pred_mc.argmax() + 1])],
                                        "pred_gibuu_vbll": [float(e[pred_gib.argmax()]), float(e[pred_gib.argmax() + 1])]},
                           "ratio_range": {"official_mc": [float(np.nanmin(r_mc)), float(np.nanmax(r_mc))], "gibuu_vbll": [float(np.nanmin(r_gib)), float(np.nanmax(r_gib))]}}
    t = out["variables"][n]["totals"]
    log(f"{n}: data {t['data']:.0f} | MC {t['pred_official_mc']:.0f} (sig {t['mc_signal']:.0f} + bkg {t['mc_background']:.0f}) data/MC {t['data_over_official_mc']:.3f} "
        f"-2lnL/ndf {gof_mc['minus2lnL']:.0f}/{gof_mc['ndf']} | GiBUU+VBLL {t['pred_gibuu_vbll']:.0f} (sig {t['gibuu_vbll_signal']:.0f}) data/pred {t['data_over_gibuu_vbll']:.3f} "
        f"-2lnL/ndf {gof_gib['minus2lnL']:.0f}/{gof_gib['ndf']} | ratio ranges MC {np.nanmin(r_mc):.2f}-{np.nanmax(r_mc):.2f}, GiBUU {np.nanmin(r_gib):.2f}-{np.nanmax(r_gib):.2f}")
json.dump(out, open(HERE / f"muon_data_vs_gibuu_{TAG}.json", "w"), indent=1)
log(f"wrote muon_data_vs_gibuu_{TAG}.json + figures")
