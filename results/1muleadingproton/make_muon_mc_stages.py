#!/usr/bin/env python
"""results/1muleadingproton, section 4: the official MC signal channel at three stages — true without cut, true with the
selection cuts, and reconstructed.

    python results/1muleadingproton/make_muon_mc_stages.py        (default pixi environment; ~6 min, one playlist at a time)

Source: MINERvA Open Data StandardMC, the 12 ME FHC playlists (4.978e21 POT), scaled to the data exposure (1.057e21 POT).
  1. true, no cut: every truth signal event (channel signal definition) with its true vertex in the tracker fiducial volume,
     from the playlist truth skims — the efficiency denominator;
  2. true, cuts applied: the reconstructed candidates passing the channel selection `minerva_ccqelike_1mu1p_v0` whose truth
     is signal in the fiducial volume, plotted in their TRUE muon kinematics — the efficiency numerator;
  3. reconstructed: the same candidates in their RECONSTRUCTED muon kinematics (`reco_p`, `reco_theta`).
Angles are with respect to the beam (true: primary lepton rotated from the detector frame; reco: the tuple's beam-frame
angle). Lower panels: the selection efficiency per bin (2 / 1) and the migration ratio (3 / 2). Unweighted CV MC.
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
from ndp.channels import load_channel, observables as obs           # noqa: E402
from ndp.channels.selections import select                          # noqa: E402
from ndp.config import load_site_config                             # noqa: E402
from ndp.events import TruthTable                                   # noqa: E402
from ndp.products import iter_reco_chunks, playlist_dir, total_pot  # noqa: E402
from scripts import plot_style as ps                                # noqa: E402

CHANNEL = "minerva_me_ccqelike_1mu1p"
T0 = time.time()
cfg = load_site_config(); ch = load_channel(CHANNEL)
pot_data, pot_mc = total_pot(cfg, ch, "data"), total_pot(cfg, ch, "mc")
scale = pot_data / pot_mc
cos17 = float(np.cos(np.radians(17.0)))
specs = {"muon_p": (np.linspace(2, 20, 37), "muon momentum [GeV/c]"),
         "muon_theta": (np.linspace(0, 17, 35), "muon angle to the beam [deg]"),
         "muon_costheta": (np.linspace(cos17, 1.0, 45), "muon cos(theta) to the beam")}
H = {n: {"true_all": np.zeros(len(e) - 1), "true_selected": np.zeros(len(e) - 1), "reco_selected": np.zeros(len(e) - 1)} for n, (e, _) in specs.items()}
n_all, n_sel, inputs = 0, 0, []
sums = {n: {"true_all": [0.0, 0.0], "true_selected": [0.0, 0.0], "reco_selected": [0.0, 0.0]} for n in specs}   # sum x, count for means


def add(n, key, x):
    H[n][key] += np.histogram(x, bins=specs[n][0])[0]
    sums[n][key][0] += float(x.sum()); sums[n][key][1] += len(x)


for label, r, pot, src in iter_reco_chunks(cfg, ch, "mc"):
    beam, pl = label.split("/", 1)
    d = playlist_dir(cfg, ch, beam, pl)
    t = TruthTable.load(d / f"truth_{pl}_skim.npz")
    fid = ch.in_phase_space(t) & ch.is_signal(t)
    tp, tth = obs.lep_p(t, "beam")[fid], obs.lep_theta_deg(t, "beam")[fid]
    add("muon_p", "true_all", tp); add("muon_theta", "true_all", tth); add("muon_costheta", "true_all", np.cos(np.radians(tth)))
    n_all += int(fid.sum()); del t
    rt = TruthTable.load(d / f"reco_{pl}_truthcols_skim.npz")
    keep = select(ch, r) & np.asarray(rt["signal_minerva_ccqelike_1mu1p"], bool) & ch.in_phase_space(rt)
    sp, sth = obs.lep_p(rt, "beam")[keep], obs.lep_theta_deg(rt, "beam")[keep]
    rp, rth = np.asarray(r["reco_p"], float)[keep], np.degrees(np.asarray(r["reco_theta"], float)[keep])
    add("muon_p", "true_selected", sp); add("muon_theta", "true_selected", sth); add("muon_costheta", "true_selected", np.cos(np.radians(sth)))
    add("muon_p", "reco_selected", rp); add("muon_theta", "reco_selected", rth); add("muon_costheta", "reco_selected", np.cos(np.radians(rth)))
    n_sel += int(keep.sum()); inputs.append({"input": label, "n_truth_signal_fiducial": int(fid.sum()), "n_selected_signal": int(keep.sum())})
    print(f"[{time.time()-T0:6.1f}s] {label}: {fid.sum()} fiducial truth signal, {keep.sum()} selected", flush=True)
    del r, rt

figs = HERE / "figs"
c1, c2, c3, ce = ps.color_for("mc_signal"), ps.okabe_ito()["blue"], ps.okabe_ito()["vermilion"], ps.color_for("efficiency")
out = {"pot_data": pot_data, "pot_mc": pot_mc, "scale_to_data_pot": scale, "n_truth_signal_fiducial": n_all, "n_selected_signal": n_sel,
       "mean_efficiency": n_sel / n_all, "inputs": inputs, "selection": ch.selection["name"], "variables": {}}
for n, (e, xlabel) in specs.items():
    ha, hs, hr = (H[n][k] * scale for k in ("true_all", "true_selected", "reco_selected"))
    with np.errstate(invalid="ignore", divide="ignore"):
        eff = np.where(ha > 0, hs / ha, np.nan); mig = np.where(hs > 0, hr / hs, np.nan)
    fig, (ax, ae, am) = plt.subplots(3, 1, figsize=(7, 8.2), sharex=True, gridspec_kw={"height_ratios": [3, 1.1, 1.1], "hspace": 0.06})
    ax.stairs(ha, e, color=c1, lw=1.6, baseline=None, label="true, no cut (fiducial truth signal)")
    ax.stairs(hs, e, color=c2, lw=1.6, baseline=None, label="true, selection cuts applied")
    ax.stairs(hr, e, color=c3, lw=1.8, baseline=0, fill=True, alpha=0.2)
    ax.stairs(hr, e, color=c3, lw=1.8, baseline=None, label="reconstructed, selection cuts applied")
    ax.set_ylabel("signal events at 1.057e21 POT (fiducial)"); ax.set_ylim(bottom=0); ax.grid(alpha=0.3); ax.legend(fontsize=8)
    ax.set_title(f"Official MC (FHC 1A-1P StandardMC), 1mu1p signal: {xlabel}", fontsize=10)
    ae.stairs(eff, e, color=ce, lw=1.6, baseline=None); ae.set_ylabel("efficiency\n(cuts / no cut)"); ae.set_ylim(0, max(0.6, np.nanmax(eff) * 1.15)); ae.grid(alpha=0.3)
    am.stairs(mig, e, color=c3, lw=1.6, baseline=None); am.axhline(1.0, color=ps.color_for("reference"), ls="--", lw=0.8)
    am.set_ylabel("reco / true\n(selected)"); am.set_ylim(0, max(2.0, np.nanmax(mig) * 1.1)); am.grid(alpha=0.3); am.set_xlabel(xlabel)
    fig.savefig(figs / f"{n}_mc_stages.png", dpi=150, bbox_inches="tight"); plt.close(fig)
    c = 0.5 * (e[:-1] + e[1:])
    v = {"edges": e.tolist(), "true_all": ha.tolist(), "true_selected": hs.tolist(), "reco_selected": hr.tolist(),
         "efficiency_per_bin": eff.tolist(), "reco_over_true_selected": mig.tolist(),
         "totals_at_data_pot": {"true_all": float(ha.sum()), "true_selected": float(hs.sum()), "reco_selected": float(hr.sum())}}
    for k in ("true_all", "true_selected", "reco_selected"):
        h = {"true_all": ha, "true_selected": hs, "reco_selected": hr}[k]; cum = np.cumsum(h) / h.sum()
        v[f"median_{k}"] = float(c[np.searchsorted(cum, 0.5)]); v[f"mean_{k}"] = sums[n][k][0] / sums[n][k][1]
        v[f"peak_bin_{k}"] = [float(e[h.argmax()]), float(e[h.argmax() + 1])]
    out["variables"][n] = v
    print(f"{n}: totals {v['totals_at_data_pot']} | medians all {v['median_true_all']:.4g} sel-true {v['median_true_selected']:.4g} reco {v['median_reco_selected']:.4g} "
          f"| reco/true range {np.nanmin(mig):.2f}-{np.nanmax(mig):.2f}", flush=True)
json.dump(out, open(HERE / "muon_mc_stages.json", "w"), indent=1)
print(f"[{time.time()-T0:6.1f}s] fiducial truth signal {n_all}, selected {n_sel}, efficiency {n_sel/n_all:.4f}; POT scale {scale:.4f}", flush=True)
