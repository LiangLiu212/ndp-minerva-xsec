#!/usr/bin/env python
"""results/1muleadingproton, section 3: migration matrices of the official MC on the same bins as the VBLL ones.

    python results/1muleadingproton/make_muon_mc_migration.py        (default pixi environment; ~10 min, one playlist at a time)

Source: the MINERvA Open Data StandardMC of the 12 ME FHC playlists (products/FHC/<pl>/reco_<pl>_mc.npz + the reco rows'
truth columns), the same rows the platform's binned responses were learned from. Events: reco candidates passing the
channel selection `minerva_ccqelike_1mu1p_v0` whose truth is signal inside the fiducial volume (the responses' numerator).
True muon momentum and beam angle from the primary lepton (rotated from the detector frame), reconstructed ones from
the cached `reco_p` and `reco_theta` (MasterAnaDev_leptonE and the beam-frame angle branches). Unweighted CV MC.
Matrices are column-normalised to P(reco bin | true bin); on the analysis grids they must reproduce the migration
counts stored in the binned responses (checked). Also draws the diagonal fractions of MC and VBLL side by side.
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
from ndp.channels import load_channel, load_measurement, observables as obs   # noqa: E402
from ndp.channels.selections import select                          # noqa: E402
from ndp.config import load_site_config                             # noqa: E402
from ndp.events import TruthTable                                   # noqa: E402
from ndp.pipeline import resolve_surrogate                          # noqa: E402
from ndp.products import iter_reco_chunks, playlist_dir             # noqa: E402
from scripts import plot_style as ps                                # noqa: E402

CHANNEL = "minerva_me_ccqelike_1mu1p"
T0 = time.time()
cfg = load_site_config(); ch = load_channel(CHANNEL)
cos17 = float(np.cos(np.radians(17.0)))
specs = {"muon_p": ("muon_p", np.linspace(2, 20, 37), "muon momentum [GeV/c]"),
         "muon_theta": ("muon_theta", np.linspace(0, 17, 35), "muon angle to the beam [deg]"),
         "muon_costheta": ("muon_costheta", np.linspace(cos17, 1.0, 45), "muon cos(theta) to the beam")}
grids = {n: np.asarray(load_measurement(ch, g).x.edges, float) for n, (g, _, _) in specs.items()}
H = {n: {"analysis_grid": np.zeros((len(grids[n]) - 1,) * 2), "fine": np.zeros((len(e) - 1,) * 2)} for n, (_, e, _) in specs.items()}
lost = {n: {"analysis_grid": np.zeros(len(grids[n]) - 1), "fine": np.zeros(len(specs[n][1]) - 1)} for n in specs}
n_sel_sig, inputs = 0, []
for label, r, pot, src in iter_reco_chunks(cfg, ch, "mc"):
    beam, pl = label.split("/", 1)
    rt = TruthTable.load(playlist_dir(cfg, ch, beam, pl) / f"reco_{pl}_truthcols_skim.npz")
    keep = select(ch, r) & np.asarray(rt["signal_minerva_ccqelike_1mu1p"], bool) & ch.in_phase_space(rt)
    tp, tth = obs.lep_p(rt, "beam")[keep], obs.lep_theta_deg(rt, "beam")[keep]
    rp, rth = np.asarray(r["reco_p"], float)[keep], np.degrees(np.asarray(r["reco_theta"], float)[keep])
    vals = {"muon_p": (tp, rp), "muon_theta": (tth, rth), "muon_costheta": (np.cos(np.radians(tth)), np.cos(np.radians(rth)))}
    for n, (xt, xr) in vals.items():
        for tag, e in (("analysis_grid", grids[n]), ("fine", specs[n][1])):
            H[n][tag] += np.histogram2d(xr, xt, bins=[e, e])[0]
            inside = (xr >= e[0]) & (xr < e[-1])
            lost[n][tag] += np.histogram(xt[~inside], bins=e)[0]
    n_sel_sig += int(keep.sum()); inputs.append({"input": label, "n_selected_signal": int(keep.sum())})
    print(f"[{time.time()-T0:6.1f}s] {label}: {keep.sum()} selected signal candidates", flush=True)
    del r, rt

vb = json.load(open(HERE / "muon_vbll_smeared.json"))["migration"]
out = {"n_selected_signal_in_fiducial": n_sel_sig, "inputs": inputs, "selection": ch.selection["name"], "matrices": {}}
figs = HERE / "figs"
for n, (g, fine, xlabel) in specs.items():
    coarse = grids[n]
    res = {}
    for tag in ("analysis_grid", "fine"):
        col = H[n][tag].sum(axis=0) + lost[n][tag]
        with np.errstate(invalid="ignore", divide="ignore"):
            res[tag] = np.where(col > 0, H[n][tag] / col, np.nan)
    Pc = res["analysis_grid"]; diag = np.diag(Pc)
    sur, sp = resolve_surrogate(ch, load_measurement(ch, g), cfg)
    same = bool(np.allclose(H[n]["analysis_grid"], sur.migration_counts))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5.6))
    im1 = a1.pcolormesh(coarse, coarse, Pc * 100.0, cmap="Blues", vmin=0, vmax=100, shading="flat")
    for j in range(len(coarse) - 1):
        for i in range(len(coarse) - 1):
            if np.isfinite(Pc[i, j]) and Pc[i, j] >= 0.005:
                a1.text(0.5 * (coarse[j] + coarse[j + 1]), 0.5 * (coarse[i] + coarse[i + 1]), f"{100 * Pc[i, j]:.0f}", ha="center", va="center",
                        fontsize=7 if len(coarse) > 10 else 8, color="white" if Pc[i, j] > 0.5 else "black")
    a1.set_xlabel(f"true {xlabel}"); a1.set_ylabel(f"reconstructed (MINERvA) {xlabel}"); a1.set_title(f"P(reco bin | true bin) [%] on the {g} grid", fontsize=10)
    fig.colorbar(im1, ax=a1, label="%")
    im2 = a2.pcolormesh(fine, fine, res["fine"] * 100.0, cmap="Blues", vmin=0, shading="flat")
    a2.plot([fine[0], fine[-1]], [fine[0], fine[-1]], color=ps.color_for("reference"), lw=0.8, ls="--")
    a2.set_xlabel(f"true {xlabel}"); a2.set_ylabel(f"reconstructed (MINERvA) {xlabel}"); a2.set_title("same, fine bins", fontsize=10)
    fig.colorbar(im2, ax=a2, label="%")
    fig.suptitle(f"Official MC migration of the selected signal (FHC 1A-1P StandardMC): {xlabel}", fontsize=11)
    fig.tight_layout(); fig.savefig(figs / f"{n}_mc_migration.png", dpi=150, bbox_inches="tight"); plt.close(fig)
    out["matrices"][n] = {"grid": g, "grid_edges": coarse.tolist(), "P_reco_given_true_percent_reco_by_true": (np.nan_to_num(Pc) * 100).round(2).tolist(),
                          "counts_reco_by_true": H[n]["analysis_grid"].astype(int).tolist(), "diagonal_percent": (np.nan_to_num(diag) * 100).round(1).tolist(),
                          "loss_fraction_per_true_bin": np.nan_to_num(lost[n]["analysis_grid"] / np.maximum(H[n]["analysis_grid"].sum(0) + lost[n]["analysis_grid"], 1)).round(4).tolist(),
                          "weighted_mean_diagonal_percent": float(100 * np.nansum(diag * H[n]["analysis_grid"].sum(0)) / H[n]["analysis_grid"].sum()),
                          "matches_binned_response_migration_counts": same, "binned_response": str(sp.relative_to(PLATFORM)),
                          "vbll_diagonal_percent": vb[n]["diagonal_percent"], "vbll_weighted_mean_diagonal_percent": vb[n]["weighted_mean_diagonal_percent"]}
    print(f"{n}: MC diagonal % {out['matrices'][n]['diagonal_percent']} (mean {out['matrices'][n]['weighted_mean_diagonal_percent']:.1f}) | VBLL {vb[n]['diagonal_percent']} "
          f"(mean {vb[n]['weighted_mean_diagonal_percent']:.1f}) | matches binned response: {same}", flush=True)

# diagonal fractions side by side
fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
for ax, (n, (g, fine, xlabel)) in zip(axes, specs.items()):
    e = grids[n]; c = 0.5 * (e[:-1] + e[1:]); wdt = np.diff(e)
    ax.stairs(np.array(out["matrices"][n]["diagonal_percent"]), e, color=ps.color_for("data"), lw=1.8, baseline=None, label="official MC")
    ax.stairs(np.array(vb[n]["diagonal_percent"]), e, color=ps.okabe_ito()["vermilion"], lw=1.8, baseline=None, label="VBLL x60_het")
    ax.set_xlabel(f"true {xlabel}"); ax.set_ylabel("P(same bin) [%]"); ax.set_ylim(0, 100); ax.grid(alpha=0.3); ax.set_title(f"{g} grid", fontsize=10)
axes[0].legend(fontsize=9)
fig.suptitle("Probability of reconstructing in the true bin: official MC vs the VBLL surrogate (selected 1mu1p signal)", fontsize=11)
fig.tight_layout(); fig.savefig(figs / "migration_diagonal_mc_vs_vbll.png", dpi=150, bbox_inches="tight"); plt.close(fig)
json.dump(out, open(HERE / "muon_mc_migration.json", "w"), indent=1)
print(f"[{time.time()-T0:6.1f}s] total selected signal {n_sel_sig}; wrote muon_mc_migration.json + figures", flush=True)
