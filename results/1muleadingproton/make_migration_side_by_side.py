#!/usr/bin/env python
"""results/1muleadingproton, section 5: official-MC and fhc6_het migration matrices side by side.

    python results/1muleadingproton/make_migration_side_by_side.py        (default pixi environment)

Reads the analysis-grid matrices already recorded in muon_mc_migration.json (official MC, selected signal of the 12
FHC playlists) and muon_vbll_fhc6.json (fhc6_het applied to the selected GiBUU signal), both column-normalised to
P(reco bin | true bin) in percent, and draws them next to each other with their difference (fhc6_het minus MC, in
percentage points). The true-bin populations differ (MINERvA StandardMC vs GiBUU), so column normalisation is what
makes the two comparable.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PLATFORM = HERE.parents[1]
sys.path.insert(0, str(PLATFORM))
import matplotlib                                                   # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                     # noqa: E402

mc = json.load(open(HERE / "muon_mc_migration.json"))["matrices"]
MODEL = sys.argv[1] if len(sys.argv) > 1 else "fhc6_het"
TAG = MODEL.replace("_het", "")
# the surrogate's GiBUU migration: muon_vbll_<tag>.json (platform-trained models, key gibuu/migration) or, for the ported
# x60_het, the section-3 file muon_vbll_smeared.json (key migration)
src = HERE / ("muon_vbll_smeared.json" if MODEL == "x60_het" else f"muon_vbll_{TAG}.json")
j = json.load(open(src)); nw = j["migration"] if MODEL == "x60_het" else j["gibuu"]["migration"]
labels = {"muon_p": "muon momentum [GeV/c]", "muon_theta": "muon angle to the beam [deg]", "muon_costheta": "muon cos(theta) to the beam"}
out = {}


def annotate(ax, e, M, small, signed=False):
    for j in range(len(e) - 1):
        for i in range(len(e) - 1):
            v = M[i, j]
            if np.isfinite(v) and abs(v) >= 0.5:
                ax.text(0.5 * (e[j] + e[j + 1]), 0.5 * (e[i] + e[i + 1]), f"{v:+.0f}" if signed else f"{v:.0f}", ha="center", va="center",
                        fontsize=6.5 if small else 8, color="white" if (not signed and v > 50) else "black")


for name, xlabel in labels.items():
    e = np.array(mc[name]["grid_edges"], float); assert np.allclose(e, nw[name]["grid_edges"])
    A = np.array(mc[name]["P_reco_given_true_percent_reco_by_true"], float); B = np.array(nw[name]["P_reco_given_true_percent_reco_by_true"], float)
    D = B - A
    small = len(e) > 10
    fig, axes = plt.subplots(1, 3, figsize=(18.5, 5.6))
    for ax, M, title, cmap, vmin, vmax, signed in ((axes[0], A, "official MC (selected signal, FHC 1A-1P)", "Blues", 0, 100, False),
                                                   (axes[1], B, f"VBLL {MODEL} (selected GiBUU signal)", "Blues", 0, 100, False),
                                                   (axes[2], D, f"{MODEL} minus MC [percentage points]", "RdBu_r", -30, 30, True)):
        im = ax.pcolormesh(e, e, M, cmap=cmap, vmin=vmin, vmax=vmax, shading="flat")
        annotate(ax, e, M, small, signed)
        ax.set_xlabel(f"true {xlabel}"); ax.set_ylabel(f"reconstructed {xlabel}"); ax.set_title(title, fontsize=10)
        fig.colorbar(im, ax=ax, label="%" if not signed else "percentage points")
    fig.suptitle(f"P(reco bin | true bin) on the {name} grid: official MC vs VBLL {MODEL}", fontsize=11)
    fig.tight_layout(); fig.savefig(HERE / "figs" / f"{name}_migration_mc_vs_{TAG}.png", dpi=150, bbox_inches="tight"); plt.close(fig)
    diag = np.diag(D); off = D - np.diag(diag)
    out[name] = {"diagonal_difference_pp": diag.round(1).tolist(), "max_abs_diagonal_difference_pp": float(np.abs(diag).max()),
                 "max_abs_offdiagonal_difference_pp": float(np.abs(off).max()),
                 "mc_diagonal_percent": mc[name]["diagonal_percent"], f"{TAG}_diagonal_percent": nw[name]["diagonal_percent"]}
    print(f"{name}: diagonal difference ({TAG} - MC, pp) {diag.round(1).tolist()}; max |off-diagonal difference| {np.abs(off).max():.1f} pp")
json.dump(out, open(HERE / f"migration_mc_vs_{TAG}.json", "w"), indent=1)
print(f"wrote figs/<grid>_migration_mc_vs_{TAG}.png + migration_mc_vs_{TAG}.json")
