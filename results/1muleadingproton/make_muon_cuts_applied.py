#!/usr/bin/env python
"""results/1muleadingproton, section 2: GiBUU truth-signal muon momentum and angle with the selection cuts applied.

    python results/1muleadingproton/make_muon_cuts_applied.py        (default pixi environment)

The cuts are the channel's reconstruction-level selection `minerva_ccqelike_1mu1p_v0` (channels/minerva_me_ccqelike_1mu1p.yaml):
tracker fiducial vertex, MINOS-matched negative muon, dead time, muon window (17 deg, 2-20 GeV/c), a contained proton
candidate with dE/dx score > 0.35 inside 90 deg / 0.4-1.3 GeV/c, no Michel electron, at most one isolated blob. A generator
sample has no reconstruction, so the cuts enter as the selection efficiency learned from the official MC: per true cell
of a grid, eps = selected truth signal / truth signal in the fiducial (the `eff` of the grid's binned response,
surrogates/<channel>/<grid>/binned_FHC_1A-1P, built from the 12 ME FHC playlists). Each signal event is weighted by the
eps of its true cell on the grid of the plotted variable (muon_p, muon_theta, muon_costheta); the lower panel shows the
resulting efficiency per fine bin. Truth signal and weights as in section 1 (beam-frame angles, data exposure).
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
from ndp.channels import load_channel, load_measurement, observables as obs   # noqa: E402
from ndp.compare.folded import expected_true_cells                  # noqa: E402
from ndp.config import load_site_config                             # noqa: E402
from ndp.events import TruthTable                                   # noqa: E402
from ndp.pipeline import resolve_surrogate                          # noqa: E402
from ndp.products import total_pot                                  # noqa: E402
from scripts import plot_style as ps                                # noqa: E402

CHANNEL = "minerva_me_ccqelike_1mu1p"
GIBUU = PLATFORM / "runs/_generator_cache/gibuu_1f8bcb950be1dcd0/truth.npz"
cfg = load_site_config(); ch = load_channel(CHANNEL); t = TruthTable.load(GIBUU)
pot = total_pot(cfg, ch, "data")
scale = expected_true_cells(ch, load_measurement(ch, "muon_p"), t, pot)["scale"]
w = t["weight"] * scale
sig = ch.in_phase_space(t) & ch.is_signal(t)
ws = w[sig]
p = obs.lep_p(t, "beam")[sig]; th = obs.lep_theta_deg(t, "beam")[sig]; cth = np.cos(np.radians(th))
cos17 = float(np.cos(np.radians(17.0)))
figs = HERE / "figs"; figs.mkdir(exist_ok=True)
c_truth, c_cut, c_eff = ps.color_for("mc_signal"), ps.okabe_ito()["blue"], ps.color_for("efficiency")
out = {"selection": ch.selection["name"], "cuts": list(ch.selection["cuts"]), "n_signal": int(sig.sum()), "signal_events_at_data_pot": float(ws.sum()),
       "efficiency_source": "eff per true cell of the grid's binned response (surrogates/<channel>/<grid>/binned_FHC_1A-1P)"}
for name, grid, x, edges, xlabel in (("muon_p", "muon_p", p, np.linspace(2, 20, 37), "true muon momentum [GeV/c]"),
                                     ("muon_theta", "muon_theta", th, np.linspace(0, 17, 35), "true muon angle to the beam [deg]"),
                                     ("muon_costheta", "muon_costheta", cth, np.linspace(cos17, 1.0, 45), "true muon cos(theta) to the beam")):
    m = load_measurement(ch, grid); sur, sp = resolve_surrogate(ch, m, cfg)
    xt, yt = m.truth_observables(ch, t); g = m.binning.digitize(xt[sig], yt[sig])
    eff = np.where(g >= 0, sur.eff[np.where(g >= 0, g, 0)], 0.0)
    h_t, _ = np.histogram(x, bins=edges, weights=ws); h_c, _ = np.histogram(x, bins=edges, weights=ws * eff)
    with np.errstate(invalid="ignore", divide="ignore"):
        e_bin = np.where(h_t > 0, h_c / h_t, np.nan)
    fig, (ax, axe) = plt.subplots(2, 1, figsize=(7, 6.4), sharex=True, gridspec_kw={"height_ratios": [3, 1.3], "hspace": 0.06})
    ax.stairs(h_t, edges, color=c_truth, lw=1.6, baseline=None, label="signal definition only (section 1)")
    ax.stairs(h_c, edges, color=c_cut, lw=1.8, baseline=0, fill=True, alpha=0.25)
    ax.stairs(h_c, edges, color=c_cut, lw=1.8, baseline=None, label=f"selection cuts applied (x eps of the {grid} grid)")
    ax.set_ylabel("signal events at 1.057e21 POT (fiducial)"); ax.set_ylim(bottom=0); ax.grid(alpha=0.3); ax.legend(fontsize=9)
    axe.stairs(e_bin, edges, color=c_eff, lw=1.6, baseline=None)
    axe.set_ylabel("selection efficiency"); axe.set_ylim(0, max(0.6, np.nanmax(e_bin) * 1.15)); axe.grid(alpha=0.3); axe.set_xlabel(xlabel)
    fig.savefig(figs / f"{name}_cuts_applied.png", dpi=150, bbox_inches="tight"); plt.close(fig)
    c = 0.5 * (edges[:-1] + edges[1:]); cum = np.cumsum(h_c) / h_c.sum()
    out[name] = {"grid": grid, "grid_edges": list(m.x.edges), "grid_eff": sur.eff.tolist(), "surrogate": str(sp.relative_to(PLATFORM)),
                 "edges": edges.tolist(), "truth_signal": h_t.tolist(), "cuts_applied": h_c.tolist(), "efficiency_per_bin": e_bin.tolist(),
                 "selected_at_data_pot": float(h_c.sum()), "mean_efficiency": float(h_c.sum() / h_t.sum()),
                 "median_after_cuts": float(c[np.searchsorted(cum, 0.5)]), "mean_after_cuts": float(np.average(x, weights=ws * eff)),
                 "peak_bin_after_cuts": [float(edges[h_c.argmax()]), float(edges[h_c.argmax() + 1])],
                 "n_true_out_of_grid": int((g < 0).sum())}
    print(name, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in out[name].items() if k in ("selected_at_data_pot", "mean_efficiency", "median_after_cuts", "mean_after_cuts", "peak_bin_after_cuts", "n_true_out_of_grid")})
    print("   grid eff:", np.round(sur.eff, 3).tolist())
json.dump(out, open(HERE / "muon_cuts_applied.json", "w"), indent=1)
print("cuts:", out["cuts"])
