#!/usr/bin/env python
"""results/1muleadingproton, section 1: GiBUU muon momentum and angle after the signal definition only.

    python results/1muleadingproton/make_muon_truth_signal.py        (default pixi environment)

Sample: the GiBUU 2025 numu CC on 12C energy-scan sample under the NuMI ME FHC flux
(runs/_generator_cache/gibuu_1f8bcb950be1dcd0/truth.npz, beam frame). The only requirement is the channel's
truth signal definition (channels/minerva_me_ccqelike_1mu1p.yaml, arXiv:2503.15047): one muon with theta < 17 deg
and 2 < p < 20 GeV/c, at least one proton with theta < 70 deg and 0.5 < p < 1.1 GeV/c, no mesons, no baryons
heavier than the neutron, no photons above 10 MeV. No reconstruction cut, efficiency or smearing is applied.
Weights: the GiBUU event weights scaled to the data exposure (1.0574e21 POT, N_nucleons 3.23e30, the channel's
flux integral), so the histograms count expected signal events in the fiducial volume.
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
from ndp.channels import load_channel, observables as obs           # noqa: E402
from ndp.compare.folded import expected_true_cells                  # noqa: E402
from ndp.channels import load_measurement                           # noqa: E402
from ndp.config import load_site_config                             # noqa: E402
from ndp.events import TruthTable                                   # noqa: E402
from ndp.products import total_pot                                  # noqa: E402
from scripts import plot_style as ps                                # noqa: E402

CHANNEL = "minerva_me_ccqelike_1mu1p"
GIBUU = PLATFORM / "runs/_generator_cache/gibuu_1f8bcb950be1dcd0/truth.npz"
cfg = load_site_config(); ch = load_channel(CHANNEL); t = TruthTable.load(GIBUU)
pot = total_pot(cfg, ch, "data")
scale = expected_true_cells(ch, load_measurement(ch, "muon_p"), t, pot)["scale"]
w = t["weight"] * scale
sig = ch.in_phase_space(t) & ch.is_signal(t)
p = obs.lep_p(t, "beam")[sig]; th = obs.lep_theta_deg(t, "beam")[sig]; ws = w[sig]
figs = HERE / "figs"; figs.mkdir(exist_ok=True)
color = ps.color_for("mc_signal")
out = {"sample": str(GIBUU.relative_to(PLATFORM)), "n_generated_cc": int(t.n), "n_signal": int(sig.sum()), "pot_data": pot,
       "weight_scale": scale, "signal_events_at_data_pot": float(ws.sum()), "generated_cc_events_at_data_pot": float(w.sum())}
cos17 = float(np.cos(np.radians(17.0)))     # the signal window's lower edge in cos(theta): 0.956305
for name, x, edges, xlabel in (("muon_p", p, np.linspace(2, 20, 37), "true muon momentum [GeV/c]"),
                               ("muon_theta", th, np.linspace(0, 17, 35), "true muon angle to the beam [deg]"),
                               ("muon_costheta", np.cos(np.radians(th)), np.linspace(cos17, 1.0, 45), "true muon cos(theta) to the beam")):
    h, _ = np.histogram(x, bins=edges, weights=ws)
    fig, ax = plt.subplots(figsize=(7, 4.8))
    ax.stairs(h, edges, color=color, lw=1.8, baseline=0, fill=True, alpha=0.25)
    ax.stairs(h, edges, color=color, lw=1.8, baseline=None, label="GiBUU 2025, 1mu1p signal definition only")
    ax.set_xlabel(xlabel); ax.set_ylabel("signal events at 1.057e21 POT (fiducial)"); ax.set_ylim(bottom=0); ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
    fig.tight_layout(); fig.savefig(figs / f"{name}_truth_signal.png", dpi=150); plt.close(fig)
    c = 0.5 * (edges[:-1] + edges[1:]); cum = np.cumsum(h) / h.sum()
    out[name] = {"edges": edges.tolist(), "counts": h.tolist(), "median": float(c[np.searchsorted(cum, 0.5)]),
                 "mean": float(np.average(x, weights=ws)), "peak_bin": [float(edges[h.argmax()]), float(edges[h.argmax() + 1])]}
json.dump(out, open(HERE / "muon_truth_signal.json", "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if not isinstance(v, dict)}, indent=1))
for name in ("muon_p", "muon_theta"):
    print(name, {k: v for k, v in out[name].items() if k not in ("edges", "counts")})
