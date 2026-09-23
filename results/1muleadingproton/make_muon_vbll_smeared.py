#!/usr/bin/env python
"""results/1muleadingproton, section 3: the VBLL surrogate smearing applied to the selected GiBUU signal.

    pixi run -e ml python results/1muleadingproton/make_muon_vbll_smeared.py

Surrogate: the ported VBLL_SurrogateModel checkpoint x60_het (surrogates/<channel>/_vbll/x60_het, plan
runs/2026-09-22_vbll_surrogate_1mu1p). For every signal event of section 1 the true muon and leading-proton 4-vectors
are rotated into the model's detector frame and K = 20 reconstructed copies are drawn from the model's predictive
distribution (ndp/surrogate/vbll.py); the copies are conditioned on the selection's kinematic windows (muon 17 deg /
2-20 GeV/c, proton 90 deg / 0.4-1.3 GeV/c), the event's weight being shared among the passing copies; each event carries
the selection efficiency of its true cell as in section 2. The plotted quantities are the smeared (reconstructed)
muon momentum and beam angle, so this is the prediction that is compared with data in reconstructed space.
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
from ndp.surrogate.vbll import truth_4vectors, synthetic_reco, reco_window_pass   # noqa: E402
from ndp.surrogate.vbll_model import load_vbll_model                # noqa: E402
from scripts import plot_style as ps                                # noqa: E402

CHANNEL = "minerva_me_ccqelike_1mu1p"
GIBUU = PLATFORM / "runs/_generator_cache/gibuu_1f8bcb950be1dcd0/truth.npz"
K, SEED = 20, 0
cfg = load_site_config(); ch = load_channel(CHANNEL); t = TruthTable.load(GIBUU)
pot = total_pot(cfg, ch, "data")
scale = expected_true_cells(ch, load_measurement(ch, "muon_p"), t, pot)["scale"]
w = t["weight"] * scale
sig = ch.in_phase_space(t) & ch.is_signal(t)
ws = w[sig]
p_t = obs.lep_p(t, "beam")[sig]; th_t = obs.lep_theta_deg(t, "beam")[sig]; cth_t = np.cos(np.radians(th_t))
cos17 = float(np.cos(np.radians(17.0)))

model = load_vbll_model(cfg.surrogates / CHANNEL / "_vbll" / "x60_het")
mu, pr = truth_4vectors(ch, t, sig, model.spec.frame)
mu_s = model.smear("muon", mu, n_samples=K, seed=SEED); pr_s = model.smear("proton", pr, n_samples=K, seed=SEED + 1)
recos = [synthetic_reco(mu_s[k], pr_s[k], model.spec.frame) for k in range(K)]
pass_k = np.stack([reco_window_pass(ch, r) for r in recos]); n_pass = pass_k.sum(0)
share = np.where(n_pass > 0, 1.0 / np.maximum(n_pass, 1), 0.0)

figs = HERE / "figs"; figs.mkdir(exist_ok=True)
c_truth, c_cut, c_sm, c_ratio = ps.color_for("mc_signal"), ps.okabe_ito()["blue"], ps.okabe_ito()["vermilion"], ps.color_for("efficiency")
out = {"model": str((cfg.surrogates / CHANNEL / "_vbll" / "x60_het").relative_to(PLATFORM)), "model_frame": model.spec.frame, "n_samples": K, "seed": SEED,
       "n_signal": int(sig.sum()), "signal_events_at_data_pot": float(ws.sum()), "window_pass_fraction": float(pass_k.mean()),
       "n_events_all_copies_out_of_window": int((n_pass == 0).sum())}
specs = (("muon_p", "muon_p", p_t, "reco_p", lambda r: r["reco_p"], np.linspace(2, 20, 37), "muon momentum [GeV/c]"),
         ("muon_theta", "muon_theta", th_t, "reco_theta", lambda r: np.degrees(r["reco_theta"]), np.linspace(0, 17, 35), "muon angle to the beam [deg]"),
         ("muon_costheta", "muon_costheta", cth_t, "reco_theta", lambda r: np.cos(r["reco_theta"]), np.linspace(cos17, 1.0, 45), "muon cos(theta) to the beam"))
for name, grid, x_true, _, f_reco, edges, xlabel in specs:
    m = load_measurement(ch, grid); sur, sp = resolve_surrogate(ch, m, cfg)
    xt, yt = m.truth_observables(ch, t); g = m.binning.digitize(xt[sig], yt[sig])
    eff = np.where(g >= 0, sur.eff[np.where(g >= 0, g, 0)], 0.0)
    h_t, _ = np.histogram(x_true, bins=edges, weights=ws)
    h_c, _ = np.histogram(x_true, bins=edges, weights=ws * eff)
    h_s = np.zeros(len(edges) - 1); h_s_noeff = np.zeros(len(edges) - 1)
    for k, r in enumerate(recos):
        ok = pass_k[k]; v = f_reco(r)
        h_s += np.histogram(v[ok], bins=edges, weights=(ws * eff * share)[ok])[0]
        h_s_noeff += np.histogram(v[ok], bins=edges, weights=(ws * share)[ok])[0]
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = np.where(h_c > 0, h_s / h_c, np.nan)
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(7, 6.4), sharex=True, gridspec_kw={"height_ratios": [3, 1.3], "hspace": 0.06})
    ax.stairs(h_t, edges, color=c_truth, lw=1.4, baseline=None, label="signal definition only (section 1, true kinematics)")
    ax.stairs(h_c, edges, color=c_cut, lw=1.6, baseline=None, label="selection cuts applied (section 2, true kinematics)")
    ax.stairs(h_s, edges, color=c_sm, lw=1.8, baseline=0, fill=True, alpha=0.2)
    ax.stairs(h_s, edges, color=c_sm, lw=1.8, baseline=None, label="VBLL smeared, cuts applied (section 3, reconstructed kinematics)")
    ax.set_ylabel("signal events at 1.057e21 POT (fiducial)"); ax.set_ylim(bottom=0); ax.grid(alpha=0.3); ax.legend(fontsize=8)
    axr.stairs(ratio, edges, color=c_ratio, lw=1.6, baseline=None); axr.axhline(1.0, color=ps.color_for("reference"), ls="--", lw=0.8)
    axr.set_ylabel("smeared / unsmeared"); axr.set_ylim(0.0, max(2.0, float(np.nanmax(ratio)) * 1.1)); axr.grid(alpha=0.3); axr.set_xlabel(xlabel)
    fig.savefig(figs / f"{name}_vbll_smeared.png", dpi=150, bbox_inches="tight"); plt.close(fig)
    c = 0.5 * (edges[:-1] + edges[1:]); cum = np.cumsum(h_s) / h_s.sum()
    out[name] = {"grid": grid, "edges": edges.tolist(), "truth_signal": h_t.tolist(), "cuts_applied": h_c.tolist(), "vbll_smeared_cuts_applied": h_s.tolist(),
                 "vbll_smeared_no_efficiency": h_s_noeff.tolist(), "ratio_smeared_over_unsmeared": ratio.tolist(),
                 "total_cuts_applied": float(h_c.sum()), "total_smeared": float(h_s.sum()), "total_smeared_no_efficiency": float(h_s_noeff.sum()),
                 "median_smeared": float(c[np.searchsorted(cum, 0.5)]), "peak_bin_smeared": [float(edges[h_s.argmax()]), float(edges[h_s.argmax() + 1])],
                 "max_ratio": float(np.nanmax(ratio)), "min_ratio": float(np.nanmin(ratio))}
    print(name, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in out[name].items() if k.startswith(("total", "median", "peak", "max_", "min_"))})
    print("   ratio per bin:", np.round(ratio, 2).tolist())

# ---- migration matrices: true vs smeared (reconstructed), from the same copies -------------------------------
# Column-normalised: P(reco bin | true bin), each true column summing to 1 over the reco bins inside the range
# (copies smeared outside the range are the column's loss and are reported). Weights: w x share of the passing
# copies, as above; the selection efficiency is constant within a true cell of the analysis grid and drops out of
# the normalisation there. Left: the analysis grid of the measurement (what the folding uses); right: fine bins.
out["migration"] = {}
for name, grid, x_true, _, f_reco, edges, xlabel in specs:
    m = load_measurement(ch, grid)
    coarse = np.asarray(m.x.edges, float)
    mats = {}
    for tag, e in (("analysis_grid", coarse), ("fine", edges)):
        H = np.zeros((len(e) - 1, len(e) - 1))          # [reco, true]
        lost = np.zeros(len(e) - 1)
        for k, r in enumerate(recos):
            ok = pass_k[k]; v = f_reco(r)[ok]; xt_ok = x_true[ok]; wk = (ws * share)[ok]
            H += np.histogram2d(v, xt_ok, bins=[e, e], weights=wk)[0]
            inside = (v >= e[0]) & (v < e[-1])
            lost += np.histogram(xt_ok[~inside], bins=e, weights=wk[~inside])[0]
        col = H.sum(axis=0) + lost
        with np.errstate(invalid="ignore", divide="ignore"):
            P = np.where(col > 0, H / col, np.nan)
        mats[tag] = {"edges": e.tolist(), "counts_reco_by_true": H, "P_reco_given_true": P, "loss_fraction_per_true_bin": np.where(col > 0, lost / col, np.nan)}
    Pc = mats["analysis_grid"]["P_reco_given_true"]
    diag = np.diag(Pc)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5.6))
    im1 = a1.pcolormesh(coarse, coarse, Pc * 100.0, cmap="Blues", vmin=0, vmax=100, shading="flat")
    for j in range(len(coarse) - 1):
        for i in range(len(coarse) - 1):
            if np.isfinite(Pc[i, j]) and Pc[i, j] >= 0.005:
                a1.text(0.5 * (coarse[j] + coarse[j + 1]), 0.5 * (coarse[i] + coarse[i + 1]), f"{100 * Pc[i, j]:.0f}", ha="center", va="center",
                        fontsize=7 if len(coarse) > 10 else 8, color="white" if Pc[i, j] > 0.5 else "black")
    a1.set_xlabel(f"true {xlabel}"); a1.set_ylabel(f"reconstructed (VBLL) {xlabel}")
    a1.set_title(f"P(reco bin | true bin) [%] on the {grid} grid", fontsize=10)
    fig.colorbar(im1, ax=a1, label="%")
    Pf = mats["fine"]["P_reco_given_true"]
    im2 = a2.pcolormesh(edges, edges, Pf * 100.0, cmap="Blues", vmin=0, shading="flat")
    a2.plot([edges[0], edges[-1]], [edges[0], edges[-1]], color=ps.color_for("reference"), lw=0.8, ls="--")
    a2.set_xlabel(f"true {xlabel}"); a2.set_ylabel(f"reconstructed (VBLL) {xlabel}"); a2.set_title("same, fine bins", fontsize=10)
    fig.colorbar(im2, ax=a2, label="%")
    fig.suptitle(f"VBLL x60_het migration of the selected GiBUU signal: {xlabel}", fontsize=11)
    fig.tight_layout(); fig.savefig(figs / f"{name}_vbll_migration.png", dpi=150, bbox_inches="tight"); plt.close(fig)
    out["migration"][name] = {"grid": grid, "grid_edges": coarse.tolist(), "P_reco_given_true_percent_reco_by_true": (np.nan_to_num(Pc) * 100).round(2).tolist(),
                              "diagonal_percent": (np.nan_to_num(diag) * 100).round(1).tolist(),
                              "loss_fraction_per_true_bin": np.nan_to_num(mats["analysis_grid"]["loss_fraction_per_true_bin"]).round(4).tolist(),
                              "weighted_mean_diagonal_percent": float(100 * np.nansum(diag * mats["analysis_grid"]["counts_reco_by_true"].sum(0)) / mats["analysis_grid"]["counts_reco_by_true"].sum())}
    print(f"{name} migration on {grid}: diagonal % {out['migration'][name]['diagonal_percent']}, weighted mean {out['migration'][name]['weighted_mean_diagonal_percent']:.1f} %, "
          f"loss/true bin {out['migration'][name]['loss_fraction_per_true_bin']}")

json.dump(out, open(HERE / "muon_vbll_smeared.json", "w"), indent=1)
print({k: v for k, v in out.items() if not isinstance(v, dict)})
