#!/usr/bin/env python
"""results/1muleadingproton, section 5: a VBLL surrogate trained in the platform with (E, px, py, pz, p, costheta) as
inputs and outputs, applied to GiBUU.

    pixi run -e ml python results/1muleadingproton/make_muon_vbll_fhc6.py

Model: surrogates/<channel>/_vbll/fhc6_het, trained by `ndp surrogate train-vbll` (ndp/surrogate/vbll_train.py) on the
selected signal pairs of the 12 ME FHC StandardMC playlists in the beam frame. This script
  1. summarises the training (curves, validation pulls / coverage / widths) from the model directory;
  2. builds the per-grid wrappers (muon_p, muon_theta, muon_costheta) for the new model and, where missing, for the
     ported x60_het model, with their closure on the official MC truth (`build_vbll`), and compares the probability of
     staying in the true bin: official MC vs x60_het vs fhc6_het;
  3. applies the new model to the GiBUU signal exactly as section 3 did (20 copies per event, conditioned on the
     reco windows, efficiency of section 2): distributions at the three stages and the migration matrices.
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
from ndp.compare.folded import expected_true_cells                  # noqa: E402
from ndp.config import load_site_config                             # noqa: E402
from ndp.events import TruthTable                                   # noqa: E402
from ndp.pipeline import resolve_surrogate                          # noqa: E402
from ndp.products import total_pot                                  # noqa: E402
from ndp.surrogate.vbll import build_vbll, truth_4vectors, synthetic_reco, reco_window_pass   # noqa: E402
from ndp.surrogate.vbll_model import load_vbll_model                # noqa: E402
from scripts import plot_style as ps                                # noqa: E402

CHANNEL = "minerva_me_ccqelike_1mu1p"
NEW = sys.argv[1] if len(sys.argv) > 1 else "fhc6_het"
OLD = "x60_het"
TAG = NEW.replace("_het", "")   # file-name tag: fhc6, fhc4, ...
GRIDS = ("muon_p", "muon_theta", "muon_costheta")
GIBUU = PLATFORM / "runs/_generator_cache/gibuu_1f8bcb950be1dcd0/truth.npz"
K, SEED = 20, 0
T0 = time.time()


def log(msg):
    print(f"[{TAG} {time.time() - T0:7.1f}s] {msg}", flush=True)


cfg = load_site_config(); ch = load_channel(CHANNEL); figs = HERE / "figs"
oi = ps.okabe_ito()
c_truth, c_cut, c_new, c_old, c_mc = ps.color_for("mc_signal"), oi["blue"], oi["vermilion"], oi["orange"], ps.color_for("data")
out = {"model": NEW}

# ---- 1. training summary ----------------------------------------------------------------------------------------
mdir = cfg.surrogates / CHANNEL / "_vbll" / NEW
js = json.loads((mdir / "surrogate.json").read_text()); spec = js["spec"]; tr = spec["meta"]["training"]; hist = js["training_history"]
out["training"] = {"inputs": spec["inputs"], "outputs": spec["outputs"], "frame": spec["frame"], "head": spec["head_type"],
                   "architecture": {k: spec[k] for k in ("d_embed", "hidden", "n_layers", "noise_prior_scale")},
                   "n_train_pairs": tr["n_train_pairs"], "n_val_pairs": tr["n_val_pairs"], "epochs_run": tr["epochs_run"], "best_val_nll": tr["best_val_nll"],
                   "train_s": tr["train_s"], "batch_size": tr["batch_size"], "validation_metrics": tr["validation_metrics"]}
fig, ax = plt.subplots(figsize=(6.5, 4))
ax.plot(hist["epoch"], hist["train_objective"], color=c_new, label="training objective (ELBO, per batch mean)")
ax.plot(hist["epoch"], hist["val_nll"], color=c_mc, label="validation predictive NLL")
ax.set_xlabel("epoch"); ax.set_ylabel("loss (normalised units)"); ax.grid(alpha=0.3); ax.legend(fontsize=9)
ax.set_title(f"{NEW}: training on {tr['n_train_pairs']} pairs, validation on {tr['n_val_pairs']}", fontsize=10)
fig.tight_layout(); fig.savefig(figs / f"{TAG}_training_curves.png", dpi=150); plt.close(fig)
log(f"training: {tr['epochs_run']} epochs, best val NLL {tr['best_val_nll']:.4f}")

# ---- 2. wrappers + closure on the official MC, three grids, both models ---------------------------------------------
ms = [load_measurement(ch, g) for g in GRIDS]
closures = {}
for name in (OLD, NEW):
    have = all((m.surrogate_root(cfg) / f"vbll_{name}" / "closure.json").exists() for m in ms)
    if not have:
        log(f"building wrappers + closure for {name} on {GRIDS} ...")
        build_vbll(ch, ms, cfg, model_name=name, n_samples=K, seed=SEED, truncate_to_reco_windows=True, log=log)
    closures[name] = {m.name: json.loads((m.surrogate_root(cfg) / f"vbll_{name}" / "closure.json").read_text()) for m in ms}
mc_mig = json.load(open(HERE / "muon_mc_migration.json"))["matrices"]
gib_mig_old = json.load(open(HERE / "muon_vbll_smeared.json"))["migration"]
out["closure"] = {}
for m in ms:
    row = {}
    for name in (OLD, NEW):
        c = closures[name][m.name]
        pred, reco = np.array(c["pred_cells"]), np.array(c["reco_cells"])
        row[name] = {"ratio_per_bin": (pred / reco).round(4).tolist(), "total_ratio": c["ratio_pred_over_reco"], "max_rel_dev_cells_gt_50": c["max_rel_dev_cells_gt_50"],
                     "n_bins_within_5pct": int(np.sum(np.abs(pred / reco - 1) <= 0.05)), "n_bins": int(len(reco))}
    out["closure"][m.name] = row
    log(f"closure {m.name}: {OLD} within 5 %: {row[OLD]['n_bins_within_5pct']}/{row[OLD]['n_bins']}, {NEW}: {row[NEW]['n_bins_within_5pct']}/{row[NEW]['n_bins']}")

# closure figure: per grid, fold/selected for both models
fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
for ax, m in zip(axes, ms):
    e = np.asarray(m.x.edges, float)
    ax.stairs(np.array(out["closure"][m.name][OLD]["ratio_per_bin"]), e, color=c_old, lw=1.6, baseline=None, label=f"{OLD} (ported, 4-vector in/out)")
    ax.stairs(np.array(out["closure"][m.name][NEW]["ratio_per_bin"]), e, color=c_new, lw=1.8, baseline=None, label=f"{NEW} (platform-trained, {len(spec['inputs'])} in / {len(spec['outputs'])} out)")
    ax.axhspan(0.95, 1.05, color=ps.color_for("reference"), alpha=0.12, lw=0); ax.axhline(1, color=ps.color_for("reference"), ls="--", lw=0.8)
    ax.set_xlabel(f"true {m.x.label}"); ax.set_ylabel("VBLL fold / selected signal (official MC)"); ax.set_ylim(0.4, 1.8); ax.grid(alpha=0.3); ax.set_title(f"{m.name} grid", fontsize=10)
axes[0].legend(fontsize=8)
fig.suptitle("Closure on the official MC: the surrogate fold of the fiducial truth signal over the selected signal, per reco bin", fontsize=11)
fig.tight_layout(); fig.savefig(figs / f"{TAG}_closure_three_grids.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ---- 3. GiBUU: three stages + migration with the new model -----------------------------------------------------------
t = TruthTable.load(GIBUU); pot = total_pot(cfg, ch, "data")
scale = expected_true_cells(ch, ms[0], t, pot)["scale"]; w = t["weight"] * scale
sig = ch.in_phase_space(t) & ch.is_signal(t); ws = w[sig]
p_t, th_t = obs.lep_p(t, "beam")[sig], obs.lep_theta_deg(t, "beam")[sig]; cth_t = np.cos(np.radians(th_t)); cos17 = float(np.cos(np.radians(17.0)))
model = load_vbll_model(mdir)
mu4, pr4 = truth_4vectors(ch, t, sig, model.spec.frame)
x_mu, x_pr = model.input_features(mu4), model.input_features(pr4)
mu_s = model.smear("muon", x_mu, n_samples=K, seed=SEED); pr_s = model.smear("proton", x_pr, n_samples=K, seed=SEED + 1)
recos = [synthetic_reco(mu_s[k], pr_s[k], model.spec.frame, model.spec.outputs) for k in range(K)]
pass_k = np.stack([reco_window_pass(ch, r) for r in recos]); n_pass = pass_k.sum(0); share = np.where(n_pass > 0, 1.0 / np.maximum(n_pass, 1), 0.0)
out["gibuu"] = {"n_signal": int(sig.sum()), "window_pass_fraction": float(pass_k.mean()), "n_events_all_copies_out_of_window": int((n_pass == 0).sum()), "variables": {}, "migration": {}}
log(f"GiBUU: {sig.sum()} signal events smeared, window pass fraction {pass_k.mean():.4f}")
specs = (("muon_p", p_t, lambda r: r["reco_p"], np.linspace(2, 20, 37), "muon momentum [GeV/c]"),
         ("muon_theta", th_t, lambda r: np.degrees(r["reco_theta"]), np.linspace(0, 17, 35), "muon angle to the beam [deg]"),
         ("muon_costheta", cth_t, lambda r: np.cos(r["reco_theta"]), np.linspace(cos17, 1.0, 45), "muon cos(theta) to the beam"))
old_json = json.load(open(HERE / "muon_vbll_smeared.json"))
for name, x_true, f_reco, edges, xlabel in specs:
    m = load_measurement(ch, name); sur, _ = resolve_surrogate(ch, m, cfg)
    xt, yt = m.truth_observables(ch, t); g = m.binning.digitize(xt[sig], yt[sig]); eff = np.where(g >= 0, sur.eff[np.where(g >= 0, g, 0)], 0.0)
    h_t = np.histogram(x_true, bins=edges, weights=ws)[0]; h_c = np.histogram(x_true, bins=edges, weights=ws * eff)[0]
    h_s = np.zeros(len(edges) - 1)
    for k, r in enumerate(recos):
        ok = pass_k[k]; h_s += np.histogram(f_reco(r)[ok], bins=edges, weights=(ws * eff * share)[ok])[0]
    h_old = np.array(old_json[name]["vbll_smeared_cuts_applied"])
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = np.where(h_c > 0, h_s / h_c, np.nan); ratio_old = np.where(h_c > 0, h_old / h_c, np.nan)
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(7, 6.4), sharex=True, gridspec_kw={"height_ratios": [3, 1.3], "hspace": 0.06})
    ax.stairs(h_t, edges, color=c_truth, lw=1.4, baseline=None, label="signal definition only (section 1, true)")
    ax.stairs(h_c, edges, color=c_cut, lw=1.6, baseline=None, label="selection cuts applied (section 2, true)")
    ax.stairs(h_old, edges, color=c_old, lw=1.3, ls="--", baseline=None, label=f"{OLD} smeared, cuts applied (section 3)")
    ax.stairs(h_s, edges, color=c_new, lw=1.8, baseline=0, fill=True, alpha=0.2)
    ax.stairs(h_s, edges, color=c_new, lw=1.8, baseline=None, label=f"{NEW} smeared, cuts applied (section 5, reconstructed)")
    ax.set_ylabel("signal events at 1.057e21 POT (fiducial)"); ax.set_ylim(bottom=0); ax.grid(alpha=0.3); ax.legend(fontsize=8)
    axr.stairs(ratio_old, edges, color=c_old, lw=1.3, ls="--", baseline=None, label=OLD); axr.stairs(ratio, edges, color=c_new, lw=1.6, baseline=None, label=NEW)
    axr.axhline(1.0, color=ps.color_for("reference"), ls="--", lw=0.8); axr.set_ylabel("smeared / unsmeared"); axr.set_ylim(0.0, max(2.0, float(np.nanmax(ratio_old)) * 1.1))
    axr.grid(alpha=0.3); axr.set_xlabel(xlabel); axr.legend(fontsize=8, loc="upper left")
    fig.savefig(figs / f"{name}_vbll_{TAG}_smeared.png", dpi=150, bbox_inches="tight"); plt.close(fig)
    c = 0.5 * (edges[:-1] + edges[1:]); cum = np.cumsum(h_s) / h_s.sum()
    out["gibuu"]["variables"][name] = {"edges": edges.tolist(), "truth_signal": h_t.tolist(), "cuts_applied": h_c.tolist(), f"{TAG}_smeared_cuts_applied": h_s.tolist(),
                                       "x60_smeared_cuts_applied": h_old.tolist(), "ratio_smeared_over_unsmeared": ratio.tolist(), "total_smeared": float(h_s.sum()),
                                       "median_smeared": float(c[np.searchsorted(cum, 0.5)]), "peak_bin_smeared": [float(edges[h_s.argmax()]), float(edges[h_s.argmax() + 1])],
                                       "ratio_range": [float(np.nanmin(ratio)), float(np.nanmax(ratio))]}
    # migration matrix (analysis grid + fine), same construction as section 3
    coarse = np.asarray(m.x.edges, float); mats = {}
    for tag, e in (("analysis_grid", coarse), ("fine", edges)):
        H = np.zeros((len(e) - 1, len(e) - 1)); lostw = np.zeros(len(e) - 1)
        for k, r in enumerate(recos):
            ok = pass_k[k]; v = f_reco(r)[ok]; xt_ok = x_true[ok]; wk = (ws * share)[ok]
            H += np.histogram2d(v, xt_ok, bins=[e, e], weights=wk)[0]
            inside = (v >= e[0]) & (v < e[-1]); lostw += np.histogram(xt_ok[~inside], bins=e, weights=wk[~inside])[0]
        col = H.sum(0) + lostw
        with np.errstate(invalid="ignore", divide="ignore"):
            mats[tag] = (np.where(col > 0, H / col, np.nan), H)
    Pc, Hc = mats["analysis_grid"]; diag = np.diag(Pc)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5.6))
    im1 = a1.pcolormesh(coarse, coarse, Pc * 100.0, cmap="Blues", vmin=0, vmax=100, shading="flat")
    for j in range(len(coarse) - 1):
        for i in range(len(coarse) - 1):
            if np.isfinite(Pc[i, j]) and Pc[i, j] >= 0.005:
                a1.text(0.5 * (coarse[j] + coarse[j + 1]), 0.5 * (coarse[i] + coarse[i + 1]), f"{100 * Pc[i, j]:.0f}", ha="center", va="center",
                        fontsize=7 if len(coarse) > 10 else 8, color="white" if Pc[i, j] > 0.5 else "black")
    a1.set_xlabel(f"true {xlabel}"); a1.set_ylabel(f"reconstructed ({NEW}) {xlabel}"); a1.set_title(f"P(reco bin | true bin) [%] on the {name} grid", fontsize=10)
    fig.colorbar(im1, ax=a1, label="%")
    im2 = a2.pcolormesh(edges, edges, mats["fine"][0] * 100.0, cmap="Blues", vmin=0, shading="flat")
    a2.plot([edges[0], edges[-1]], [edges[0], edges[-1]], color=ps.color_for("reference"), lw=0.8, ls="--")
    a2.set_xlabel(f"true {xlabel}"); a2.set_ylabel(f"reconstructed ({NEW}) {xlabel}"); a2.set_title("same, fine bins", fontsize=10); fig.colorbar(im2, ax=a2, label="%")
    fig.suptitle(f"VBLL {NEW} migration of the selected GiBUU signal: {xlabel}", fontsize=11)
    fig.tight_layout(); fig.savefig(figs / f"{name}_vbll_{TAG}_migration.png", dpi=150, bbox_inches="tight"); plt.close(fig)
    out["gibuu"]["migration"][name] = {"grid_edges": coarse.tolist(), "P_reco_given_true_percent_reco_by_true": (np.nan_to_num(Pc) * 100).round(2).tolist(),
                                       "diagonal_percent": (np.nan_to_num(diag) * 100).round(1).tolist(),
                                       "weighted_mean_diagonal_percent": float(100 * np.nansum(diag * Hc.sum(0)) / Hc.sum()),
                                       "mc_diagonal_percent": mc_mig[name]["diagonal_percent"], "mc_weighted_mean_diagonal_percent": mc_mig[name]["weighted_mean_diagonal_percent"],
                                       "x60_diagonal_percent": gib_mig_old[name]["diagonal_percent"], "x60_weighted_mean_diagonal_percent": gib_mig_old[name]["weighted_mean_diagonal_percent"]}
    log(f"{name}: smeared total {h_s.sum():.0f}, ratio range {np.nanmin(ratio):.2f}-{np.nanmax(ratio):.2f}; diagonal % MC {mc_mig[name]['weighted_mean_diagonal_percent']:.1f} / "
        f"{OLD} {gib_mig_old[name]['weighted_mean_diagonal_percent']:.1f} / {NEW} {out['gibuu']['migration'][name]['weighted_mean_diagonal_percent']:.1f}")

# diagonal comparison: MC vs x60_het vs fhc6_het
fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
for ax, (name, _, _, _, xlabel) in zip(axes, specs):
    mg = out["gibuu"]["migration"][name]; e = np.array(mg["grid_edges"])
    ax.stairs(np.array(mg["mc_diagonal_percent"]), e, color=c_mc, lw=1.8, baseline=None, label="official MC")
    ax.stairs(np.array(mg["x60_diagonal_percent"]), e, color=c_old, lw=1.4, ls="--", baseline=None, label=f"VBLL {OLD}")
    ax.stairs(np.array(mg["diagonal_percent"]), e, color=c_new, lw=1.8, baseline=None, label=f"VBLL {NEW}")
    ax.set_xlabel(f"true {xlabel}"); ax.set_ylabel("P(same bin) [%]"); ax.set_ylim(0, 100); ax.grid(alpha=0.3); ax.set_title(f"{name} grid", fontsize=10)
axes[0].legend(fontsize=8)
fig.suptitle("Probability of reconstructing in the true bin: official MC vs the two VBLL surrogates (selected 1mu1p signal, GiBUU for the surrogates)", fontsize=10)
fig.tight_layout(); fig.savefig(figs / f"migration_diagonal_mc_vs_vbll_{TAG}.png", dpi=150, bbox_inches="tight"); plt.close(fig)
json.dump(out, open(HERE / f"muon_vbll_{TAG}.json", "w"), indent=1)
log(f"wrote muon_vbll_{TAG}.json + figures")
