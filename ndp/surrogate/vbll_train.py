"""Train a VBLL detector surrogate on the platform's own paired MC (ml environment: torch + vbll).

    pixi run -e ml python -m ndp surrogate train-vbll --channel minerva_me_ccqelike_1mu1p --name fhc6_het \\
        --inputs E,px,py,pz,p,costheta --outputs E,px,py,pz,p,costheta --frame beam

Pairs: for every reconstructed candidate passing the channel selection whose truth is signal inside the fiducial
volume (the binned responses' numerator, all playlist products), the true and reconstructed 4-vectors of the muon
and of the leading proton in the signal window, in `frame`, MeV — the same pairing as the platform's binned
responses (truth proton = leading in the window, chosen in the channel frame). Features (inputs and outputs) are
any of E, px, py, pz, p = |p|, costheta = pz/p in that frame (`vbll_model.FEATURES`); the reconstructed p and
costheta come from the cached `reco_p` / `reco_theta` (muon) and `reco_proton_p` / `reco_proton_theta_beam` (proton).

Model and objective: the VBLL_SurrogateModel architecture (embedding + MLP + one vbll head per particle + a
particle-type classification head weighted by lambda_cls) trained with Adam on the VBLL ELBO. For the heteroscedastic
head the training loss carries the shape fix of VBLL_SurrogateModel/code/vbll_patches.py (vbll 0.4.9's
`HetRegression._get_train_loss_fn` sums the trace term over the wrong axis); the predictive is vbll's own.
Early stopping on the validation predictive NLL (vbll's val_loss_fn, 20 noise samples). Every seed is pinned.
Saved in the platform's model format (arrays.npz + surrogate.json) under surrogates/<channel>/_vbll/<name>/.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from ..channels.observables import native_frame, frame_rotation_angle, rotate_about_x
from ..channels.signal import leading_proton
from ..events import TruthTable
from ..io import cheap_fingerprint, timestamp, git_state
from .vbll_model import VBLLSpec, VBLLSurrogateModel, build_module, features_from_4vectors, PARTICLES, FEATURES, PARTICLE_INDEX

MEV = 1e3


# ---- pairs from the playlist products ---------------------------------------------------------------------------
def collect_pairs(cfg, channel, frame: str = "beam", log=print) -> dict:
    """True and reconstructed (E, px, py, pz) [MeV, `frame`] of the muon and the leading proton for every selected
    signal candidate in the fiducial volume, one playlist at a time. Returns {"truth_muon", "reco_muon",
    "truth_proton", "reco_proton": (n, 4), "playlist": (n,) labels, "inputs": [...], "frame"}."""
    from ..channels.selections import select
    from ..products import iter_reco_chunks, playlist_dir, has_products
    from ..adapters.minerva_anatuple import cache_tag
    T, R, TP, RP, labels, inputs = [], [], [], [], [], []
    t0 = time.time()
    for label, r, pot, src in iter_reco_chunks(cfg, channel, "mc"):
        if has_products(channel):
            beam, pl = label.split("/", 1)
            rt_path = playlist_dir(cfg, channel, beam, pl) / f"reco_{pl}_truthcols_skim.npz"
        else:
            rt_path = cfg.require("data_dir") / "cache" / f"reco_{label}_truthcols.npz"
        rt = TruthTable.load(rt_path)
        keep = select(channel, r) & channel.is_signal(rt) & channel.in_phase_space(rt)
        ang_t = frame_rotation_angle(native_frame(rt), frame)
        px, py, pz = rotate_about_x(rt["lep_px"][keep], rt["lep_py"][keep], rt["lep_pz"][keep], ang_t)
        T.append(np.stack([rt["lep_E"][keep], px, py, pz], 1) * MEV)
        lp = leading_proton(rt, channel.frame, (channel.signal or {}).get("proton"))
        ang_p = frame_rotation_angle(channel.frame, frame)
        ppx, ppy, ppz = rotate_about_x(lp["px"][keep], lp["py"][keep], lp["pz"][keep], ang_p)
        TP.append(np.stack([lp["E"][keep], ppx, ppy, ppz], 1) * MEV)
        ang_r = frame_rotation_angle("detector", frame)                    # cached components are detector-frame
        mx, my, mz = rotate_about_x(np.asarray(r["reco_mu_px"], float)[keep], np.asarray(r["reco_mu_py"], float)[keep], np.asarray(r["reco_mu_pz"], float)[keep], ang_r)
        R.append(np.stack([np.asarray(r["reco_E_mu"], float)[keep], mx, my, mz], 1) * MEV)
        qx, qy, qz = rotate_about_x(np.asarray(r["reco_proton_px"], float)[keep], np.asarray(r["reco_proton_py"], float)[keep], np.asarray(r["reco_proton_pz"], float)[keep], ang_r)
        RP.append(np.stack([np.asarray(r["reco_proton_E"], float)[keep], qx, qy, qz], 1) * MEV)
        labels.append(np.full(int(keep.sum()), label)); inputs.append({"input": label, "source": src, "n_pairs": int(keep.sum())})
        log(f"  [{time.time() - t0:6.1f}s] {label}: {int(keep.sum())} pairs")
        del r, rt
    out = {"truth_muon": np.concatenate(T).astype(np.float32), "reco_muon": np.concatenate(R).astype(np.float32),
           "truth_proton": np.concatenate(TP).astype(np.float32), "reco_proton": np.concatenate(RP).astype(np.float32),
           "playlist": np.concatenate(labels), "frame": frame, "inputs": inputs, "channel": channel.name, "selection": channel.selection.get("name")}
    bad = ~np.isfinite(np.concatenate([out[k] for k in ("truth_muon", "reco_muon", "truth_proton", "reco_proton")], axis=1)).all(axis=1)
    if bad.any():
        log(f"  dropping {int(bad.sum())} pairs with a non-finite component")
        for k in ("truth_muon", "reco_muon", "truth_proton", "reco_proton", "playlist"):
            out[k] = out[k][~bad]
    out["n_pairs"] = int(len(out["truth_muon"]))
    return out


def save_pairs(pairs: dict, path: Path) -> Path:
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    arrays = {k: pairs[k] for k in ("truth_muon", "reco_muon", "truth_proton", "reco_proton", "playlist")}
    meta = {k: v for k, v in pairs.items() if k not in arrays}
    np.savez_compressed(path, __meta__=json.dumps(meta), **arrays)
    return path


def load_pairs(path: Path) -> dict:
    z = np.load(path, allow_pickle=False)
    out = {k: z[k] for k in z.files if k != "__meta__"}
    out.update(json.loads(str(z["__meta__"])))
    return out


# ---- the patched heteroscedastic training loss (VBLL_SurrogateModel/code/vbll_patches.py, vbll 0.4.9) ----------------
def het_train_loss(head, x, y):
    """-ELBO of a vbll.HetRegression head for features x and targets y, with the trace term summed over the
    output dimension (missing in vbll 0.4.9) and grad_correction reduced to (batch,)."""
    import torch
    from vbll.layers.regression import expected_gaussian_kl, gaussian_kl
    W, M = head.W, head.M
    log_noise_cov = head.log_noise(x, M)
    expect_sigma_inv = torch.exp(-log_noise_cov.mean + 0.5 * log_noise_cov.scale ** 2)
    expect_log_sigma = log_noise_cov.mean
    grad_correction = ((expect_sigma_inv.detach()) ** head.grad_correction).mean(-1)
    err = y - (W.mean @ x[..., None]).squeeze(-1)
    mse_term = (err.pow(2) * expect_sigma_inv).sum(-1)
    logdet_term = expect_log_sigma.sum(-1)
    trace_term = W.covariance_weighted_inner_prod(x.unsqueeze(-2)[..., None]).sum(-1)
    total_elbo = -0.5 * torch.mean(grad_correction * (mse_term + logdet_term + trace_term))
    kl_term_ll = torch.mean(grad_correction * expected_gaussian_kl(W, head.prior_scale, expect_sigma_inv))
    kl_term_noise = torch.mean(grad_correction * gaussian_kl(M, head.noise_prior_scale))
    total_elbo -= head.regularization_weight * (kl_term_ll + kl_term_noise)
    return -total_elbo


# ---- training ------------------------------------------------------------------------------------------------------
def _normaliser(x: dict, y: dict, names_in, names_out) -> dict:
    norm = {"input": {}, "output": {}}
    for p in PARTICLES:
        norm["input"][p] = {n: [float(x[p][:, i].mean()), float(x[p][:, i].std())] for i, n in enumerate(names_in)}
        norm["output"][p] = {n: [float(y[p][:, i].mean()), float(y[p][:, i].std())] for i, n in enumerate(names_out)}
    return norm


def train_vbll(pairs: dict, *, inputs=FEATURES, outputs=FEATURES, head_type: str = "het", d_embed: int = 8, hidden: int = 64,
               n_layers: int = 3, wishart_scale: float = 1.0, prior_scale: float = 1.0, dof: float = 1.0, noise_prior_scale: float = 0.01,
               lr: float = 1e-3, epochs: int = 80, patience: int = 10, batch_size: int = 512, lambda_cls: float = 0.01,
               val_fraction: float = 0.2, seed: int = 42, threads: int | None = None, log=print) -> tuple[VBLLSpec, dict, dict]:
    """Returns (spec with the train-split normaliser, best state_dict, history + validation metrics)."""
    import torch
    if threads:
        torch.set_num_threads(int(threads))
    torch.manual_seed(seed); rng = np.random.default_rng(seed)
    inputs, outputs = tuple(inputs), tuple(outputs)
    n = pairs["n_pairs"]
    perm = rng.permutation(n); n_val = int(round(val_fraction * n)); val_idx, tr_idx = perm[:n_val], perm[n_val:]
    X = {"muon": features_from_4vectors(pairs["truth_muon"], inputs), "proton": features_from_4vectors(pairs["truth_proton"], inputs)}
    Y = {"muon": features_from_4vectors(pairs["reco_muon"], outputs), "proton": features_from_4vectors(pairs["reco_proton"], outputs)}
    norm = _normaliser({p: X[p][tr_idx] for p in PARTICLES}, {p: Y[p][tr_idx] for p in PARTICLES}, inputs, outputs)
    spec = VBLLSpec(d_embed=d_embed, hidden=hidden, n_layers=n_layers, head_type=head_type, parameterization="diagonal",
                    wishart_scale=wishart_scale, prior_scale=prior_scale, dof=dof, noise_prior_scale=noise_prior_scale,
                    n_train_per_particle=int(len(tr_idx)), inputs=inputs, outputs=outputs, units="MeV", frame=pairs["frame"], normaliser=norm)
    model = VBLLSurrogateModel(spec, build_module(spec))

    def tensors(idx):
        xs, ys, ts = [], [], []
        for p in PARTICLES:
            xs.append(model.normalise_input(p, X[p][idx])); ys.append((Y[p][idx] - model._stats("output", p)[0]) / (model._stats("output", p)[1] + np.float32(1e-8)))
            ts.append(np.full(len(idx), PARTICLE_INDEX[p], np.int64))
        return (torch.from_numpy(np.concatenate(xs).astype(np.float32)), torch.from_numpy(np.concatenate(ys).astype(np.float32)),
                torch.from_numpy(np.concatenate(ts)))

    x_tr, y_tr, t_tr = tensors(tr_idx); x_va, y_va, t_va = tensors(val_idx)
    mod = model.module.train()
    opt = torch.optim.Adam(mod.parameters(), lr=lr)
    g = torch.Generator().manual_seed(seed)

    def objective(xb, yb, tb):
        h = mod.features(tb, xb)
        loss = lambda_cls * mod.cls_head(h).train_loss_fn(tb)
        for p, pidx in PARTICLE_INDEX.items():
            m = tb == pidx
            if not m.any():
                continue
            head = mod.reg_heads[p]
            loss = loss + (het_train_loss(head, h[m], yb[m]) if head_type == "het" else head(h[m]).train_loss_fn(yb[m]))
        return loss

    @torch.no_grad()
    def val_nll(xb, yb, tb, chunk=65536):
        mod.eval(); tot, cnt = 0.0, 0
        for lo in range(0, len(xb), chunk):
            h = mod.features(tb[lo:lo + chunk], xb[lo:lo + chunk]); yy = yb[lo:lo + chunk]; tt = tb[lo:lo + chunk]
            for p, pidx in PARTICLE_INDEX.items():
                m = tt == pidx
                if m.any():
                    tot += float(mod.reg_heads[p](h[m]).val_loss_fn(yy[m])) * int(m.sum()); cnt += int(m.sum())
        mod.train(); return tot / max(cnt, 1)

    history = {"epoch": [], "train_objective": [], "val_nll": [], "epoch_s": []}
    best, best_state, bad, n_tr = float("inf"), None, 0, len(x_tr)
    t0 = time.time()
    for ep in range(1, epochs + 1):
        te = time.time(); order = torch.randperm(n_tr, generator=g); run, nb = 0.0, 0
        for lo in range(0, n_tr, batch_size):
            b = order[lo:lo + batch_size]
            opt.zero_grad(); loss = objective(x_tr[b], y_tr[b], t_tr[b]); loss.backward(); opt.step()
            run += float(loss); nb += 1
        v = val_nll(x_va, y_va, t_va)
        history["epoch"].append(ep); history["train_objective"].append(run / nb); history["val_nll"].append(v); history["epoch_s"].append(round(time.time() - te, 1))
        improved = v < best - 1e-4
        if improved:
            best, bad = v, 0; best_state = {k: t.detach().clone() for k, t in mod.state_dict().items()}
        else:
            bad += 1
        log(f"epoch {ep:3d} | train objective {run / nb:.4f} | val NLL {v:.4f}{' *' if improved else ''} | {time.time() - te:.0f} s")
        if bad >= patience:
            log(f"early stopping at epoch {ep} (best val NLL {best:.4f})"); break
    mod.load_state_dict(best_state); mod.eval()
    metrics = validation_metrics(model, X, Y, val_idx, seed=seed)
    info = {"history": history, "best_val_nll": best, "epochs_run": len(history["epoch"]), "train_s": round(time.time() - t0, 1),
            "n_train_pairs": int(len(tr_idx)), "n_val_pairs": int(n_val), "seed": seed, "batch_size": batch_size, "lr": lr, "lambda_cls": lambda_cls,
            "patience": patience, "validation_metrics": metrics}
    return spec, best_state, info


@np.errstate(invalid="ignore", divide="ignore")
def validation_metrics(model: VBLLSurrogateModel, X: dict, Y: dict, idx, seed: int = 0) -> dict:
    """evaluate.py's numbers on the validation split, in normalised units: pull mean/std, sigma CV, coverage."""
    import torch
    from scipy import stats as spstats
    out = {}
    for p in PARTICLES:
        x = model.normalise_input(p, X[p][idx]); mu, sd = model._stats("output", p); y = (Y[p][idx] - mu) / (sd + np.float32(1e-8))
        torch.manual_seed(seed)
        mean, sig = model.predictive_normalised(p, x)
        r = mean - y; pull = r / (sig + 1e-8)
        out[p] = {"n": int(len(x)),
                  "pull_mean_std": {n: [float(pull[:, i].mean()), float(pull[:, i].std())] for i, n in enumerate(model.spec.outputs)},
                  "sigma_cv": {n: float(sig[:, i].std() / (sig[:, i].mean() + 1e-8)) for i, n in enumerate(model.spec.outputs)},
                  "coverage": {f"{a:.2f}": float(np.mean(np.abs(r) < spstats.norm.ppf((1 + a) / 2) * sig)) for a in (0.68, 0.90, 0.95)},
                  "residual_robust_sigma_physical": {n: float(1.4826 * np.median(np.abs(r[:, i] - np.median(r[:, i]))) * (sd[i] + 1e-8)) for i, n in enumerate(model.spec.outputs)},
                  "predicted_sigma_median_physical": {n: float(np.median(sig[:, i]) * (sd[i] + 1e-8)) for i, n in enumerate(model.spec.outputs)}}
    return out


def save_model(out_dir: Path, spec: VBLLSpec, state: dict, info: dict, pairs_path: Path | None, extra_meta: dict | None = None) -> Path:
    import torch, vbll  # noqa: F401
    spec.meta.update({"source": "trained in ndp-platform (ndp/surrogate/vbll_train.py) on the channel's selected signal pairs",
                      "trained": timestamp(), "platform_git": git_state(Path(__file__).resolve().parents[2]),
                      "training_pairs": {"path": str(pairs_path) if pairs_path else None, "fingerprint": cheap_fingerprint(pairs_path) if pairs_path else None,
                                         "n_train": info["n_train_pairs"], "n_val": info["n_val_pairs"]},
                      "training": {k: v for k, v in info.items() if k != "history"},
                      "versions": {"torch": torch.__version__, "vbll": "0.4.9"}, **(extra_meta or {})})
    out = VBLLSurrogateModel.save_dir(out_dir, spec, state, extra_json={"training_history": info["history"]})
    return out
