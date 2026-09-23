"""Event-level detector surrogate from a ported VBLL checkpoint (kind `vbll_event`).

    reco cell  =  sum over signal events of  w · eff(true cell) · 1[reco cell of the VBLL-smeared event]   (+ background)

One smearing of the true (muon, leading-proton) 4-vectors serves every measurement grid: the K
reconstructed 4-vectors drawn per signal event are turned into the cached reco columns (`reco_p`,
`reco_theta`, `reco_pT`, `reco_pz`, `reco_E_mu`, `reco_mu_px/py/pz`, `reco_proton_p/E/T`,
`reco_proton_px/py/pz`, `reco_proton_theta_beam`) and each grid's reco observable is evaluated with the
same `ndp.channels.reco_observables` functions the data goes through. What stays from the grid's binned
response: the efficiency of the true cell (`eff`), the MC background and the MC-stat variance; only the
migration matrix P is replaced.

Frame: the model works in the frame it was trained in (`spec.frame`; detector for the x60_het port).
Truth momenta are rotated from the table's native frame into it; the smeared moments are stored as
detector-frame component columns (what the cache holds, what the TKI observables rotate themselves) and
as beam-frame angle columns (`reco_theta`, `reco_proton_theta_beam`, like the tuple's angle branches).
Units: tables in GeV, the model in MeV.

Conditioning on the reconstruction windows: the efficiency already contains the loss from the reco
kinematic cuts (muon window, proton window of `selection.params`), so the smeared reco is conditioned on
passing them — per event the weight w·eff is shared among the samples that pass
(`truncate_to_reco_windows`, default on). An event none of whose K samples pass contributes nothing and
its weight is reported in the fold info. Without truncation every sample carries w·eff/K and the
out-of-window samples are lost to the grid.

torch is only touched through ndp.surrogate.vbll_model (lazily), so the default environment imports this module.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from .base import Surrogate
from .binned import BinnedResponse
from ..channels.binning import Binning
from ..channels.observables import native_frame, frame_rotation_angle, rotate_about_x
from ..channels.signal import leading_proton
from ..config import REPO_ROOT
from ..events import TruthTable, M_MU, M_P
from ..io import timestamp

GEV = 1e-3          # MeV -> GeV
MEV = 1e3           # GeV -> MeV
RECO_COLUMNS = ("reco_p", "reco_theta", "reco_pT", "reco_pz", "reco_E_mu", "reco_mu_px", "reco_mu_py", "reco_mu_pz",
                "reco_proton_p", "reco_proton_E", "reco_proton_T", "reco_proton_px", "reco_proton_py", "reco_proton_pz",
                "reco_proton_theta_beam")


# ---- truth -> model inputs --------------------------------------------------------------------------
def truth_4vectors(channel, t: TruthTable, mask: np.ndarray, frame: str) -> tuple[np.ndarray, np.ndarray]:
    """(muon, leading proton) true 4-vectors (E, px, py, pz) in MeV in `frame`, for the masked events.

    The muon is the primary lepton; the proton is the leading one inside the channel's signal window,
    chosen in the channel frame (the window is defined there) and then rotated into `frame`."""
    ang = frame_rotation_angle(native_frame(t), frame)
    px, py, pz = rotate_about_x(t["lep_px"][mask], t["lep_py"][mask], t["lep_pz"][mask], ang)
    mu = np.stack([t["lep_E"][mask], px, py, pz], axis=1) * MEV
    lp = leading_proton(t, channel.frame, (channel.signal or {}).get("proton"))
    ang_p = frame_rotation_angle(channel.frame, frame)
    ppx, ppy, ppz = rotate_about_x(lp["px"][mask], lp["py"][mask], lp["pz"][mask], ang_p)
    pr = np.stack([lp["E"][mask], ppx, ppy, ppz], axis=1) * MEV
    return mu.astype(np.float32), pr.astype(np.float32)


# ---- model outputs -> synthetic reco columns ----------------------------------------------------------
def synthetic_reco(mu_mev: np.ndarray, pr_mev: np.ndarray, model_frame: str) -> dict:
    """Cached-column dictionary (GeV, detector-frame components, beam-frame angles) from smeared 4-vectors
    whose momentum components are in `model_frame`. Energies are rebuilt from the smeared 3-momentum and
    the particle mass (the model's E component is not used); the muon angle is the 3D beam-frame angle as
    the tuple's `muon_thetaX/Y` give it; the proton angle is the beam-frame angle like `MasterAnaDev_proton_theta`."""
    mu = np.asarray(mu_mev, np.float64) * GEV; pr = np.asarray(pr_mev, np.float64) * GEV
    to_det = frame_rotation_angle(model_frame, "detector"); to_beam = frame_rotation_angle(model_frame, "beam")
    mx, my, mz = rotate_about_x(mu[:, 1], mu[:, 2], mu[:, 3], to_det)
    bx, by, bz = rotate_about_x(mu[:, 1], mu[:, 2], mu[:, 3], to_beam)
    p_mu = np.sqrt(mx * mx + my * my + mz * mz)
    with np.errstate(invalid="ignore", divide="ignore"):
        th_mu = np.arccos(np.clip(bz / np.where(p_mu > 0, p_mu, np.nan), -1.0, 1.0))
    px, py, pz = rotate_about_x(pr[:, 1], pr[:, 2], pr[:, 3], to_det)
    qx, qy, qz = rotate_about_x(pr[:, 1], pr[:, 2], pr[:, 3], to_beam)
    p_p = np.sqrt(px * px + py * py + pz * pz)
    with np.errstate(invalid="ignore", divide="ignore"):
        th_p = np.arccos(np.clip(qz / np.where(p_p > 0, p_p, np.nan), -1.0, 1.0))
    E_p = np.sqrt(p_p * p_p + M_P * M_P)
    return {"reco_p": p_mu, "reco_theta": th_mu, "reco_pT": p_mu * np.sin(th_mu), "reco_pz": p_mu * np.cos(th_mu),
            "reco_E_mu": np.sqrt(p_mu * p_mu + M_MU * M_MU), "reco_mu_px": mx, "reco_mu_py": my, "reco_mu_pz": mz,
            "reco_proton_p": p_p, "reco_proton_E": E_p, "reco_proton_T": E_p - M_P,
            "reco_proton_px": px, "reco_proton_py": py, "reco_proton_pz": pz, "reco_proton_theta_beam": th_p}


def reco_window_pass(channel, r: dict) -> np.ndarray:
    """The selection's kinematic windows (MuonWindow, ProtonWindow of `selection.params`) on synthetic reco
    columns — the same predicates as the adapter's cutflow; NaN fails."""
    params = (channel.selection or {}).get("params") or {}
    mu, pr = params.get("muon"), params.get("proton")
    ok = np.ones(len(r["reco_p"]), bool)
    with np.errstate(invalid="ignore"):
        if mu:
            ok &= (r["reco_theta"] < np.deg2rad(mu["theta_max_deg"])) & (r["reco_p"] > mu["p_min_gev"]) & (r["reco_p"] < mu["p_max_gev"])
        if pr:
            ok &= (r["reco_proton_theta_beam"] < np.deg2rad(pr["theta_max_deg"])) & (r["reco_proton_p"] > pr["p_min_gev"]) & (r["reco_proton_p"] < pr["p_max_gev"])
    return ok


# ---- the surrogate ----------------------------------------------------------------------------------
class VBLLEventSurrogate(Surrogate):
    kind = "vbll_event"

    def __init__(self, binned: BinnedResponse, model, measurement, channel, n_samples: int = 20, seed: int = 0,
                 truncate_to_reco_windows: bool = True, meta: dict | None = None):
        if binned.binning != measurement.binning:
            raise ValueError(f"binned response grid != measurement {measurement.name} grid")
        super().__init__(binned.binning, meta)
        self.binned, self.model, self.measurement, self.channel = binned, model, measurement, channel
        self.n_samples, self.seed, self.truncate = int(n_samples), int(seed), bool(truncate_to_reco_windows)
        self.last_fold_info: dict = {}
        self.meta.setdefault("kind_note", "migration = VBLL event-level smearing; efficiency, background, MC-stat variance = the grid's binned response")
        self.meta.update({"measurement": measurement.name, "channel": channel.name, "n_samples": self.n_samples, "seed": self.seed,
                          "truncate_to_reco_windows": self.truncate, "model_frame": model.spec.frame, "model_head": model.spec.head_type,
                          "model_units": model.spec.units})

    # ---- construction / persistence -------------------------------------------------------------
    @classmethod
    def from_parts(cls, cfg, channel, measurement, binned: BinnedResponse, binned_path: str | Path | None,
                   model_name: str = "x60_het", n_samples: int = 20, seed: int = 0, truncate_to_reco_windows: bool = True,
                   model_dir: str | Path | None = None) -> "VBLLEventSurrogate":
        from .vbll_model import load_vbll_model
        mdir = Path(model_dir) if model_dir else cfg.surrogates / channel.name / "_vbll" / model_name
        model = load_vbll_model(mdir)
        meta = {"model_dir": _rel(mdir), "binned_dir": _rel(binned_path) if binned_path else None,
                "model_provenance": {k: (model.spec.meta.get(k)) for k in ("source", "checkpoint") if k in model.spec.meta},
                "binned_training_counts": binned.meta.get("training_counts"), "built": timestamp()}
        return cls(binned, model, measurement, channel, n_samples=n_samples, seed=seed,
                   truncate_to_reco_windows=truncate_to_reco_windows, meta=meta)

    def _arrays(self) -> dict:
        return {"eff": self.binned.eff}

    @classmethod
    def _from_arrays(cls, binning: Binning, arrays: dict, meta: dict) -> "VBLLEventSurrogate":
        from .base import load_surrogate
        from .vbll_model import load_vbll_model
        from ..channels import load_channel, load_measurement
        channel = load_channel(meta["channel"])
        measurement = load_measurement(channel, meta["measurement"])
        binned = load_surrogate(_abs(meta["binned_dir"]))
        model = load_vbll_model(_abs(meta["model_dir"]))
        return cls(binned, model, measurement, channel, n_samples=meta.get("n_samples", 20), seed=meta.get("seed", 0),
                   truncate_to_reco_windows=meta.get("truncate_to_reco_windows", True), meta=meta)

    # ---- delegation to the binned response ----------------------------------------------------------
    @property
    def eff(self) -> np.ndarray:
        return self.binned.eff

    def fold(self, true_cells: np.ndarray) -> np.ndarray:
        raise TypeError("vbll_event smears truth events (fold_table), not true-cell populations; "
                        "the grid's binned fold is available as .binned.fold")

    def fold_eff_only(self, true_cells: np.ndarray) -> np.ndarray:
        return self.binned.fold_eff_only(true_cells)

    def background(self, pot_data: float) -> np.ndarray:
        return self.binned.background(pot_data)

    def background_by_category(self, pot_data: float):
        return self.binned.background_by_category(pot_data)

    def fold_variance(self, true_cells: np.ndarray) -> np.ndarray:
        return self.binned.fold_variance(true_cells)

    def diagnostics(self) -> dict:
        d = {"kind": self.kind, "model_dir": self.meta.get("model_dir"), "binned_dir": self.meta.get("binned_dir"),
             "n_samples": self.n_samples, "seed": self.seed, "truncate_to_reco_windows": self.truncate,
             "model_frame": self.model.spec.frame, "closure": self.meta.get("closure")}
        d.update({f"binned_{k}": v for k, v in self.binned.diagnostics().items() if k in ("eff_population_weighted_mean", "n_den", "n_num")})
        return d

    # ---- the fold ----------------------------------------------------------------------------------------
    def fold_table(self, channel, t: TruthTable, weights=None, seed: int | None = None, mask=None) -> np.ndarray:
        """Reco cells of the signal (in phase space) events of `t`: per event w·eff(true cell), spread over the
        reco cells of its K smeared copies. `weights`: one per event of `t` (default t['weight'])."""
        cells, info = fold_tables([self], t, weights=weights, seed=self.seed if seed is None else seed, mask=mask)
        self.last_fold_info = info | {"per_measurement": info["per_measurement"].get(self.measurement.name)}
        return cells[self.measurement.name]


def _rel(p) -> str | None:
    if p is None:
        return None
    p = Path(p).resolve()
    try:
        return str(p.relative_to(REPO_ROOT))
    except ValueError:
        return str(p)


def _abs(p) -> Path:
    p = Path(p)
    return p if p.is_absolute() else REPO_ROOT / p


def fold_tables(surrogates: list, t: TruthTable, weights=None, seed: int = 0, mask=None, block: int = 200_000) -> tuple[dict, dict]:
    """Fold one truth table through several vbll_event surrogates (same model, channel, K, truncation) with a
    single smearing. Returns ({measurement: reco cells}, info)."""
    s0 = surrogates[0]
    channel, model, K, truncate = s0.channel, s0.model, s0.n_samples, s0.truncate
    for s in surrogates[1:]:
        if s.model is not model or s.channel.name != channel.name or s.n_samples != K or s.truncate != truncate:
            raise ValueError("fold_tables: every surrogate must share the model, channel, n_samples and truncation")
    if mask is None:
        mask = channel.in_phase_space(t) & channel.is_signal(t)
    mask = np.asarray(mask, bool)
    w_all = np.asarray(t["weight"] if weights is None else weights, float)
    if w_all.shape != (t.n,):
        raise ValueError("weights must be one per event of the table")
    w = w_all[mask]
    n = int(mask.sum())
    frame = model.spec.frame
    if frame is None:
        raise ValueError("the VBLL model's frame is not pinned (spec.frame is null)")
    mu, pr = truth_4vectors(channel, t, mask, frame)
    per_m = []
    for s in surrogates:
        x, y = s.measurement.truth_observables(channel, t)
        g = s.binning.digitize(x[mask], y[mask])
        inside = g >= 0
        eff_w = np.where(inside, s.binned.eff[np.where(inside, g, 0)], 0.0) * w
        per_m.append((s, eff_w, inside))
    cells = {s.measurement.name: np.zeros(s.binning.n_cells) for s in surrogates}
    info = {"n_events": n, "sum_w": float(w.sum()), "n_samples": K, "truncate_to_reco_windows": truncate, "model_frame": frame,
            "per_measurement": {s.measurement.name: {"n_true_out_of_grid": int((~ins).sum()), "sum_w_eff": float(ew.sum()),
                                                     "sum_w_eff_lost_all_samples_out_of_window": 0.0} for s, ew, ins in per_m},
            "n_events_all_samples_out_of_window": 0, "sample_pass_fraction": None}
    t0 = time.time(); n_pass_total = 0
    params = channel.observable_params
    for b, lo in enumerate(range(0, n, block)):
        sl = slice(lo, min(lo + block, n)); nb = sl.stop - sl.start
        mu_s = model.smear("muon", mu[sl], n_samples=K, seed=seed * 1_000_003 + 2 * b)
        pr_s = model.smear("proton", pr[sl], n_samples=K, seed=seed * 1_000_003 + 2 * b + 1)
        recos = [synthetic_reco(mu_s[k], pr_s[k], frame) for k in range(K)]
        if truncate:
            pass_k = np.stack([reco_window_pass(channel, r) for r in recos])
            n_pass = pass_k.sum(axis=0)
            share = np.where(n_pass > 0, 1.0 / np.maximum(n_pass, 1), 0.0)
            info["n_events_all_samples_out_of_window"] += int((n_pass == 0).sum())
            for s, eff_w, _ in per_m:
                info["per_measurement"][s.measurement.name]["sum_w_eff_lost_all_samples_out_of_window"] += float(eff_w[sl][n_pass == 0].sum())
        else:
            pass_k = np.ones((K, nb), bool); share = np.full(nb, 1.0 / K)
        n_pass_total += int(pass_k.sum())
        for k, r in enumerate(recos):
            use = pass_k[k]
            if not use.any():
                continue
            for s, eff_w, _ in per_m:
                xr, yr = s.measurement.reco_observables(r, params=params)
                ww = eff_w[sl] * share
                h, _, _ = s.binning.histogram(xr[use], yr[use], ww[use])
                cells[s.measurement.name] += h
    info["sample_pass_fraction"] = float(n_pass_total / max(K * n, 1))
    info["elapsed_s"] = round(time.time() - t0, 1)
    for s in surrogates:
        info["per_measurement"][s.measurement.name]["sum_folded"] = float(cells[s.measurement.name].sum())
    return cells, info


# ---- building the per-measurement wrappers + their closure on the training MC -------------------------------
def _iter_truth(cfg, channel):
    """The official MC truth, one playlist skim (or legacy truth cache) at a time."""
    from ..products import has_products, playlists, playlist_dir, legacy_files
    from ..adapters.minerva_anatuple import cache_tag
    if has_products(channel):
        for beam, pl in playlists(channel, "mc"):
            p = playlist_dir(cfg, channel, beam, pl) / f"truth_{pl}_skim.npz"
            if not p.exists():
                raise FileNotFoundError(f"missing playlist product {p}")
            yield f"{beam}/{pl}", TruthTable.load(p)
        return
    cache = cfg.require("data_dir") / "cache"
    for fn in legacy_files(channel, "mc"):
        yield cache_tag(fn), TruthTable.load(cache / f"truth_{cache_tag(fn)}.npz")


def closure_target(binned: BinnedResponse) -> np.ndarray | None:
    """Selected signal with true cell inside the grid, per reco cell = the binned fold of the training truth."""
    return None if binned.migration_counts is None else binned.migration_counts.sum(axis=1)


def build_vbll(channel, measurements: list, cfg, model_name: str = "x60_het", n_samples: int = 20, seed: int = 0,
               truncate_to_reco_windows: bool = True, out_root: Path | None = None, closure: bool = True, log=print) -> list[dict]:
    """Wrap each measurement's binned response with the VBLL model, save `surrogates/<channel>/<measurement>/vbll_<model>/`,
    and (closure=True) fold the official MC truth signal once for all of them and compare with each grid's selected
    signal per reco cell (the binned fold of the same truth, exact by construction)."""
    from ..pipeline import resolve_surrogate
    wrappers, paths = [], []
    for m in measurements:
        sur, sp = resolve_surrogate(channel, m, cfg)
        if sur is None or not isinstance(sur, BinnedResponse):
            raise FileNotFoundError(f"measurement {m.name}: no binned response (ndp surrogate build --kind binned first)")
        w = VBLLEventSurrogate.from_parts(cfg, channel, m, sur, sp, model_name=model_name, n_samples=n_samples, seed=seed,
                                          truncate_to_reco_windows=truncate_to_reco_windows)
        if len(wrappers) and w.model is not wrappers[0].model:      # one model instance for all (fold_tables requires identity)
            w.model = wrappers[0].model
        wrappers.append(w); paths.append((Path(out_root) if out_root else m.surrogate_root(cfg)) / f"vbll_{model_name}")
    results = []
    if closure:
        cells = {w.measurement.name: np.zeros(w.binning.n_cells) for w in wrappers}
        infos = []
        for label, t in _iter_truth(cfg, channel):
            c, info = fold_tables(wrappers, t, seed=seed)
            for k in cells:
                cells[k] += c[k]
            infos.append({"input": label, **{k: v for k, v in info.items() if k != "per_measurement"}})
            log(f"  folded {label}: {info['n_events']} signal events in {info['elapsed_s']} s (sample pass fraction {info['sample_pass_fraction']:.3f})")
            del t
    for w, p in zip(wrappers, paths):
        res = {"measurement": w.measurement.name, "path": str(p)}
        if closure:
            target = closure_target(w.binned)
            pred = cells[w.measurement.name]
            dev = pred - target
            big = target > 50
            clo = {"n_pred": float(pred.sum()), "n_reco": float(target.sum()), "ratio_pred_over_reco": float(pred.sum() / target.sum()) if target.sum() else None,
                   "max_abs_dev": float(np.abs(dev).max()), "max_rel_dev_cells_gt_50": float(np.max(np.abs(dev[big]) / target[big])) if big.any() else None,
                   "pred_cells": pred.tolist(), "reco_cells": target.tolist(), "binned_fold_cells": w.binned.fold(w.binned.den_counts).tolist() if w.binned.den_counts is not None else None,
                   "compared_with": "selected signal in the true phase space with true cell inside the grid (migration column sums)",
                   "inputs": infos, "built": timestamp()}
            w.meta["closure"] = {k: v for k, v in clo.items() if k not in ("pred_cells", "reco_cells", "binned_fold_cells", "inputs")}
            res["closure"] = clo
        w.save(p)
        if closure:
            (p / "closure.json").write_text(json.dumps(res["closure"], indent=1))
            log(f"{w.measurement.name}: saved {p}; closure pred/reco = {res['closure']['ratio_pred_over_reco']:.4f}, "
                f"max rel dev (cells > 50) = {res['closure']['max_rel_dev_cells_gt_50']}")
        else:
            log(f"{w.measurement.name}: saved {p}")
        results.append(res)
    return results
