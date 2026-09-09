"""Learn a measurement's detector surrogates from the experiment's paired truth/reco MC.

The surrogate answers "what does the detector do to signal in this measurement's true
cells"; it is learned once per (channel, measurement) from the cached MC tables and
certified by closure: folding the training MC's own truth must reproduce its selected
reco counts cell by cell (exact for the binned response).

    den : signal & in true phase space (Truth tree)              -> true cells
    num : reco-selected signal & in true phase space              -> true cell, reco cell
    bkg : every other reco-selected candidate in the reco grid    -> reco cell (per POT)
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ..events import TruthTable
from ..io import cheap_fingerprint, timestamp
from .binned import BinnedResponse
from .parametric import SmearingSurrogate


def load_training_cache(cfg, channel) -> dict:
    """The cached MC tables of the channel's first reco_mc_file (built by `ndp data cache`)."""
    from ..adapters.minerva_anatuple import cache_tag, load_reco_cache
    data_dir = cfg.require("data_dir")
    cache = data_dir / "cache"
    mc_file = channel.data["reco_mc_files"][0]
    tag = cache_tag(mc_file)
    truth = TruthTable.load(cache / f"truth_{tag}.npz")
    reco = load_reco_cache(cache / f"reco_{tag}.npz")
    reco_truth = TruthTable.load(cache / f"reco_{tag}_truthcols.npz")
    return {"truth": truth, "reco": reco, "reco_truth": reco_truth, "tag": tag, "mc_file": data_dir / mc_file,
            "pot_mc": truth.norm.pot}


def training_arrays(channel, measurement, truth: TruthTable, reco: dict, reco_truth: TruthTable) -> dict:
    """The fit() keyword arrays for both surrogate kinds, plus the background sample."""
    sig_den = channel.is_signal(truth) & channel.in_phase_space(truth)
    xd, yd = measurement.truth_observables(channel, truth)
    passed = np.asarray(reco["passed"], bool)
    sig_rt = channel.is_signal(reco_truth) & channel.in_phase_space(reco_truth)
    sig_num = passed & sig_rt
    xn, yn = measurement.truth_observables(channel, reco_truth)
    xr, yr = measurement.reco_observables(reco)
    bkg = passed & ~sig_rt                                       # non-signal + out-of-phase-space signal
    return {"kw": dict(x_true_den=xd[sig_den], y_true_den=yd[sig_den], x_true_num=xn[sig_num], y_true_num=yn[sig_num],
                       x_reco_num=xr[sig_num], y_reco_num=yr[sig_num]),
            "bkg": dict(x_reco_bkg=xr[bkg], y_reco_bkg=yr[bkg]),
            "counts": {"n_den": int(sig_den.sum()), "n_num": int(sig_num.sum()), "n_bkg": int(bkg.sum()), "n_passed": int(passed.sum())},
            "reco_selected": (xr[passed], yr[passed])}


def build_surrogates(channel, measurement, cfg, kinds=("binned", "parametric"), n_samples: int = 20,
                     out_root: Path | None = None, log=print) -> list[dict]:
    """Fit, save (with a response figure) and closure-check the requested surrogate kinds."""
    from ..compare import plots
    tc = load_training_cache(cfg, channel)
    arrs = training_arrays(channel, measurement, tc["truth"], tc["reco"], tc["reco_truth"])
    b = measurement.binning
    root = Path(out_root) if out_root else measurement.surrogate_root(cfg)
    meta = {"channel": channel.name, "measurement": measurement.name, "measurement_spec": measurement.to_dict(),
            "built": timestamp(), "training_mc": cheap_fingerprint(tc["mc_file"]), "pot_mc": tc["pot_mc"],
            "generator": tc["truth"].meta.get("generator"), "selection": channel.selection.get("name"),
            "phase_space": channel.phase_space, "signal": channel.signal, "training_counts": arrs["counts"]}
    results = []
    for kind in kinds:
        if kind == "binned":
            sur = BinnedResponse.fit(b, **arrs["kw"], **arrs["bkg"], pot_mc=tc["pot_mc"], meta=meta)
            out = root / f"binned_{tc['tag']}"; title = f"binned response, MC {tc['tag']}, {measurement.name}"
        elif kind == "parametric":
            sur = SmearingSurrogate.fit(b, **arrs["kw"], mode="diff,ratio", n_samples=n_samples, meta=meta)
            out = root / f"parametric_{tc['tag']}"; title = f"parametric smearing, MC {tc['tag']}, {measurement.name}"
        else:
            raise ValueError(f"unknown surrogate kind {kind!r}")
        sur.save(out)
        plots.response_figure(sur, out / "response.png", title)
        clo = closure(channel, measurement, sur, tc["truth"], tc["reco"], tc["reco_truth"])
        res = {"kind": kind, "path": str(out), "diagnostics": sur.diagnostics(), "closure": clo}
        results.append(res)
        log(f"{kind}: saved {out}; closure max|pred-reco| = {clo['max_abs_dev']:.3g} on {clo['n_reco']:.0f} selected in-grid events")
    return results


def closure(channel, measurement, surrogate, truth: TruthTable, reco: dict, reco_truth: TruthTable | None = None) -> dict:
    """Fold the training truth and compare with the training reco counts (exact for binned)."""
    sumw, _, _, _ = measurement.truth_cells(channel, truth)
    pot = truth.norm.pot
    if surrogate.kind == "binned_response":
        pred = surrogate.fold(sumw) + surrogate.background(pot)
    else:  # event-level smearing: fold the truth events themselves
        mask = channel.in_phase_space(truth) & channel.is_signal(truth)
        x, y = measurement.truth_observables(channel, truth)
        pred = surrogate.fold_events(x[mask], y[mask], truth["weight"][mask], np.random.default_rng(0))
    passed = np.asarray(reco["passed"], bool)
    xr, yr = measurement.reco_observables(reco)
    target = "all selected candidates"
    if surrogate.kind != "binned_response" and reco_truth is not None:
        # a smearing model predicts signal only: compare with the selected signal (in phase space)
        passed = passed & channel.is_signal(reco_truth) & channel.in_phase_space(reco_truth)
        target = "selected signal in the true phase space"
    h, _, _ = measurement.binning.histogram(xr[passed], yr[passed])
    dev = pred - h
    return {"n_pred": float(pred.sum()), "n_reco": float(h.sum()), "max_abs_dev": float(np.abs(dev).max()),
            "exact": bool(np.allclose(pred, h, atol=1e-6)) if surrogate.kind == "binned_response" else None,
            "ratio_pred_over_reco": float(pred.sum() / h.sum()) if h.sum() else None, "compared_with": target,
            "max_rel_dev_cells_gt_50": float(np.max(np.abs(dev[h > 50]) / h[h > 50])) if (h > 50).any() else None}
