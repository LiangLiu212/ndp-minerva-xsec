"""Folded-space comparison: surrogate-smeared model vs reconstructed, selected data counts.

    N_true[j]   expected signal events in true cell j at the data exposure
    N_reco[i]   = sum_j P[i,j] eff[j] N_true[j] + background[i]
    data[i]     selected data candidates in reco cell i

Goodness of fit is Poisson: Baker–Cousins -2 ln lambda summed over cells, plus a Pearson
chi2 with the surrogate's MC-statistical variance added to the denominator. No unfolding,
no regularisation — the model is pushed through the detector, not the data pulled back.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ..channels import ChannelSpec
from ..events import TruthTable
from ..surrogate.base import Surrogate


def expected_true_cells(channel: ChannelSpec, t: TruthTable, pot_data: float, *, phi_per_pot=None, n_nucleons=None) -> dict:
    sumw, sumw2, n_out, mask = channel.truth_cells(t)
    norm = t.norm
    if norm.kind == "pot":
        scale = pot_data / float(norm.pot)
        how = f"N_true = N_mc * POT_data/POT_mc ({scale:.5g})"
    elif norm.kind == "xsec_per_nucleon":
        phi = phi_per_pot if phi_per_pot is not None else channel.normalization["phi_per_pot_cm2"]
        nn = n_nucleons if n_nucleons is not None else channel.normalization["n_nucleons"]
        scale = float(norm.xsec_per_unit_weight) * float(nn) * float(phi) * pot_data
        how = f"N_true = sigma_cell * N_nuc({nn:.3g}) * Phi({phi:.3g}) * POT({pot_data:.4g})"
    else:
        raise ValueError("shape-only sample cannot predict an event rate")
    return {"N_true": sumw * scale, "var": sumw2 * scale ** 2, "scale": scale, "how": how,
            "n_signal_in_ps": int(mask.sum()), "n_out_of_grid": n_out}


def data_reco_cells(channel: ChannelSpec, cfg, reco_cache: Path | None = None) -> dict:
    """Selected data candidates in reco cells (from the cached selection or the AnaTuple)."""
    from ..adapters.minerva_anatuple import read_reco, read_pot
    data_dir = cfg.require("data_dir")
    files = channel.data["reco_data_files"]
    cells = np.zeros(channel.binning.n_cells)
    n_sel = n_out = 0
    pot = 0.0
    for fn in files:
        path = data_dir / fn
        cache = (reco_cache or data_dir / "cache") / f"reco_{Path(fn).stem.replace('MasterAnaDev_data_AnaTuple_run000', 'data')}.npz"
        if cache.exists():
            z = np.load(cache)
            passed, pT, pz = z["passed"], z["reco_pT"], z["reco_pz"]
        else:
            r = read_reco(path, is_mc=False)
            passed, pT, pz = r["passed"], r["reco"]["pT"], r["reco"]["pz"]
        h, _, out = channel.binning.histogram(pT[passed], pz[passed])
        cells += h; n_sel += int(passed.sum()); n_out += out
        pot += read_pot(path)["pot_used"]
    return {"cells": cells, "n_selected": n_sel, "n_out_of_grid": n_out, "pot": pot, "files": files}


def poisson_gof(data: np.ndarray, pred: np.ndarray, var_mc: np.ndarray | None = None) -> dict:
    data = np.asarray(data, float); pred = np.asarray(pred, float)
    use = pred > 0
    d, m = data[use], pred[use]
    with np.errstate(divide="ignore", invalid="ignore"):
        ll = np.where(d > 0, d * np.log(d / m), 0.0)
    m2lnl = float(2.0 * np.sum(m - d + ll))
    denom = m + (var_mc[use] if var_mc is not None else 0.0)
    pearson = float(np.sum((d - m) ** 2 / denom))
    ndf = int(use.sum())
    dropped = float(data[~use].sum())
    return {"minus2lnL": m2lnl, "pearson_chi2": pearson, "ndf": ndf, "minus2lnL_per_ndf": m2lnl / ndf if ndf else None,
            "pearson_per_ndf": pearson / ndf if ndf else None, "n_cells_used": ndf,
            "data_in_cells_with_zero_prediction": dropped}


def compare_folded(channel: ChannelSpec, t: TruthTable, surrogate: Surrogate, data: dict, *,
                   phi_per_pot=None, n_nucleons=None, use_events: bool = False, rng=None) -> dict:
    exp = expected_true_cells(channel, t, data["pot"], phi_per_pot=phi_per_pot, n_nucleons=n_nucleons)
    if use_events and hasattr(surrogate, "sample_reco"):
        mask = channel.in_phase_space(t) & channel.is_signal(t)
        x, y = channel.observables(t)
        pred_sig = surrogate.fold_events(x[mask], y[mask], t["weight"][mask] * exp["scale"], rng)
    else:
        pred_sig = surrogate.fold(exp["N_true"])
    bkg = surrogate.background(data["pot"])
    pred = pred_sig + bkg
    var_mc = surrogate.fold_variance(exp["N_true"])
    gof = poisson_gof(data["cells"], pred, var_mc)
    b = channel.binning
    return {
        "pred_cells": pred, "pred_signal_cells": pred_sig, "bkg_cells": bkg, "data_cells": data["cells"],
        "var_mc_cells": var_mc, "N_true_cells": exp["N_true"], "expected": {k: v for k, v in exp.items() if k not in ("N_true", "var")},
        "totals": {"data": float(data["cells"].sum()), "pred": float(pred.sum()), "pred_signal": float(pred_sig.sum()),
                   "bkg": float(bkg.sum()), "ratio_data_over_pred": float(data["cells"].sum() / pred.sum()) if pred.sum() else None},
        "gof": gof,
        "projections": {
            "x": {"edges": list(b.x_edges), "data": b.project(data["cells"], "x", False).tolist(), "pred": b.project(pred, "x", False).tolist()},
            "y": {"edges": list(b.y_edges), "data": b.project(data["cells"], "y", False).tolist(), "pred": b.project(pred, "y", False).tolist()},
        },
        "pot_data": data["pot"], "n_data_selected": data["n_selected"],
    }
