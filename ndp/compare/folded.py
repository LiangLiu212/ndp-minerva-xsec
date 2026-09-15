"""Folded-space comparison: surrogate-smeared model vs reconstructed, selected data counts.

    N_true[j]   expected signal events in true cell j at the data exposure
    N_reco[i]   = sum_j P[i,j] eff[j] N_true[j] + background[i]
    data[i]     selected data candidates in reco cell i

The cells are the *measurement's*: any observable pair the analyst declared, at truth level
for N_true and at reco level for the data. Goodness of fit is Poisson: Baker–Cousins
-2 ln lambda summed over cells, plus a Pearson chi2 with the surrogate's MC-statistical
variance added to the denominator. No unfolding, no regularisation — the model is pushed
through the detector, not the data pulled back (forward folding).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ..channels import ChannelSpec, Measurement
from ..events import TruthTable
from ..surrogate.base import Surrogate


def expected_true_cells(channel: ChannelSpec, measurement: Measurement, t: TruthTable, pot_data: float, *,
                        phi_per_pot=None, n_nucleons=None) -> dict:
    sumw, sumw2, n_out, mask = measurement.truth_cells(channel, t)
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


def data_reco_cells(channel: ChannelSpec, measurement: Measurement, cfg, reco_cache: Path | None = None) -> dict:
    """Selected data candidates in the measurement's reco cells, from the playlist products or the
    per-file caches (ndp.products); the POT comes with the tables, never from the AnaTuples."""
    from ..products import load_reco, has_products, playlists
    from ..channels.selections import select
    r, pot, sources = load_reco(cfg, channel, "data")
    passed = select(channel, r)
    x, y = measurement.reco_observables(r, params=channel.observable_params)
    cells, _, n_out = measurement.binning.histogram(x[passed], y[passed])
    files = [f"{b}/{p}" for b, p in playlists(channel, "data")] if has_products(channel) else list(channel.data.get("reco_data_files", []))
    return {"cells": cells, "n_selected": int(passed.sum()), "n_out_of_grid": int(n_out), "pot": float(pot),
            "files": files, "sources": sources}


def data_reco_cells_multi(channel: ChannelSpec, measurements: list, cfg) -> dict:
    """`data_reco_cells` for several measurements in one pass over the data inputs (one playlist at a time)."""
    from ..channels.selections import select
    from ..products import iter_reco_chunks
    cells = {m.name: np.zeros(m.binning.n_cells) for m in measurements}
    n_sel = 0; n_out = {m.name: 0 for m in measurements}; pot = 0.0; sources = []
    for label, r, p, src in iter_reco_chunks(cfg, channel, "data"):
        sel = select(channel, r)
        n_sel += int(sel.sum()); pot += float(p); sources.append(src)
        for m in measurements:
            x, y = m.reco_observables(r, params=channel.observable_params)
            h, _, no = m.binning.histogram(x[sel], y[sel])
            cells[m.name] += h; n_out[m.name] += int(no)
        del r
    return {m.name: {"cells": cells[m.name], "n_selected": n_sel, "n_out_of_grid": n_out[m.name], "pot": pot,
                     "files": list(sources), "sources": list(sources)} for m in measurements}


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


def compare_folded(channel: ChannelSpec, measurement: Measurement, t: TruthTable, surrogate: Surrogate, data: dict, *,
                   phi_per_pot=None, n_nucleons=None, use_events: bool = False, rng=None, folding: str = "full",
                   truth_weights=None) -> dict:
    """folding = "full" (efficiency x migration), "eff_only" (efficiency, no migration: truth bins taken as reco bins)
    or "weighted" (the truth events carry their own per-event efficiency weights `truth_weights`, e.g. the factorised
    ansatz; histogrammed in truth bins and added to the background)."""
    if surrogate.binning != measurement.binning:
        raise ValueError(f"surrogate was built on a different grid ({surrogate.binning.x_name} x {surrogate.binning.y_name}, "
                         f"{surrogate.binning.n_cells} cells) than measurement {measurement.name} ({measurement.binning.n_cells} cells)")
    exp = expected_true_cells(channel, measurement, t, data["pot"], phi_per_pot=phi_per_pot, n_nucleons=n_nucleons)
    if folding == "eff_only":
        pred_sig = surrogate.fold_eff_only(exp["N_true"]); folding_how = "true cells x efficiency (no migration)"
    elif folding == "weighted":
        if truth_weights is None:
            raise ValueError("folding='weighted' needs truth_weights (one per truth event)")
        sumw, _, _, mask = measurement.truth_cells(channel, t, weights=np.asarray(truth_weights) * t["weight"])
        pred_sig = sumw * exp["scale"]; folding_how = "truth events x per-event efficiency weights (no migration)"
    elif use_events and hasattr(surrogate, "sample_reco"):
        mask = channel.in_phase_space(t) & channel.is_signal(t)
        x, y = measurement.truth_observables(channel, t)
        pred_sig = surrogate.fold_events(x[mask], y[mask], t["weight"][mask] * exp["scale"], rng)
        folding_how = "event-level smearing of the truth events"
    else:
        pred_sig = surrogate.fold(exp["N_true"])
        folding_how = "true cells x response"
    bkg = surrogate.background(data["pot"])
    pred = pred_sig + bkg
    var_mc = surrogate.fold_variance(exp["N_true"])
    gof = poisson_gof(data["cells"], pred, var_mc)
    b = measurement.binning
    out = {
        "pred_cells": pred, "pred_signal_cells": pred_sig, "bkg_cells": bkg, "data_cells": data["cells"],
        "var_mc_cells": var_mc, "N_true_cells": exp["N_true"], "expected": {k: v for k, v in exp.items() if k not in ("N_true", "var")},
        "folding": folding_how, "folding_mode": folding,
        "totals": {"data": float(data["cells"].sum()), "pred": float(pred.sum()), "pred_signal": float(pred_sig.sum()),
                   "bkg": float(bkg.sum()), "ratio_data_over_pred": float(data["cells"].sum() / pred.sum()) if pred.sum() else None},
        "gof": gof,
        "projections": {
            "x": {"edges": list(b.x_edges), "data": b.project(data["cells"], "x", False).tolist(), "pred": b.project(pred, "x", False).tolist(),
                  "bkg": b.project(bkg, "x", False).tolist(), "label": measurement.x.axis_label("reco"), "log": measurement.x.log},
            "y": {"edges": list(b.y_edges), "data": b.project(data["cells"], "y", False).tolist(), "pred": b.project(pred, "y", False).tolist(),
                  "bkg": b.project(bkg, "y", False).tolist(), "label": measurement.y.axis_label("reco"), "log": measurement.y.log},
        },
        "is_1d": measurement.is_1d,
        "pot_data": data["pot"], "n_data_selected": data["n_selected"], "n_data_out_of_grid": data["n_out_of_grid"],
    }
    return out
