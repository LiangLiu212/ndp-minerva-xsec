"""Unfolded-space comparison: model d2sigma per cell vs the published result + covariance.

The model's cross section is built from its truth events (signal, in phase space, binned
in the channel's true cells) with the sample's own normalisation, then handed to the
MINERvA benchmark engine's covariance-aware chi2 (total, shape-with-profiled-norm, norm
offset). Nothing is unfolded here — the *data* were unfolded by the experiment; the model
is compared in the space the experiment published.
"""
from __future__ import annotations

import numpy as np

from ..channels import ChannelSpec
from ..events import TruthTable
from .minerva_bridge import PaperRelease


def check_basis(channel: ChannelSpec, rel: PaperRelease) -> None:
    b = channel.binning
    if b.n_cells != rel.n_cells:
        raise ValueError(f"channel has {b.n_cells} cells, release {rel.n_cells}")
    if rel.pt_edges is not None and (not np.allclose(b.x_edges, rel.pt_edges) or not np.allclose(b.y_edges, rel.pz_edges)):
        raise ValueError("channel edges differ from the release edges")
    fam = {"ipt*n_pz + ipz": "ix*n_y + iy", "ipz*n_pt + ipt": "iy*n_x + ix"}
    if fam.get(rel.formula, rel.formula) != b.formula:
        raise ValueError(f"channel cell formula {b.formula!r} != release {rel.formula!r}")


def xsec_vector_from_truth(channel: ChannelSpec, t: TruthTable, *, phi_per_pot: float | None = None,
                           n_nucleons: float | None = None) -> dict:
    """d2sigma/(dx dy) per cell [cm^2/GeV^2/nucleon] with its MC-statistical variance."""
    sumw, sumw2, n_out, mask = channel.truth_cells(t)
    areas = channel.binning.areas()
    norm = t.norm
    if norm.kind == "xsec_per_nucleon":
        scale = float(norm.xsec_per_unit_weight)
        how = "sigma_cell = sum(w) * sigma_avg/sum(all w)"
    elif norm.kind == "pot":
        phi = phi_per_pot if phi_per_pot is not None else channel.normalization["phi_per_pot_cm2"]
        nn = n_nucleons if n_nucleons is not None else channel.normalization["n_nucleons"]
        scale = 1.0 / (float(norm.pot) * float(phi) * float(nn))
        how = f"sigma_cell = N_true / (POT_mc={norm.pot:.4g} * Phi={phi:.3g} * N_nuc={nn:.3g})"
    else:
        raise ValueError("a shape-only sample has no absolute cross section; supply a normalisation")
    sigma_cell = sumw * scale
    var_cell = sumw2 * scale ** 2
    return {"vec": sigma_cell / areas, "var": var_cell / areas ** 2, "sigma_cell": sigma_cell,
            "n_signal_in_ps": int(mask.sum()), "sumw_in_grid": float(sumw.sum()), "n_out_of_grid": n_out,
            "sigma_total_phase_space_cm2": float(sigma_cell.sum()), "normalisation": how}


def score_unfolded(channel: ChannelSpec, rel: PaperRelease, vec: np.ndarray, var: np.ndarray | None,
                   label: str) -> dict:
    """Score `vec` against the release. Two rows when the model has MC-stat variance: the
    paper covariance alone (what the shipped curves get) and paper + diag(model stat)."""
    check_basis(channel, rel)
    rows = []
    s = rel.compare(vec)
    rows.append({"label": label, "denominator": "paper_total", **s})
    if var is not None and np.any(var > 0):
        cov = rel.cov_total + np.diag(var)
        s2 = rel.compare(vec, cov=cov, mask=rel.mask & np.isfinite(vec))
        rows.append({"label": label, "denominator": "paper_total+model_stat", **s2})
    return {"rows": rows, "n_shared_cells": int(rel.mask.sum())}


def shipped_ranking(rel: PaperRelease) -> list[dict]:
    """Every generator curve the release ships, scored the same way (context for the user)."""
    out = []
    for name in rel.shipped_models():
        s = rel.compare(rel.shipped_curve(name))
        out.append({"label": name, "denominator": "paper_total", **s})
    out.sort(key=lambda r: r["chi2_total_per_ndf"])
    return out
