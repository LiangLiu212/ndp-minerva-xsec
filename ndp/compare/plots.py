"""Figures for both comparison modes. Self-contained matplotlib; ratio panel on every comparison."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA_COLOR, MODEL_COLORS = "#08519c", ["#a63603", "#006d2c", "#54278f", "#a50f15", "#636363"]


def _ratio_axes(n_panels: int, figsize):
    fig, axes = plt.subplots(2, n_panels, figsize=figsize, gridspec_kw={"height_ratios": [3, 1]}, sharex="col", squeeze=False)
    return fig, axes


def _step(ax, edges, y, **kw):
    ax.step(edges, np.append(y, y[-1]), where="post", **kw)


def unfolded_projections(rel, vecs: dict, out: Path, title: str, x_label="muon p_T [GeV/c]", y_label="muon p_|| [GeV/c]") -> Path:
    """dsigma/dpT and dsigma/dp|| projections of published data (with sqrt(diag cov)) and models."""
    fig, axes = _ratio_axes(2, (13, 6.5))
    for k, (axis, lab, edges) in enumerate((("pt", x_label, rel.pt_edges), ("pz", y_label, rel.pz_edges))):
        data_vec = np.where(rel.mask, rel.data * rel.areas, 0.0)                      # back to sigma per cell
        cov_cell = rel.cov_total * np.outer(rel.areas, rel.areas)
        d1, c1 = rel.project_1d(data_vec, cov_cell, axis)
        err = np.sqrt(np.maximum(np.diag(c1), 0))
        centres = 0.5 * (edges[:-1] + edges[1:]); widths = np.diff(edges)
        ax, axr = axes[0, k], axes[1, k]
        ax.errorbar(centres, d1, xerr=widths / 2, yerr=err, fmt="o", ms=4, color=DATA_COLOR, label=f"data (arXiv:{rel.arxiv})", capsize=2)
        for i, (name, vec) in enumerate(vecs.items()):
            m1, _ = rel.project_1d(np.where(rel.mask, vec * rel.areas, 0.0), None, axis)
            c = MODEL_COLORS[i % len(MODEL_COLORS)]
            _step(ax, edges, m1, color=c, lw=1.6, label=name)
            with np.errstate(invalid="ignore", divide="ignore"):
                _step(axr, edges, np.where(d1 > 0, m1 / d1, np.nan), color=c, lw=1.4)
        axr.fill_between(edges, 1 - np.append(err / np.where(d1 > 0, d1, np.nan), np.nan), 1 + np.append(err / np.where(d1 > 0, d1, np.nan), np.nan),
                         step="post", color=DATA_COLOR, alpha=0.15, lw=0)
        axr.axhline(1, color="k", lw=0.8, ls="--"); axr.set_ylim(0.5, 1.5)
        axr.set_xlabel(lab); ax.set_ylabel("d$\\sigma$/dx  [cm$^2$/(GeV/c)/nucleon]"); axr.set_ylabel("model / data")
        if axis == "pz":
            ax.set_xscale("log"); axr.set_xscale("log")
        ax.legend(fontsize=8)
    fig.suptitle(title); fig.tight_layout(); fig.savefig(out, dpi=120); plt.close(fig)
    return out


def folded_projections(res: dict, out: Path, title: str, x_label="reco muon p_T [GeV/c]", y_label="reco muon p_|| [GeV/c]") -> Path:
    fig, axes = _ratio_axes(2, (13, 6.5))
    for k, (axis, lab) in enumerate((("x", x_label), ("y", y_label))):
        p = res["projections"][axis]
        edges = np.asarray(p["edges"]); d = np.asarray(p["data"]); m = np.asarray(p["pred"])
        centres = 0.5 * (edges[:-1] + edges[1:]); widths = np.diff(edges)
        ax, axr = axes[0, k], axes[1, k]
        ax.errorbar(centres, d, xerr=widths / 2, yerr=np.sqrt(d), fmt="o", ms=4, color=DATA_COLOR, label=f"data ({res['n_data_selected']} selected)", capsize=2)
        _step(ax, edges, m, color=MODEL_COLORS[0], lw=1.6, label="model → surrogate (signal + bkg)")
        with np.errstate(invalid="ignore", divide="ignore"):
            _step(axr, edges, np.where(m > 0, d / m, np.nan), color=DATA_COLOR, lw=1.4)
            axr.fill_between(edges, 1 - np.append(np.sqrt(d) / np.where(m > 0, m, np.nan), np.nan),
                             1 + np.append(np.sqrt(d) / np.where(m > 0, m, np.nan), np.nan), step="post", color=DATA_COLOR, alpha=0.15, lw=0)
        axr.axhline(1, color="k", lw=0.8, ls="--"); axr.set_ylim(0.5, 1.5)
        axr.set_xlabel(lab); ax.set_ylabel("selected events / bin"); axr.set_ylabel("data / model")
        if axis == "y":
            ax.set_xscale("log"); axr.set_xscale("log")
        ax.legend(fontsize=8)
    g = res["gof"]
    fig.suptitle(f"{title}\n-2lnL/ndf = {g['minus2lnL']:.1f}/{g['ndf']}   Pearson χ²/ndf = {g['pearson_chi2']:.1f}/{g['ndf']}   "
                 f"data/pred = {res['totals']['ratio_data_over_pred']:.3f}", fontsize=10)
    fig.tight_layout(); fig.savefig(out, dpi=120); plt.close(fig)
    return out


def cell_ratio_map(binning, num: np.ndarray, den: np.ndarray, out: Path, title: str, vmin=0.5, vmax=1.5) -> Path:
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where(den > 0, num / den, np.nan)
    grid = binning.to_grid(r)   # [n_x, n_y]
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    im = ax.imshow(grid.T, origin="lower", aspect="auto", cmap="coolwarm", vmin=vmin, vmax=vmax)
    ax.set_xlabel(f"{binning.x_name} bin"); ax.set_ylabel(f"{binning.y_name} bin"); ax.set_title(title)
    fig.colorbar(im, ax=ax, label="ratio"); fig.tight_layout(); fig.savefig(out, dpi=120); plt.close(fig)
    return out


def response_figure(surrogate, out: Path, title: str) -> Path:
    b = surrogate.binning
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    eff = b.to_grid(surrogate.eff)
    im0 = axes[0].imshow(eff.T, origin="lower", aspect="auto", vmin=0, vmax=1)
    axes[0].set_title("efficiency (true cells)"); axes[0].set_xlabel(f"{b.x_name} bin"); axes[0].set_ylabel(f"{b.y_name} bin")
    fig.colorbar(im0, ax=axes[0])
    if hasattr(surrogate, "P"):
        with np.errstate(divide="ignore"):
            im1 = axes[1].imshow(np.log10(surrogate.P + 1e-4), origin="lower", aspect="auto")
        axes[1].set_title("log10 P(reco cell | true cell)"); axes[1].set_xlabel("true cell"); axes[1].set_ylabel("reco cell")
        fig.colorbar(im1, ax=axes[1])
    else:
        axes[1].plot(surrogate.ratio_std[:, 0], ".", label=f"{b.x_name} reco/true std"); axes[1].plot(surrogate.ratio_std[:, 1], ".", label=f"{b.y_name} reco/true std")
        axes[1].set_xlabel("true cell"); axes[1].set_ylabel("resolution"); axes[1].legend(); axes[1].set_title("smearing widths")
    fig.suptitle(title); fig.tight_layout(); fig.savefig(out, dpi=120); plt.close(fig)
    return out
