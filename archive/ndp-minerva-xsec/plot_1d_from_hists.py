"""1D p_T and p_parallel distributions + migration matrices from hists.npz.

Post-processing of a plot_2d_ptpl.py run — no streaming. Definitional note:
all 1D quantities are projections of the in-grid 2D objects, i.e. they are
defined WITHIN the 2D phase space (1.5 <= p_par < 60, 0 <= p_T < 4.5 GeV/c)
exactly like the paper's 1D results (projections of the 2D grid). Events
outside the grid in either variable are excluded (0.75% of selected signal,
see the parent run's summary.json).

Outputs: pt_1d.png / pl_1d.png (data vs POT-scaled MC, with ratio panel),
migration_pt.png (14x14), migration_pl.png (16x16), all row-normalized
P(reco|true) for selected signal; plus a small JSON with the projections.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

from runlog_tools import (RunLog, add_label, args_to_inputs, default_outdir,
                          make_parser)
from xsec import binning


def project_migration(m224):
    """(224,224) gid-indexed migration -> (14,14) pT and (16,16) p|| blocks."""
    m4 = m224.reshape(binning.N_PT, binning.N_PL, binning.N_PT, binning.N_PL)
    return m4.sum(axis=(1, 3)), m4.sum(axis=(0, 2))


def overlay_1d(edges, data, mc_scaled, xlabel, logx, path, title):
    centers = 0.5 * (edges[:-1] + edges[1:])
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(7.5, 6.5), sharex=True,
                                  height_ratios=[3, 1])
    ax.stairs(mc_scaled, edges, color="C1", lw=1.5,
              label="MC selected x POT ratio (unweighted CV)")
    ax.errorbar(centers, data, yerr=np.sqrt(data), fmt="ko", ms=3.5,
                label="data")
    ax.set_ylabel("selected events / bin")
    ax.set_yscale("log")
    ax.set_title(title)
    ax.legend()
    ratio = np.divide(data, mc_scaled, out=np.zeros_like(mc_scaled),
                      where=mc_scaled > 0)
    rerr = np.divide(np.sqrt(data), mc_scaled, out=np.zeros_like(mc_scaled),
                     where=mc_scaled > 0)
    axr.errorbar(centers, ratio, yerr=rerr, fmt="ko", ms=3.5)
    axr.axhline(1.0, color="C1", lw=1)
    axr.set_ylabel("data / MC")
    axr.set_xlabel(xlabel)
    axr.set_ylim(0.7, 1.3)
    if logx:
        ax.set_xscale("log")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return ratio


def plot_migration_1d(m, edges, axis_label, path, title, logscale=False):
    """Migration on PHYSICAL axes: cell sizes show the real bin widths."""
    rows = m.sum(axis=1, keepdims=True)
    norm = np.divide(m, rows, out=np.zeros_like(m), where=rows > 0)
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    # norm[true, reco]: rows are y (true), columns x (reco)
    pcm = ax.pcolormesh(edges, edges, np.ma.masked_equal(norm, 0.0),
                        norm=LogNorm(vmin=1e-4, vmax=1.0), cmap="viridis")
    if logscale:
        ax.set_xscale("log")
        ax.set_yscale("log")
    ax.set_xticks(edges)
    ax.set_yticks(edges)
    ax.set_xticklabels([f"{e:g}" for e in edges], rotation=90, fontsize=7)
    ax.set_yticklabels([f"{e:g}" for e in edges], fontsize=7)
    ax.minorticks_off()
    for e in edges:
        ax.axhline(e, color="w", lw=0.2, alpha=0.4)
        ax.axvline(e, color="w", lw=0.2, alpha=0.4)
    ax.set_xlabel(f"reco {axis_label} [GeV/c]")
    ax.set_ylabel(f"true {axis_label} [GeV/c]")
    ax.set_title(title)
    fig.colorbar(pcm, ax=ax, label="row-normalized P(reco | true)")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return norm


def main():
    parser = make_parser("1D pT/p|| distributions + migrations from a "
                         "plot_2d_ptpl.py hists.npz (no streaming).")
    parser.add_argument("--hists", required=True,
                        help="hists.npz from a plot_2d_ptpl.py run")
    parser.add_argument("--outdir", default=None)
    add_label(parser)
    args = parser.parse_args()

    outdir = Path(args.outdir) if args.outdir else default_outdir(__file__)
    outdir.mkdir(parents=True, exist_ok=True)

    with RunLog(__file__, "1D projections + migrations from 2D hists",
                inputs=args_to_inputs(args)) as log:
        f = np.load(args.hists)
        pot_scale = float(f["pot_data"]) / float(f["pot_mc"])
        data2d, mc2d = f["data2d"], f["mc2d_selected"]
        mig_pt, mig_pl = project_migration(f["migration"])

        data_pt, data_pl = data2d.sum(axis=1), data2d.sum(axis=0)
        mc_pt, mc_pl = mc2d.sum(axis=1) * pot_scale, mc2d.sum(axis=0) * pot_scale

        r_pt = overlay_1d(binning.PT_EDGES_GEV, data_pt, mc_pt,
                          r"$p_{T,\mu}$ [GeV/c]", False, outdir / "pt_1d.png",
                          "Playlist 1A, selected events (grid-projected)")
        r_pl = overlay_1d(binning.PL_EDGES_GEV, data_pl, mc_pl,
                          r"$p_{\parallel,\mu}$ [GeV/c]", True, outdir / "pl_1d.png",
                          "Playlist 1A, selected events (grid-projected)")
        plot_migration_1d(mig_pt, binning.PT_EDGES_GEV, r"$p_T$",
                          outdir / "migration_pt.png",
                          r"$p_T$ migration — selected signal")
        plot_migration_1d(mig_pl, binning.PL_EDGES_GEV, r"$p_\parallel$",
                          outdir / "migration_pl.png",
                          r"$p_\parallel$ migration — selected signal",
                          logscale=True)

        diag_pt = np.diag(mig_pt / np.maximum(mig_pt.sum(axis=1, keepdims=True), 1))
        diag_pl = np.diag(mig_pl / np.maximum(mig_pl.sum(axis=1, keepdims=True), 1))
        summary = {
            "pot_scale_data_over_mc": pot_scale,
            "data_total_in_grid": float(data2d.sum()),
            "mc_scaled_total_in_grid": float(mc2d.sum() * pot_scale),
            "data_over_mc_overall": float(data2d.sum() / (mc2d.sum() * pot_scale)),
            "ratio_pt_bins": [round(float(x), 4) for x in r_pt],
            "ratio_pl_bins": [round(float(x), 4) for x in r_pl],
            "migration_diagonal_pt": [round(float(x), 4) for x in diag_pt],
            "migration_diagonal_pl": [round(float(x), 4) for x in diag_pl],
        }
        (outdir / "summary_1d.json").write_text(json.dumps(summary, indent=2))

        log.out("outdir", str(outdir))
        log.out("plots", [str(outdir / p) for p in
                          ("pt_1d.png", "pl_1d.png", "migration_pt.png",
                           "migration_pl.png")])
        log.out("data_over_mc_overall", summary["data_over_mc_overall"])
        log.out("summary_1d", str(outdir / "summary_1d.json"))
        print(json.dumps({k: v for k, v in summary.items()
                          if not k.startswith(("ratio", "migration"))}, indent=2))


if __name__ == "__main__":
    main()
