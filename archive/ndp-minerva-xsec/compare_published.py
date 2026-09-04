"""Compare an extract_xsec.py result against the published anc 2D cross section.

Loads the paper's data_result table (cm²/(GeV/c)², width-divided) and our
xsec.npz, and reports per-cell ratios over the 205 reported cells + the
integrated cross section. Plots: ratio map and per-cell overlay.
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
from extract_xsec import cell_areas


def load_anc(path):
    """anc data_result -> (dsigma[224], stat[224], total[224]) by GlobalID."""
    dsig = np.zeros(binning.N_CELLS)
    stat = np.zeros(binning.N_CELLS)
    tot = np.zeros(binning.N_CELLS)
    with open(path) as fh:
        fh.readline()
        for line in fh:
            p = line.replace(",", " ").split()
            if len(p) < 5:
                continue
            pl_bin, pt_bin = int(p[0]), int(p[1])
            gid = (pt_bin - 1) * binning.N_PL + (pl_bin - 1)
            dsig[gid], stat[gid], tot[gid] = float(p[2]), float(p[3]), float(p[4])
    return dsig, stat, tot


def main():
    parser = make_parser("Compare extracted d2sigma vs the published anc table.")
    parser.add_argument("--xsec", required=True, help="xsec.npz from extract_xsec.py")
    parser.add_argument("--published", default="config/published.json")
    parser.add_argument("--outdir", default=None)
    add_label(parser)
    args = parser.parse_args()

    outdir = Path(args.outdir) if args.outdir else default_outdir(__file__)
    outdir.mkdir(parents=True, exist_ok=True)

    pub = json.loads(Path(args.published).read_text())
    anc_path = Path(pub["anc_dir"]) / pub["key_files"]["data_result_2d"]
    anc_dsig, anc_stat, anc_tot = load_anc(anc_path)
    x = np.load(args.xsec)
    ours, ours_err = x["dsigma"], x["dsigma_err"]
    areas = cell_areas()

    with RunLog(__file__, "compare extracted vs published 2D xsec",
                inputs=args_to_inputs(args)) as log:
        rep = anc_dsig > 0                       # 205 reported cells
        ratio = np.divide(ours, anc_dsig, out=np.zeros_like(ours), where=rep)
        r = ratio[rep]
        pull = np.divide(ours - anc_dsig, np.hypot(ours_err, anc_tot),
                         out=np.zeros_like(ours), where=rep)[rep]
        integ_ours = float((ours * areas)[rep].sum())
        integ_anc = float((anc_dsig * areas)[rep].sum())

        summary = {
            "n_reported_cells": int(rep.sum()),
            "n_ours_positive_in_reported": int((ours[rep] > 0).sum()),
            "ratio_median": float(np.median(r[r > 0])),
            "ratio_p16_p84": [float(np.percentile(r[r > 0], 16)),
                              float(np.percentile(r[r > 0], 84))],
            "integrated_ours": integ_ours,
            "integrated_published": integ_anc,
            "integrated_ratio": integ_ours / integ_anc,
            "median_pull_vs_total_unc": float(np.median(np.abs(pull[pull != 0]))),
        }

        # ratio map (pt x pl)
        rmap = np.full((binning.N_PT, binning.N_PL), np.nan)
        for gid in np.where(rep)[0]:
            rmap[gid // binning.N_PL, gid % binning.N_PL] = ratio[gid]
        fig, ax = plt.subplots(figsize=(8, 5))
        pcm = ax.pcolormesh(binning.PL_EDGES_GEV, binning.PT_EDGES_GEV,
                            np.ma.masked_invalid(rmap), vmin=0.7, vmax=1.3,
                            cmap="RdBu_r")
        ax.set_xscale("log")
        ax.set_xlabel(r"$p_{\parallel,\mu}$ [GeV/c]"); ax.set_ylabel(r"$p_{T,\mu}$ [GeV/c]")
        ax.set_title(f"extracted / published (1A, stat-only) — median {summary['ratio_median']:.3f}")
        fig.colorbar(pcm, ax=ax, label="ours / published")
        fig.tight_layout(); fig.savefig(outdir / "ratio_map.png", dpi=160); plt.close(fig)

        # per-cell overlay (sorted by published value)
        order = np.argsort(anc_dsig[rep])[::-1]
        fig, ax = plt.subplots(figsize=(11, 4))
        idx = np.arange(rep.sum())
        ax.errorbar(idx, ours[rep][order], yerr=ours_err[rep][order], fmt="ko",
                    ms=2.5, label="extracted (1A, stat-only)")
        ax.errorbar(idx, anc_dsig[rep][order], yerr=anc_tot[rep][order], fmt="r.",
                    ms=3, alpha=0.6, label="published (total unc.)")
        ax.set_yscale("log"); ax.set_xlabel("reported cell (sorted by published d2sigma)")
        ax.set_ylabel(r"d$^2\sigma$/d$p_T$d$p_\parallel$ [cm$^2$/(GeV/c)$^2$]")
        ax.legend(); ax.set_title("playlist 1A vs published")
        fig.tight_layout(); fig.savefig(outdir / "overlay.png", dpi=160); plt.close(fig)

        (outdir / "summary.json").write_text(json.dumps(summary, indent=2))
        log.out("outdir", str(outdir))
        log.out("plots", [str(outdir / "ratio_map.png"), str(outdir / "overlay.png")])
        log.out("summary", summary)
        print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
