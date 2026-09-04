"""Before/after comparison of two plot_2d_ptpl.py runs (unweighted vs weighted).

Overlays data against POT-scaled MC for both weight modes and shows the
data/MC ratio panel — the direct picture of what the MnvTune v1 weights do.
Pure post-processing of two hists.npz; no streaming.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from runlog_tools import (RunLog, add_label, args_to_inputs, default_outdir,
                          make_parser)
from xsec import binning


def overlay(ax, axr, edges, data, mc_u, mc_w, logx):
    c = 0.5 * (edges[:-1] + edges[1:])
    ax.stairs(mc_u, edges, color="0.6", lw=1.5, label="MC unweighted CV")
    ax.stairs(mc_w, edges, color="C3", lw=1.7, label="MC MnvTune v1 (weighted)")
    ax.errorbar(c, data, yerr=np.sqrt(data), fmt="ko", ms=3.5, label="data")
    ax.set_yscale("log")
    ax.legend(fontsize=8)
    ru = np.divide(data, mc_u, out=np.zeros_like(mc_u), where=mc_u > 0)
    rw = np.divide(data, mc_w, out=np.zeros_like(mc_w), where=mc_w > 0)
    axr.plot(c, ru, "o", color="0.6", ms=4, label="data/MC unweighted")
    axr.plot(c, rw, "s", color="C3", ms=4, label="data/MC weighted")
    axr.axhline(1.0, color="k", lw=0.8)
    axr.axhline(1.118, color="C0", lw=0.8, ls="--", label="paper data/TuneV1 = 1.118")
    axr.set_ylim(0.7, 1.5)
    axr.legend(fontsize=7, ncol=2)
    if logx:
        ax.set_xscale("log")
        axr.set_xscale("log")


def main():
    parser = make_parser("Before/after (unweighted vs weighted) data-vs-MC comparison.")
    parser.add_argument("--unweighted", required=True, help="hists.npz, weight_mode=none")
    parser.add_argument("--weighted", required=True, help="hists.npz, weight_mode=cv")
    parser.add_argument("--outdir", default=None)
    add_label(parser)
    args = parser.parse_args()

    outdir = Path(args.outdir) if args.outdir else default_outdir(__file__)
    outdir.mkdir(parents=True, exist_ok=True)

    with RunLog(__file__, "weight-mode comparison", inputs=args_to_inputs(args)) as log:
        u = np.load(args.unweighted)
        w = np.load(args.weighted)
        psu = float(u["pot_data"] / u["pot_mc"])
        psw = float(w["pot_data"] / w["pot_mc"])
        data = u["data2d"]
        mc_u = u["mc2d_selected"] * psu
        mc_w = w["mc2d_selected"] * psw

        fig, axes = plt.subplots(2, 2, figsize=(13, 7), sharex="col",
                                 height_ratios=[3, 1])
        overlay(axes[0, 0], axes[1, 0], binning.PT_EDGES_GEV,
                data.sum(1), mc_u.sum(1), mc_w.sum(1), False)
        overlay(axes[0, 1], axes[1, 1], binning.PL_EDGES_GEV,
                data.sum(0), mc_u.sum(0), mc_w.sum(0), True)
        axes[1, 0].set_xlabel(r"$p_{T,\mu}$ [GeV/c]")
        axes[1, 1].set_xlabel(r"$p_{\parallel,\mu}$ [GeV/c]")
        axes[0, 0].set_ylabel("selected events / bin")
        axes[1, 0].set_ylabel("data / MC")
        fig.suptitle("Playlist 1A: unweighted CV vs MnvTune v1 — selected event rate")
        fig.tight_layout()
        fig.savefig(outdir / "weight_compare_1d.png", dpi=160)
        plt.close(fig)

        summary = {
            "data_over_mc_total_unweighted": float(data.sum() / mc_u.sum()),
            "data_over_mc_total_weighted": float(data.sum() / mc_w.sum()),
            "net_weight_on_selected_mc": float(w["mc2d_selected"].sum() / u["mc2d_selected"].sum()),
            "paper_data_over_tunev1": 1.118,
        }
        (outdir / "summary.json").write_text(json.dumps(summary, indent=2))
        log.out("outdir", str(outdir))
        log.out("plot", str(outdir / "weight_compare_1d.png"))
        log.out("summary", summary)
        print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
