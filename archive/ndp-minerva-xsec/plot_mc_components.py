"""MC composition plots from a plot_2d_ptpl.py hists.npz: total / signal / background.

Total    = all reco-selected MC events;
Signal   = selected AND true nu_mu CC (is_signal);
Background = selected AND NOT is_signal (NC hadron fakes + wrong-sign/other CC).

Raw MC counts (no POT scaling — this shows the MC's own composition); the
background is ~0.24% of the total, hence log scales everywhere.

Outputs: mc_components_1d.png (pT and p|| grid projections, overlaid) and
mc_components_2d.png (three 2D maps on a common color scale).
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


def main():
    parser = make_parser("MC total/signal/background composition plots from hists.npz.")
    parser.add_argument("--hists", required=True)
    parser.add_argument("--outdir", default=None)
    add_label(parser)
    args = parser.parse_args()

    outdir = Path(args.outdir) if args.outdir else default_outdir(__file__)
    outdir.mkdir(parents=True, exist_ok=True)

    with RunLog(__file__, "MC components: total/signal/background",
                inputs=args_to_inputs(args)) as log:
        f = np.load(args.hists)
        comps = {"total (selected)": (f["mc2d_selected"], "red"),
                 "signal": (f["mc2d_signal"], "C0"),
                 "background": (f["mc2d_bkg"], "black")}
        pot_mc = float(f["pot_mc"])

        # ---- 1D projections -------------------------------------------------
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
        for (axis_label, edges, proj_axis, logx), ax in zip(
                [("$p_{T,\\mu}$ [GeV/c]", binning.PT_EDGES_GEV, 1, False),
                 ("$p_{\\parallel,\\mu}$ [GeV/c]", binning.PL_EDGES_GEV, 0, True)],
                axes):
            for name, (h, color) in comps.items():
                ax.stairs(h.sum(axis=proj_axis), edges, color=color, lw=1.6,
                          label=name)
            ax.set_yscale("log")
            if logx:
                ax.set_xscale("log")
            ax.set_xlabel(axis_label)
            ax.set_ylabel("selected MC events / bin")
            ax.legend(fontsize=9)
        fig.suptitle(f"MC composition, playlist 1A (grid-projected, "
                     f"unweighted CV, {pot_mc:.3e} POT)")
        fig.tight_layout()
        fig.savefig(outdir / "mc_components_1d.png", dpi=160)
        plt.close(fig)

        # ---- 2D maps on a common color scale --------------------------------
        vmax = comps["total (selected)"][0].max()
        fig, axes = plt.subplots(1, 3, figsize=(16, 4.4), sharey=True)
        for ax, (name, (h, _)) in zip(axes, comps.items()):
            pcm = ax.pcolormesh(binning.PL_EDGES_GEV, binning.PT_EDGES_GEV,
                                np.ma.masked_equal(h, 0.0),
                                norm=LogNorm(vmin=1.0, vmax=vmax), cmap="viridis")
            ax.set_xscale("log")
            ax.set_xlabel(r"$p_{\parallel,\mu}$ [GeV/c]")
            ax.set_title(f"{name}: {int(h.sum())} in grid")
        axes[0].set_ylabel(r"$p_{T,\mu}$ [GeV/c]")
        fig.colorbar(pcm, ax=axes, label="selected MC events / cell", pad=0.01)
        fig.suptitle("MC composition, playlist 1A (unweighted CV)")
        fig.savefig(outdir / "mc_components_2d.png", dpi=160)
        plt.close(fig)

        frac = float(f["mc2d_bkg"].sum() / f["mc2d_selected"].sum())
        log.out("outdir", str(outdir))
        log.out("plots", [str(outdir / "mc_components_1d.png"),
                          str(outdir / "mc_components_2d.png")])
        log.out("background_fraction_in_grid", frac)
        print(json.dumps({"background_fraction_in_grid": frac}, indent=2))


if __name__ == "__main__":
    main()
