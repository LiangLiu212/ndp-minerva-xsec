"""Paper-Fig.2-style event-rate panels from a plot_2d_ptpl.py hists.npz.

Replicates the presentation of arXiv:2106.16210 Fig. 2 (the
"multiplier-panel" event-rate comparison):

  - pt view: one panel per p_T bin, events vs p_parallel (log x);
  - pl view: one panel per p_parallel bin, events vs p_T (linear x);
  - y = events / (Delta p_T * Delta p_par), in units of 1e5 per (GeV/c)^2,
    each panel scaled by a "x m" multiplier so all panels share one y-range.

Differences vs the paper, by construction of the current pipeline state:
data vs TOTAL MC (POT-scaled, **unweighted CV** — not MINERvA Tune v1) plus
the MC background component; no interaction-channel breakdown (that needs a
re-stream with mc_intType / W / Q^2 truth branches).
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

YSCALE = 1e5  # y axis shows value/1e5, labelled "(x10^-5)"


def panel_arrays(h, axis_bin, view):
    """Counts of one panel and the matching (other-axis) bin widths/areas."""
    dpt = np.diff(binning.PT_EDGES_GEV)
    dpl = np.diff(binning.PL_EDGES_GEV)
    if view == "pt":            # panel = pT bin, x axis = p||
        vals = h[axis_bin, :]
        area = dpt[axis_bin] * dpl
        edges = binning.PL_EDGES_GEV
    else:                       # panel = p|| bin, x axis = pT
        vals = h[:, axis_bin]
        area = dpt * dpl[axis_bin]
        edges = binning.PT_EDGES_GEV
    return vals / area, np.sqrt(vals) / area, edges


def make_view(view, data2d, mc2d, bkg2d, outpath, title):
    n_panels = binning.N_PT if view == "pt" else binning.N_PL
    pedges = binning.PT_EDGES_GEV if view == "pt" else binning.PL_EDGES_GEV
    sym = r"p_{t}" if view == "pt" else r"p_{\|}"

    dens_d, errs, dens_m, dens_b, edges = [], [], [], [], None
    for i in range(n_panels):
        d, e, edges = panel_arrays(data2d, i, view)
        m, _, _ = panel_arrays(mc2d, i, view)
        b, _, _ = panel_arrays(bkg2d, i, view)
        dens_d.append(d), errs.append(e), dens_m.append(m), dens_b.append(b)

    ref = max(max(d.max(), m.max()) for d, m in zip(dens_d, dens_m))
    mults = [round(ref / max(d.max(), m.max()), 1) if max(d.max(), m.max()) > 0
             else 1.0 for d, m in zip(dens_d, dens_m)]

    fig, axes = plt.subplots(4, 4, figsize=(13, 8.5), sharey=True)
    axes = axes.ravel()
    centers = 0.5 * (edges[:-1] + edges[1:])
    for i in range(n_panels):
        ax, m = axes[i], mults[i]
        ax.stairs(dens_m[i] * m / YSCALE, edges, color="red", lw=1.8,
                  label="MC total (unweighted CV)")
        ax.stairs(dens_b[i] * m / YSCALE, edges, color="black", lw=1.0,
                  label="MC background")
        ax.errorbar(centers, dens_d[i] * m / YSCALE, yerr=errs[i] * m / YSCALE,
                    fmt="ks", ms=2.5, label="data")
        lo, hi = pedges[i], pedges[i + 1]
        ax.text(0.97, 0.95, rf"${lo:.2f} < {sym} < {hi:.2f}$",
                transform=ax.transAxes, ha="right", va="top", fontsize=8)
        if m != 1.0:
            ax.text(0.97, 0.80, rf"$\times\,{m:g}$", transform=ax.transAxes,
                    ha="right", va="top", fontsize=8)
        if view == "pt":
            ax.set_xscale("log")
            ax.set_xlim(edges[0], edges[-1])
            ax.set_xticks([2, 4, 10, 20, 40, 60])
            ax.set_xticklabels(["2", "4", "10", "20", "40", "60"], fontsize=7)
        else:
            ax.set_xlim(edges[0], edges[-1])
            ax.set_xticks([0, 1, 2, 3, 4])
            ax.tick_params(labelsize=7)
    for j in range(n_panels, 16):
        axes[j].axis("off")
    handles, labels = axes[0].get_legend_handles_labels()
    if n_panels < 16:
        axes[n_panels].legend(handles, labels, loc="center", fontsize=10,
                              frameon=False)
    else:
        fig.legend(handles, labels, loc="upper center", ncol=3, fontsize=9,
                   frameon=False, bbox_to_anchor=(0.5, 0.955))
    fig.supxlabel("Muon Longitudinal Momentum (GeV/c)" if view == "pt"
                  else "Muon Transverse Momentum (GeV/c)", fontsize=12)
    fig.supylabel(r"Events ($\times 10^{-5}$) per (GeV/c)$^2$", fontsize=12)
    fig.suptitle(title, fontsize=11)
    axes[0].set_ylim(0, ref * 1.15 / YSCALE)
    fig.tight_layout(rect=(0.01, 0.0, 1, 0.93 if n_panels == 16 else 0.97))
    fig.savefig(outpath, dpi=160)
    plt.close(fig)
    return mults


def main():
    parser = make_parser("Fig.2-style event-rate multiplier panels from hists.npz.")
    parser.add_argument("--hists", required=True)
    parser.add_argument("--outdir", default=None)
    add_label(parser)
    args = parser.parse_args()

    outdir = Path(args.outdir) if args.outdir else default_outdir(__file__)
    outdir.mkdir(parents=True, exist_ok=True)

    with RunLog(__file__, "Fig.2-style event-rate panels", inputs=args_to_inputs(args)) as log:
        f = np.load(args.hists)
        scale = float(f["pot_data"]) / float(f["pot_mc"])
        data2d = f["data2d"]
        mc2d = f["mc2d_selected"] * scale
        bkg2d = f["mc2d_bkg"] * scale

        m_pt = make_view("pt", data2d, mc2d, bkg2d, outdir / "evtrate_pt_panels.png",
                         "Playlist 1A selected events — panels per $p_T$ bin "
                         "(cf. arXiv:2106.16210 Fig. 2)")
        m_pl = make_view("pl", data2d, mc2d, bkg2d, outdir / "evtrate_pl_panels.png",
                         "Playlist 1A selected events — panels per $p_\\parallel$ bin "
                         "(cf. arXiv:2106.16210 Fig. 2)")

        log.out("outdir", str(outdir))
        log.out("plots", [str(outdir / "evtrate_pt_panels.png"),
                          str(outdir / "evtrate_pl_panels.png")])
        log.out("multipliers", {"pt_view": m_pt, "pl_view": m_pl})
        print(json.dumps({"pt_view_multipliers": m_pt,
                          "pl_view_multipliers": m_pl}, indent=2))


if __name__ == "__main__":
    main()
