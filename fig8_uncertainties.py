"""Fig-8-style fractional-uncertainty breakdown of the 2D cross-section budget.

Reproduces arXiv:2106.16210 Fig. 8 from an assembled cov_total.npz: two
small-multiple grids of per-(p_T, p_||)-cell fractional uncertainty —
  errors-pt : 14 panels (one per p_T bin),  x = p_||  (log)
  errors-pz : 16 panels (one per p_|| bin), x = p_T   (bin index)

Curves map our groups onto the paper's seven:
  Total / Statistical / Flux / Models / Muon Reconstruction / Normalization.
Hadronic Response is not yet built -> omitted (flagged in the title).
Normalization is the exact flat band (target-nucleon count, fully correlated).
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from runlog_tools import (RunLog, add_label, args_to_inputs, default_outdir,
                          make_parser)
from xsec import binning


def fractional(cov, sig, rep):
    """Per-cell fractional uncertainty sqrt(diag)/sigma, 0 where not reported."""
    f = np.zeros_like(sig)
    dg = np.clip(np.diag(cov), 0.0, None)
    f[rep] = np.sqrt(dg[rep]) / sig[rep]
    return f


def main():
    parser = make_parser("Fig-8-style fractional-uncertainty breakdown of the 2D budget.")
    parser.add_argument("--cov", required=True,
                        help="assembled cov_total.npz (per-group covariances)")
    parser.add_argument("--xsec", required=True,
                        help="xsec.npz with the CV dsigma (fractional denominator)")
    parser.add_argument("--norm-frac", type=float, default=0.014,
                        help="target-nucleon-count normalization band (flat, fully correlated)")
    parser.add_argument("--outdir", default=None)
    add_label(parser)
    args = parser.parse_args()

    outdir = Path(args.outdir) if args.outdir else default_outdir(__file__)
    outdir.mkdir(parents=True, exist_ok=True)

    cov_path, xs_path = Path(args.cov).resolve(), Path(args.xsec).resolve()
    with RunLog(__file__, "Fig-8 fractional-uncertainty breakdown",
                inputs={**args_to_inputs(args),
                        "cov": str(cov_path), "xsec": str(xs_path)}) as log:
        N_PT, N_PL = binning.N_PT, binning.N_PL
        PT_E, PL_E = binning.PT_EDGES_GEV, binning.PL_EDGES_GEV
        PL_CEN = 0.5 * (PL_E[:-1] + PL_E[1:])

        d = np.load(cov_path)
        sig = np.load(xs_path)["dsigma"]
        rep = sig > 0                  # reported cells; rest are acceptance gaps

        cov_models = (d["cov_genie"] + d["cov_twop2h"] + d["cov_rpa"]
                      + d["cov_geniervx1pi"])
        f_flux = fractional(d["cov_flux"], sig, rep)
        f_muon = fractional(d["cov_energyscale"], sig, rep)   # full muon-reco group
        f_models = fractional(cov_models, sig, rep)
        f_stat = fractional(d["cov_stat"], sig, rep)          # 1A only -> inflated
        f_norm = np.where(rep, args.norm_frac, 0.0)
        # total = assembled cov_total + the exact normalization band, in quadrature.
        # Only Hadronic Response (~1%) is still missing from this black curve.
        f_total = np.sqrt(fractional(d["cov_total"], sig, rep) ** 2 + f_norm ** 2)

        curves = [
            ("Total Uncertainty",   f_total,  dict(color="black",   ls="-",  lw=1.7)),
            ("Statistical (1A)",     f_stat,   dict(color="0.45",    ls="--", lw=1.5)),
            ("Flux",                 f_flux,   dict(color="#9b7fd4", ls="-",  lw=1.6)),
            ("Models",               f_models, dict(color="#e41a1c", ls="-",  lw=1.6)),
            ("Muon Reconstruction",  f_muon,   dict(color="#ff8c00", ls="-",  lw=1.6)),
            ("Normalization",        f_norm,   dict(color="#1f47d6", ls="-",  lw=1.6)),
        ]
        repm = rep.reshape(N_PT, N_PL)

        def masked(fvec):
            """(N_PT, N_PL) with NaN in acceptance gaps so lines break like the paper."""
            m = fvec.reshape(N_PT, N_PL).astype(float).copy()
            m[~repm] = np.nan
            return m

        cmats = [(lab, masked(fv), st) for lab, fv, st in curves]
        handles = [plt.Line2D([0], [0], **st) for _, _, st in curves]
        labels = [lab for lab, _, _ in curves]
        sup = ("Fractional uncertainties of  d$^2\\sigma$/d$p_T$d$p_\\parallel$  "
               "— playlist 1A, Fig. 8 format\n"
               "Hadronic Response not yet built (omitted);  Statistical is 1A-only "
               "(inflated vs full dataset)")

        def style_axis(ax):
            ax.set_ylim(0, 0.30)
            ax.set_yticks([0, 0.1, 0.2, 0.3])
            ax.tick_params(labelsize=8)

        # ---- Figure A: errors-pt (panels by p_T, x = p_||, log) ------------
        figA, axA = plt.subplots(4, 4, figsize=(13, 9), sharex=True, sharey=True)
        axA = axA.ravel()
        for ipt in range(N_PT):
            ax = axA[ipt]
            for _, m, st in cmats:
                ax.plot(PL_CEN, m[ipt], **st)
            ax.set_xscale("log")
            ax.set_xlim(1.5, 60)
            ax.set_xticks([4, 10, 20, 40, 60])
            ax.set_xticklabels(["4", "10", "20", "40", "60"])
            ax.get_xaxis().set_minor_formatter(plt.NullFormatter())
            style_axis(ax)
            ax.text(0.96, 0.93, f"{PT_E[ipt]:.2f} < $p_T$ < {PT_E[ipt + 1]:.2f}",
                    transform=ax.transAxes, ha="right", va="top", fontsize=8)
        for k in range(N_PT, 16):
            axA[k].axis("off")
        axA[15].legend(handles, labels, loc="center", fontsize=10, frameon=False)
        figA.text(0.5, 0.04, "Muon Longitudinal Momentum (GeV/c)", ha="center", fontsize=13)
        figA.text(0.04, 0.5, "Fractional uncertainty", va="center",
                  rotation="vertical", fontsize=13)
        figA.suptitle(sup, fontsize=11)
        figA.subplots_adjust(left=0.08, right=0.98, top=0.90, bottom=0.09,
                             hspace=0.06, wspace=0.06)
        pt_png = outdir / "errors_pt_fig8.png"
        figA.savefig(pt_png, dpi=160)
        plt.close(figA)

        # ---- Figure B: errors-pz (panels by p_||, x = p_T, bin index) ------
        xb = np.arange(N_PT)
        pt_tick_val = [0.40, 1.00, 1.50, 2.50, 4.50]
        pt_tick_pos = [int(np.where(np.isclose(PT_E, v))[0][0]) - 0.5 for v in pt_tick_val]
        figB, axB = plt.subplots(4, 4, figsize=(13, 9), sharex=True, sharey=True)
        axB = axB.ravel()
        for ipl in range(N_PL):
            ax = axB[ipl]
            for _, m, st in cmats:
                ax.plot(xb, m[:, ipl], **st)
            ax.set_xlim(-0.5, 13.5)
            ax.set_xticks(pt_tick_pos)
            ax.set_xticklabels([f"{v:.1f}" for v in pt_tick_val])
            style_axis(ax)
            ax.text(0.96, 0.93, f"{PL_E[ipl]:.1f} < $p_\\parallel$ < {PL_E[ipl + 1]:.1f}",
                    transform=ax.transAxes, ha="right", va="top", fontsize=8)
        figB.legend(handles, labels, loc="lower center", ncol=6, fontsize=10,
                    frameon=False, bbox_to_anchor=(0.5, 0.0))
        figB.text(0.5, 0.05, "Muon Transverse Momentum (GeV/c)", ha="center", fontsize=13)
        figB.text(0.04, 0.5, "Fractional uncertainty", va="center",
                  rotation="vertical", fontsize=13)
        figB.suptitle(sup, fontsize=11)
        figB.subplots_adjust(left=0.08, right=0.98, top=0.90, bottom=0.11,
                             hspace=0.06, wspace=0.06)
        pz_png = outdir / "errors_pz_fig8.png"
        figB.savefig(pz_png, dpi=160)
        plt.close(figB)

        medians = {lab: float(np.median(fv[rep])) for lab, fv, _ in curves}
        log.out("figures", [str(pt_png), str(pz_png)])
        log.out("median_frac", medians)
        log.out("n_reported", int(rep.sum()))
        print("median frac uncertainty over %d reported cells:" % rep.sum())
        for lab, mv in medians.items():
            print(f"  {lab:22s} {mv * 100:5.2f} %")
        print("wrote:", pt_png, "and", pz_png)


if __name__ == "__main__":
    main()
