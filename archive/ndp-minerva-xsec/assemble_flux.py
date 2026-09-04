"""Assemble the shape-resolved flux covariance and validate vs the anc cov_flux.

Per PPFX universe u, the cross section uses universe u's flux everywhere: the
per-event weight U_u/Φ_cv (already in the streamed ingredients) AND its own
integrated flux Φ_u in the denominator (FluxCV.universe_integrals). The flux
covariance is the ν-e-constraint-WEIGHTED sample covariance of the per-universe
cross sections about their weighted mean — MnvVertErrorBand::CalcCovMx with the
constraint SetUnivWgt (xsec.systematics.sample_covariance(..., weights=w)).

This supersedes the S2 normalization model (a flat 3.23 %); it resolves the
per-cell shape and the off-diagonal correlations. Validated against
cov_ptpl_*_flux.txt. Offline — consumes the streamed flux_universes.npz.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from runlog_tools import (RunLog, add_label, args_to_inputs, default_outdir,
                          make_parser)
from xsec import binning, targets
from xsec import systematics as sx
from xsec.flux import FluxCV
from extract_xsec import extract


def load_anc_cov(path):
    C = np.zeros((binning.N_CELLS, binning.N_CELLS))
    with open(path) as f:
        f.readline()
        for line in f:
            p = line.replace(",", " ").split()
            if len(p) >= 3:
                C[int(p[0]), int(p[1])] = float(p[2])
    return C


def main():
    parser = make_parser("Shape-resolved flux covariance; validate vs anc cov_flux.")
    parser.add_argument("--flux-universes", required=True, help="flux_universes.npz")
    parser.add_argument("--ingredients", required=True)
    parser.add_argument("--xsec", required=True)
    parser.add_argument("--published", default="config/published.json")
    parser.add_argument("--playlist", default="minervame1A")
    parser.add_argument("--outdir", default=None)
    add_label(parser)
    args = parser.parse_args()

    outdir = Path(args.outdir) if args.outdir else default_outdir(__file__)
    outdir.mkdir(parents=True, exist_ok=True)

    fu = np.load(args.flux_universes)
    ing = np.load(args.ingredients)
    sig = np.load(args.xsec)["dsigma"]
    n_univ = int(fu["n_universes"])

    flux = FluxCV(args.playlist)
    phi = flux.universe_integrals(0.0, 100.0)[:n_univ]      # per-universe Φ_u (ν/m²/POT)
    w = flux.constraint_weights[:n_univ]                    # ν-e constraint weights
    n_nuc, _ = targets.tracker_n_nucleons()

    with RunLog(__file__, "shape-resolved flux covariance", inputs=args_to_inputs(args)) as log:
        sigs = np.empty((n_univ, binning.N_CELLS))
        for u in range(n_univ):
            m = fu["migration"][u]
            sigs[u], _ = extract(ing["data_reco"], fu["bkg"][u], m,
                                 eff_num=m.sum(0), eff_denom=fu["eff_denom"][u],
                                 pot_data=float(ing["pot_data"]), pot_mc=float(ing["pot_mc"]),
                                 flux_integral_m2=float(phi[u]), n_nucleons=n_nuc)
        cov_flux = sx.sample_covariance(sigs, weights=w)    # CalcCovMx + SetUnivWgt

        pub = json.loads(Path(args.published).read_text())
        anc = load_anc_cov(Path(pub["anc_dir"]) /
                           "cov_ptpl_minerva_inclusive_6GeV_flux.txt")
        rep = sig > 0
        of = sx.fractional_error(cov_flux, sig)
        af = sx.fractional_error(anc, sig)
        r = np.divide(of, af, out=np.zeros_like(of), where=af > 1e-3)
        # off-diagonal structure: compare correlation matrices on reported cells
        cc_ours = sx.correlation_matrix(cov_flux)[np.ix_(rep, rep)]
        cc_anc = sx.correlation_matrix(anc)[np.ix_(rep, rep)]
        iu = np.triu_indices(rep.sum(), k=1)
        corr_of_corr = float(np.corrcoef(cc_ours[iu], cc_anc[iu])[0, 1])

        np.savez(outdir / "cov_flux.npz", cov_flux=cov_flux, frac_ours=of, frac_anc=af)

        def med(a):
            return float(np.median(a[rep]))
        summary = {
            "n_universes": n_univ,
            "flux_frac_median_ours": med(of),
            "flux_frac_median_anc": med(af),
            "ours_over_anc_median": float(np.median(r[r > 0])),
            "offdiagonal_corr_agreement": corr_of_corr,
            "norm_model_was": 0.0323,
        }
        (outdir / "summary.json").write_text(json.dumps(summary, indent=2))

        fig, ax = plt.subplots(figsize=(5.5, 5.5))
        ax.plot([0, 0.08], [0, 0.08], "k--", lw=0.8)
        ax.scatter(af[rep] * 100, of[rep] * 100, s=10, alpha=0.6)
        ax.set_xlabel("published flux unc. [%]"); ax.set_ylabel("our flux unc. [%]")
        ax.set_title(f"shape-resolved flux ({n_univ} PPFX univ)\n"
                     f"median ours/anc {summary['ours_over_anc_median']:.2f}")
        fig.tight_layout(); fig.savefig(outdir / "flux_scatter.png", dpi=160); plt.close(fig)

        log.out("outdir", str(outdir)); log.out("summary", summary)
        print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
