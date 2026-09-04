"""Assemble the muon-reconstruction covariance from the band universes.

Each band is a ±1σ pair → pair_covariance(σ⁺,σ⁻). Two products:
  * cov_energyscale = Muon_Energy_MINERvA ⊕ Muon_Energy_MINOS — validated
    against the anc cov_energyscale file (the only muon piece with a dedicated
    anc covariance);
  * cov_muon = all eight bands — the full Fig 8 "Muon Reconstruction" category,
    folded into Cov_total by assemble_total.py --muon.

eff_denom is the CV one (every muon band is a reco-side shift); eff_num per
universe is the migration column sum (re-selection changes it for angle bands).
Offline — consumes the streamed muon_universes.npz.
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

ENERGY_SCALE_BANDS = ["Muon_Energy_MINERvA", "Muon_Energy_MINOS"]


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
    parser = make_parser("Assemble muon-reconstruction covariance; validate energy scale vs anc.")
    parser.add_argument("--muon", required=True, help="muon_universes.npz")
    parser.add_argument("--ingredients", required=True)
    parser.add_argument("--xsec", required=True)
    parser.add_argument("--published", default="config/published.json")
    parser.add_argument("--playlist", default="minervame1A")
    parser.add_argument("--outdir", default=None)
    add_label(parser)
    args = parser.parse_args()

    outdir = Path(args.outdir) if args.outdir else default_outdir(__file__)
    outdir.mkdir(parents=True, exist_ok=True)

    mu = np.load(args.muon)
    ing = np.load(args.ingredients)
    sig = np.load(args.xsec)["dsigma"]
    bands = [b for b in mu["bands"]]
    flux_m2 = FluxCV(args.playlist).integral(0.0, 100.0)
    n_nuc, _ = targets.tracker_n_nucleons()
    kw = dict(pot_data=float(ing["pot_data"]), pot_mc=float(ing["pot_mc"]),
              flux_integral_m2=flux_m2, n_nucleons=n_nuc)

    def band_cov(band):
        out = {}
        for tag in ("plus", "minus"):
            m = mu[f"{band}__mig__{tag}"]
            out[tag], _ = extract(ing["data_reco"], mu[f"{band}__bkg__{tag}"], m,
                                  eff_num=m.sum(0), eff_denom=ing["eff_denom"], **kw)
        return sx.pair_covariance(out["plus"], out["minus"])

    with RunLog(__file__, "assemble muon-reco covariance", inputs=args_to_inputs(args)) as log:
        covs = {b: band_cov(b) for b in bands}
        cov_es = sx.total_covariance([covs[b] for b in ENERGY_SCALE_BANDS])
        cov_muon = sx.total_covariance(list(covs.values()))

        pub = json.loads(Path(args.published).read_text())
        anc = load_anc_cov(Path(pub["anc_dir"]) /
                           "cov_ptpl_minerva_inclusive_6GeV_energyscale.txt")
        rep = sig > 0
        es_f = sx.fractional_error(cov_es, sig)
        anc_f = sx.fractional_error(anc, sig)
        r = np.divide(es_f, anc_f, out=np.zeros_like(es_f), where=anc_f > 1e-3)

        np.savez(outdir / "cov_muon.npz", cov_muon=cov_muon, cov_energyscale=cov_es,
                 **{f"cov_{b}": covs[b] for b in bands})

        def med(a):
            return float(np.median(a[rep]))
        summary = {
            "bands_frac_median": {b: med(sx.fractional_error(covs[b], sig)) for b in bands},
            "energy_scale_frac_median_ours": med(es_f),
            "energy_scale_frac_median_anc": med(anc_f),
            "energy_scale_ours_over_anc_median": float(np.median(r[r > 0])),
            "muon_reco_total_frac_median": med(sx.fractional_error(cov_muon, sig)),
        }
        (outdir / "summary.json").write_text(json.dumps(summary, indent=2))

        fig, ax = plt.subplots(figsize=(5.5, 5.5))
        ax.plot([0, 0.12], [0, 0.12], "k--", lw=0.8)
        ax.scatter(anc_f[rep] * 100, es_f[rep] * 100, s=10, alpha=0.6)
        ax.set_xlabel("published energy-scale unc. [%]")
        ax.set_ylabel("our energy-scale unc. [%]")
        ax.set_title("muon energy scale (MINERvA ⊕ MINOS)\n"
                     f"median ours/anc {summary['energy_scale_ours_over_anc_median']:.2f}")
        fig.tight_layout(); fig.savefig(outdir / "energyscale_scatter.png", dpi=160); plt.close(fig)

        log.out("outdir", str(outdir)); log.out("summary", summary)
        print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
