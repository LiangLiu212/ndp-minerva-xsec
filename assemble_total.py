"""S5 — assemble Cov_total from all systematic groups + stat; compare to anc.

Groups: flux normalization (S2 model), muon energy scale (S4 shifted), GENIE
(S5, 28 knobs × ±1σ), plus the data-statistical covariance (S3). Each group is
built from its per-universe cross sections (extract per universe), summed by
systematics.total_covariance, and the diagonal is compared to the published
cov_total. Offline (consumes the streamed npz products).
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
    parser = make_parser("Assemble Cov_total and compare to published cov_total.")
    parser.add_argument("--ingredients", required=True)
    parser.add_argument("--xsec", required=True)
    parser.add_argument("--energyscale", required=True)
    parser.add_argument("--stat-cov", required=True)
    parser.add_argument("--genie", required=True)
    parser.add_argument("--flux-frac", type=float, default=0.0323)
    parser.add_argument("--published", default="config/published.json")
    parser.add_argument("--playlist", default="minervame1A")
    parser.add_argument("--outdir", default=None)
    add_label(parser)
    args = parser.parse_args()

    outdir = Path(args.outdir) if args.outdir else default_outdir(__file__)
    outdir.mkdir(parents=True, exist_ok=True)

    ing = np.load(args.ingredients)
    sig = np.load(args.xsec)["dsigma"]
    flux_m2 = FluxCV(args.playlist).integral(0.0, 100.0)
    n_nuc, _ = targets.tracker_n_nucleons()
    kw = dict(pot_data=float(ing["pot_data"]), pot_mc=float(ing["pot_mc"]),
              flux_integral_m2=flux_m2, n_nucleons=n_nuc)

    with RunLog(__file__, "assemble Cov_total", inputs=args_to_inputs(args)) as log:
        # flux (normalization model, S2)
        cov_flux = sx.normalization_covariance(sig, args.flux_frac)

        # energy scale (S4)
        es = np.load(args.energyscale)
        sp, _ = extract(ing["data_reco"], es["bkg_plus"], es["migration_plus"],
                        eff_num=ing["eff_num"], eff_denom=ing["eff_denom"], **kw)
        sm, _ = extract(ing["data_reco"], es["bkg_minus"], es["migration_minus"],
                        eff_num=ing["eff_num"], eff_denom=ing["eff_denom"], **kw)
        cov_es = sx.pair_covariance(sp, sm)

        # stat (S3)
        cov_stat = np.load(args.stat_cov)["cov_stat"]

        # GENIE (S5): per knob ±1σ pair, summed
        g = np.load(args.genie)
        nk = len(g["knobs"])
        cov_genie = np.zeros((binning.N_CELLS, binning.N_CELLS))
        for i in range(nk):
            mp, mm = g["migration"][2 * i], g["migration"][2 * i + 1]
            sp_i, _ = extract(ing["data_reco"], g["bkg"][2 * i], mp,
                              eff_num=mp.sum(0), eff_denom=g["eff_denom"][2 * i], **kw)
            sm_i, _ = extract(ing["data_reco"], g["bkg"][2 * i + 1], mm,
                              eff_num=mm.sum(0), eff_denom=g["eff_denom"][2 * i + 1], **kw)
            cov_genie += sx.pair_covariance(sp_i, sm_i)

        cov_total = sx.total_covariance([cov_flux, cov_es, cov_genie], cov_stat)

        pub = json.loads(Path(args.published).read_text())
        anc = load_anc_cov(Path(pub["anc_dir"]) /
                           "cov_ptpl_minerva_inclusive_6GeV_total.txt")
        rep = sig > 0
        ours = sx.fractional_error(cov_total, sig)
        ancf = sx.fractional_error(anc, sig)
        # the anc total is full-dataset; scale our stat down to compare totals
        # fairly, OR report the systematic-only comparison too.
        cov_sys = sx.total_covariance([cov_flux, cov_es, cov_genie])
        ours_sys = sx.fractional_error(cov_sys, sig)

        np.savez(outdir / "cov_total.npz", cov_total=cov_total, cov_flux=cov_flux,
                 cov_energyscale=cov_es, cov_genie=cov_genie, cov_stat=cov_stat,
                 frac_total=ours, frac_anc=ancf)

        def med(a, m=rep):
            return float(np.median(a[m]))
        r = np.divide(ours, ancf, out=np.zeros_like(ours), where=ancf > 1e-3)
        summary = {
            "groups_frac_median": {
                "flux": med(sx.fractional_error(cov_flux, sig)),
                "energy_scale": med(sx.fractional_error(cov_es, sig)),
                "genie": med(sx.fractional_error(cov_genie, sig)),
                "stat_1A": med(sx.fractional_error(cov_stat, sig)),
            },
            "total_frac_median_ours": med(ours),
            "systematic_only_frac_median_ours": med(ours_sys),
            "total_frac_median_anc": med(ancf),
            "ours_over_anc_median": float(np.median(r[r > 0])),
            "note": "ours total includes 1A stat (~4.7%, larger than full-dataset); "
                    "missing groups: 2p2h, RPA, MINOS-eff, Geant-hadron, beam angle, "
                    "muon resolution, and the flux shape term.",
        }
        (outdir / "summary.json").write_text(json.dumps(summary, indent=2))

        fig, ax = plt.subplots(figsize=(5.5, 5.5))
        ax.plot([0, 0.15], [0, 0.15], "k--", lw=0.8)
        ax.scatter(ancf[rep] * 100, ours[rep] * 100, s=10, alpha=0.6)
        ax.set_xlabel("published total unc. [%]"); ax.set_ylabel("our total unc. [%]")
        ax.set_title(f"per-cell total uncertainty\nmedian ours/anc {summary['ours_over_anc_median']:.2f}")
        fig.tight_layout(); fig.savefig(outdir / "total_scatter.png", dpi=160); plt.close(fig)

        log.out("outdir", str(outdir)); log.out("summary", summary)
        print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
