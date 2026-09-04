"""S4 — assemble Cov_energyscale from shifted ingredients and compare to anc.

σ⁺, σ⁻ from the ±δ shifted migration+bkg (make_energyscale.py) with the CV
eff_num/eff_denom/data; Cov_energyscale = pair_covariance(σ⁺, σ⁻). Offline.
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
from compare_published import load_anc  # reuse the anc data_result parser


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
    parser = make_parser("Assemble Cov_energyscale and compare to anc.")
    parser.add_argument("--ingredients", required=True)
    parser.add_argument("--energyscale", required=True)
    parser.add_argument("--xsec", required=True)
    parser.add_argument("--published", default="config/published.json")
    parser.add_argument("--playlist", default="minervame1A")
    parser.add_argument("--outdir", default=None)
    add_label(parser)
    args = parser.parse_args()

    outdir = Path(args.outdir) if args.outdir else default_outdir(__file__)
    outdir.mkdir(parents=True, exist_ok=True)

    ing = np.load(args.ingredients)
    es = np.load(args.energyscale)
    sig = np.load(args.xsec)["dsigma"]
    flux_m2 = FluxCV(args.playlist).integral(0.0, 100.0)
    n_nuc, _ = targets.tracker_n_nucleons()
    kw = dict(eff_num=ing["eff_num"], eff_denom=ing["eff_denom"],
              pot_data=float(ing["pot_data"]), pot_mc=float(ing["pot_mc"]),
              flux_integral_m2=flux_m2, n_nucleons=n_nuc)

    with RunLog(__file__, "assemble Cov_energyscale", inputs=args_to_inputs(args)) as log:
        sp, _ = extract(ing["data_reco"], es["bkg_plus"], es["migration_plus"], **kw)
        sm, _ = extract(ing["data_reco"], es["bkg_minus"], es["migration_minus"], **kw)
        cov = sx.pair_covariance(sp, sm)

        pub = json.loads(Path(args.published).read_text())
        anc_cov = load_anc_cov(Path(pub["anc_dir"]) /
                               "cov_ptpl_minerva_inclusive_6GeV_energyscale.txt")
        rep = sig > 0
        frac = sx.fractional_error(cov, sig)
        ancf = sx.fractional_error(anc_cov, sig)
        r = np.divide(frac, ancf, out=np.zeros_like(frac), where=ancf > 1e-3)

        np.savez(outdir / "cov_energyscale.npz", cov=cov, frac=frac, delta=float(es["delta"]))
        summary = {
            "delta": float(es["delta"]),
            "ours_frac_median": float(np.median(frac[rep])),
            "anc_frac_median": float(np.median(ancf[rep])),
            "ours_over_anc_median": float(np.median(r[r > 0])),
        }
        (outdir / "summary.json").write_text(json.dumps(summary, indent=2))

        fig, ax = plt.subplots(figsize=(5.5, 5.5))
        ax.plot([0, 0.1], [0, 0.1], "k--", lw=0.8)
        ax.scatter(ancf[rep] * 100, frac[rep] * 100, s=10, alpha=0.6)
        ax.set_xlabel("published energy-scale unc. [%]")
        ax.set_ylabel("our energy-scale unc. [%]")
        ax.set_title(f"per-cell muon energy-scale unc. (δ={float(es['delta'])*100:.2f}%)\n"
                     f"median ratio {summary['ours_over_anc_median']:.2f}")
        ax.set_xlim(0, 9); ax.set_ylim(0, 9)
        fig.tight_layout(); fig.savefig(outdir / "energyscale_scatter.png", dpi=160); plt.close(fig)

        log.out("outdir", str(outdir))
        log.out("cov_energyscale", str(outdir / "cov_energyscale.npz"))
        log.out("summary", summary)
        print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
