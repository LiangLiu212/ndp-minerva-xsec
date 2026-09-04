"""S3 — statistical covariance of the 2D cross section via data-Poisson toys.

Throws N Poisson replicas of the slot-space data histogram, runs the full E4
chain per toy (bkg fixed; only the data fluctuates), and takes the sample
covariance of the toy cross sections. This is the data-statistical covariance
including the D'Agostini iteration feedback that E4's analytic `var`
(U²·data_var, final-iteration only) approximates.

Offline — operates on the ingredients.npz; no streaming.
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
from extract_xsec import extract, cell_areas


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
    parser = make_parser("Statistical covariance via data-Poisson toys (offline).")
    parser.add_argument("--ingredients", required=True)
    parser.add_argument("--published", default="config/published.json")
    parser.add_argument("--playlist", default="minervame1A")
    parser.add_argument("--n-toys", type=int, default=1000)
    parser.add_argument("--n-iter", type=int, default=10)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--outdir", default=None)
    add_label(parser)
    args = parser.parse_args()

    outdir = Path(args.outdir) if args.outdir else default_outdir(__file__)
    outdir.mkdir(parents=True, exist_ok=True)

    ing = np.load(args.ingredients)
    pot_data, pot_mc = float(ing["pot_data"]), float(ing["pot_mc"])
    flux_m2 = FluxCV(args.playlist).integral(0.0, 100.0)
    n_nuc, _ = targets.tracker_n_nucleons()
    kw = dict(migration=ing["migration"], eff_num=ing["eff_num"],
              eff_denom=ing["eff_denom"], pot_data=pot_data, pot_mc=pot_mc,
              flux_integral_m2=flux_m2, n_nucleons=n_nuc, n_iter=args.n_iter)

    with RunLog(__file__, f"stat covariance: {args.n_toys} toys, {args.playlist}",
                inputs={**args_to_inputs(args), "flux_integral_m2": flux_m2}) as log:
        # CV result + E4 analytic error (data Poisson var = counts)
        dsigma, dsigma_err = extract(ing["data_reco"], ing["bkg"],
                                     data_var=np.abs(ing["data_reco"]), **kw)

        # Poisson toys (only the data fluctuates)
        rng = np.random.default_rng(args.seed)
        sig = np.zeros((args.n_toys, binning.N_CELLS))
        for t in range(args.n_toys):
            data_toy = rng.poisson(np.clip(ing["data_reco"], 0, None)).astype(np.float64)
            sig[t], _ = extract(data_toy, ing["bkg"], **kw)

        cov_stat = sx.sample_covariance(sig)
        toy_err = np.sqrt(np.clip(np.diag(cov_stat), 0, None))

        # comparison vs the published anc stat covariance (POT-scaled)
        pub = json.loads(Path(args.published).read_text())
        anc_cov = load_anc_cov(Path(pub["anc_dir"]) /
                               "cov_ptpl_minerva_inclusive_6GeV_stat.txt")
        anc_stat_err = np.sqrt(np.clip(np.diag(anc_cov), 0, None))
        # anc stat is full-dataset (small); our 1A stat is LARGER by
        # sqrt(POT_full / POT_1A). Scale the anc UP to the 1A exposure.
        pot_scale = np.sqrt(pub["reference_scalars"]["pot_e20"] * 1e20 / pot_data)
        rep = dsigma > 0
        anc_scaled = anc_stat_err * pot_scale          # anc stat at 1A exposure

        np.savez(outdir / "cov_stat.npz", cov_stat=cov_stat, dsigma=dsigma,
                 toy_err=toy_err, analytic_err=dsigma_err, pot_data=pot_data)

        def med_ratio(a, b, m):
            r = np.divide(a, b, out=np.zeros_like(a), where=(b > 0) & m)
            return float(np.median(r[r > 0]))

        summary = {
            "n_toys": args.n_toys, "playlist": args.playlist,
            "n_reported_cells": int(rep.sum()),
            "toy_vs_analytic_median": med_ratio(toy_err, dsigma_err, rep),
            "toy_frac_stat_median": float(np.median((toy_err / dsigma)[rep])),
            "anc_scaled_frac_median": float(np.median((anc_scaled / dsigma)[rep])),
            "toy_vs_ancscaled_median": med_ratio(toy_err, anc_scaled, rep),
            "pot_scale_full_over_1A": float(pot_scale),
        }
        (outdir / "summary.json").write_text(json.dumps(summary, indent=2))

        # stat-uncertainty map (per-cell fractional)
        frac = np.full((binning.N_PT, binning.N_PL), np.nan)
        with np.errstate(invalid="ignore", divide="ignore"):
            frac_cell = np.where(dsigma > 0, toy_err / dsigma, np.nan)
        for g in np.where(rep)[0]:
            frac[g // binning.N_PL, g % binning.N_PL] = frac_cell[g]
        fig, ax = plt.subplots(figsize=(8, 5))
        pcm = ax.pcolormesh(binning.PL_EDGES_GEV, binning.PT_EDGES_GEV,
                            np.ma.masked_invalid(frac) * 100, cmap="viridis")
        ax.set_xscale("log"); ax.set_xlabel(r"$p_{\parallel,\mu}$ [GeV/c]")
        ax.set_ylabel(r"$p_{T,\mu}$ [GeV/c]")
        ax.set_title(f"toy statistical uncertainty (1A) [%] — {args.n_toys} toys")
        fig.colorbar(pcm, ax=ax, label="fractional stat. unc. [%]")
        fig.tight_layout(); fig.savefig(outdir / "stat_map.png", dpi=160); plt.close(fig)

        log.out("outdir", str(outdir))
        log.out("cov_stat", str(outdir / "cov_stat.npz"))
        log.out("summary", summary)
        print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
