"""2D D'Agostini unfolding via the actual RooUnfold library (PyROOT), the same
method as MnvUnfold::UnfoldHisto2D:

    RooUnfoldResponse response(reco, true, migration);
    RooUnfoldBayes    bayes(&response, data, n_iter);
    h_unfold = bayes.Hreco();   cov = bayes.Ereco();

RooUnfold unfolds a 2D measurement by flattening the (p_T, p_||) grid to a 1D
global-bin index internally (RooUnfoldResponse::FindBin = FindFixBin-1,
RooUnfoldResponse.cxx:536). We do that flattening explicitly onto the pipeline's
288 count-conserving slots (under/overflow included) and hand RooUnfold the
flattened histograms. This is identical math to the TH2D interface but avoids
RooUnfold's acknowledged-buggy 2D-overflow path (the "TODO this doesn't work for
overflows" at RooUnfoldResponse.cxx:540 — the reason MnvUnfold's UnfoldHisto2D
runs with UseOverflow off).

Validation, all on playlist 1A:
  * MC self-closure: unfold the MC reco -> recover MC truth (~1e-15).
  * cross-check vs our pure-Python xsec.unfold.dagostini_unfold (identical
    inputs) -> machine precision.
  * carry both through the identical efficiency + normalization to the 224-cell
    d2sigma and compare.

Requires the RooUnfold lib (build once, no env install needed):
    pixi run bash external/roounfold/build_roounfold.sh
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import ROOT

from runlog_tools import (RunLog, add_label, args_to_inputs, default_outdir,
                          make_parser)
from xsec import binning, targets
from xsec.flux import FluxCV
from xsec.unfold import dagostini_unfold
import extract_xsec

NG = binning.N_SLOTS               # 288 count-conserving slots


def th1_slots(name, vec):
    """TH1D(288) with bin i+1 = slot i."""
    h = ROOT.TH1D(name, "", NG, 0, NG)
    for s, v in enumerate(vec):
        h.SetBinContent(s + 1, float(v))
    return h


def migration_th2(name, M):
    """TH2D(288 x 288) response with [reco_slot+1, true_slot+1] = M[reco, true]."""
    h = ROOT.TH2D(name, "", NG, 0, NG, NG, 0, NG)
    rs, ts = np.nonzero(M)
    for r, t in zip(rs, ts):
        h.SetBinContent(int(r) + 1, int(t) + 1, float(M[r, t]))
    return h


_KEEP = []   # hold ROOT objects so PyROOT's GC doesn't free them under RooUnfold


def roo_unfold(resp, data_slots, n_iter, tag):
    """RooUnfoldBayes -> (unfolded slot vector[288], Ereco diagonal[288], nrows).

    Reads everything into numpy before the bayes object can be collected.
    """
    dh = th1_slots("data_" + tag, data_slots)
    ROOT.SetOwnership(dh, False)
    bayes = ROOT.RooUnfoldBayes(resp, dh, int(n_iter))
    bayes.SetVerbose(0)
    h = bayes.Hreco()
    unf = np.array([h.GetBinContent(s + 1) for s in range(NG)])
    ec = bayes.Ereco()
    n = int(ec.GetNrows())
    cov_diag = np.array([ec[s][s] for s in range(NG)]) if n == NG else np.zeros(NG)
    _KEEP.append((dh, bayes))
    return unf, cov_diag, n


def main():
    parser = make_parser("2D D'Agostini unfolding via RooUnfold (MnvUnfold::UnfoldHisto2D).")
    parser.add_argument("--ingredients", required=True, help="make_ingredients .npz")
    parser.add_argument("--lib", default="external/roounfold/libRooUnfoldMin.so")
    parser.add_argument("--n-iter", type=int, default=10)
    parser.add_argument("--playlist", default="minervame1A")
    parser.add_argument("--outdir", default=None)
    add_label(parser)
    args = parser.parse_args()

    outdir = Path(args.outdir) if args.outdir else default_outdir(__file__)
    outdir.mkdir(parents=True, exist_ok=True)

    if ROOT.gSystem.Load(args.lib) < 0 or not hasattr(ROOT, "RooUnfoldBayes"):
        raise SystemExit(f"RooUnfold not loadable from {args.lib}; build it: "
                         "pixi run bash external/roounfold/build_roounfold.sh")
    ROOT.gErrorIgnoreLevel = ROOT.kError

    ing = np.load(args.ingredients)
    migration = ing["migration"].astype(np.float64)
    eff_num, eff_denom = ing["eff_num"], ing["eff_denom"]
    data_reco, bkg = ing["data_reco"].astype(np.float64), ing["bkg"]
    pot_data, pot_mc = float(ing["pot_data"]), float(ing["pot_mc"])

    measured = migration.sum(axis=1)        # MC reco projection (row sums)
    truth = migration.sum(axis=0)           # MC reco-selected truth (== eff_num)
    data_slot = data_reco - (pot_data / pot_mc) * bkg

    with RunLog(__file__, "RooUnfold 2D unfolding (flattened slots)",
                inputs={**args_to_inputs(args),
                        "ingredients": str(Path(args.ingredients).resolve())}) as log:
        hm, ht, hmig = (th1_slots("meas", measured), th1_slots("true", truth),
                        migration_th2("mig", migration))
        for h in (hm, ht, hmig):
            ROOT.SetOwnership(h, False)
        resp = ROOT.RooUnfoldResponse(hm, ht, hmig)

        # (1) MC self-closure: unfold the MC reco -> should recover MC truth
        closure_unf, _, _ = roo_unfold(resp, measured, args.n_iter, "closure")
        self_closure = float(np.max(np.abs(closure_unf - truth)) / max(truth.max(), 1.0))

        # (2) unfold the real (bkg-subtracted) data
        roo_unf, cov_slot, ecov_n = roo_unfold(resp, data_slot, args.n_iter, "data")

        # (3) cross-check vs our pure-Python D'Agostini, identical inputs
        our_unf, _, _ = dagostini_unfold(data_slot, migration, prior=eff_num,
                                         n_iter=args.n_iter, data_var=np.abs(data_reco))
        m = our_unf > 1e-9 * our_unf.max()
        frac = np.abs(roo_unf[m] - our_unf[m]) / our_unf[m]
        xcheck = {"median": float(np.median(frac)), "p95": float(np.percentile(frac, 95)),
                  "max": float(np.max(frac))}

        # (4) carry both through the identical efficiency + normalization to d2sigma
        flux_m2 = FluxCV(args.playlist).integral(0.0, 100.0)
        n_nuc, _ = targets.tracker_n_nucleons()
        eff = np.divide(eff_num, eff_denom, out=np.zeros_like(eff_num), where=eff_denom > 0)
        norm = 1e4 / flux_m2 / (n_nuc * pot_data)
        areas = extract_xsec.cell_areas()
        eff_corr = np.divide(roo_unf, eff, out=np.zeros_like(roo_unf), where=eff > 0)
        roo_xsec = binning.to_measurement(eff_corr * norm) / areas
        our_xsec, _ = extract_xsec.extract(data_reco, bkg, migration, eff_num, eff_denom,
                                           pot_data, pot_mc, flux_m2, n_nuc, n_iter=args.n_iter)
        rep = our_xsec > 0
        xsec_ratio = np.divide(roo_xsec, our_xsec, out=np.ones_like(roo_xsec), where=rep)
        roo_int = float((roo_xsec[rep] * areas[rep]).sum())
        our_int = float((our_xsec[rep] * areas[rep]).sum())

        # RooUnfold's full Ereco -> per-cell fractional unfolding stat uncertainty
        roo_frac = np.divide(np.sqrt(np.clip(cov_slot, 0, None)), roo_unf,
                             out=np.zeros_like(cov_slot), where=roo_unf > 0)
        frac_unc_med = float(np.median(binning.to_measurement(roo_frac)[rep]))

        summary = {
            "n_iter": args.n_iter,
            "self_closure_max_frac": self_closure,
            "xcheck_vs_python_unfold": xcheck,
            "xsec_median_ratio_roo_over_ours": float(np.median(xsec_ratio[rep])),
            "xsec_max_abs_dev": float(np.max(np.abs(xsec_ratio[rep] - 1.0))),
            "integrated_roo_over_ours": roo_int / our_int,
            "unfold_frac_unc_median_cell": frac_unc_med,
            "ereco_nrows": ecov_n,
        }
        np.savez(outdir / "unfold_roounfold.npz",
                 roo_unfolded_slot=roo_unf, our_unfolded_slot=our_unf,
                 roo_xsec=roo_xsec, our_xsec=our_xsec, cov_slot=cov_slot)
        (outdir / "summary.json").write_text(json.dumps(summary, indent=2))

        fig, ax = plt.subplots(1, 2, figsize=(10, 4.6))
        ax[0].plot(our_unf[m], roo_unf[m], ".", ms=4, alpha=0.5)
        lim = [0, our_unf[m].max() * 1.05]
        ax[0].plot(lim, lim, "k--", lw=0.8); ax[0].set_xlim(lim); ax[0].set_ylim(lim)
        ax[0].set_xlabel("our dagostini_unfold [slots]"); ax[0].set_ylabel("RooUnfold [slots]")
        ax[0].set_title(f"unfolded distribution\nmedian frac diff {xcheck['median']:.1e}")
        ax[1].plot(our_xsec[rep], roo_xsec[rep], ".", ms=4, alpha=0.5)
        lim2 = [0, our_xsec[rep].max() * 1.05]
        ax[1].plot(lim2, lim2, "k--", lw=0.8); ax[1].set_xlim(lim2); ax[1].set_ylim(lim2)
        ax[1].set_xlabel("our d2sigma [224 cells]"); ax[1].set_ylabel("RooUnfold d2sigma")
        ax[1].set_title(f"cross section\nmedian ratio {np.median(xsec_ratio[rep]):.4f}")
        fig.tight_layout(); fig.savefig(outdir / "roounfold_vs_ours.png", dpi=150); plt.close(fig)

        log.out("outdir", str(outdir))
        log.out("self_closure_max_frac", self_closure)
        log.out("summary", summary)
        print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
