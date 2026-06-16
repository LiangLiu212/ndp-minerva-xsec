"""E4 — absolute 2D cross section d²σ/dp_T dp_∥ from an ingredients.npz.

Chain (ExtractCrossSection.cpp:124-130 order), in count-conserving slot space:
  1. background subtraction:  N_sub = data_reco - (POT_data/POT_mc) * bkg
  2. D'Agostini unfold:       reco -> true (prior = MC reco-selected truth = eff_num)
  3. efficiency correction:   / (eff_num / eff_denom)
  4. flux:                    / Φ_int          (Φ in ν/m²/POT)
  5. normalization:           / (N_nucleons * POT_data)
  6. units:                   × 1e4            (m² -> cm²)
  7. project to the 224 measurement cells, then / bin area (Δp_T·Δp_∥)

Stat-only errors: data Poisson propagated through unfold + efficiency (the
0.2 % background's MC stat term is neglected; systematics are the M4 stage).

MC self-closure: feeding the full MC reco (signal+bkg) as pseudo-data at MC POT
must reproduce the directly-normalized MC truth cross section (= effDenom
normalized) — exact, because the migration conserves counts.
"""
import json
from pathlib import Path

import numpy as np

from runlog_tools import (RunLog, add_label, args_to_inputs, default_outdir,
                          make_parser)
from xsec import binning, targets
from xsec.flux import FluxCV
from xsec.unfold import dagostini_unfold


def cell_areas():
    """(224,) bin area Δp_T·Δp_∥ per GlobalID (tutorial edges, GeV²/c²)."""
    gids = np.arange(binning.N_CELLS)
    pt_bin = gids // binning.N_PL          # 0..13
    pl_bin = gids % binning.N_PL           # 0..15
    dpt = np.diff(binning.PT_EDGES_GEV)[pt_bin]
    dpl = np.diff(binning.PL_EDGES_GEV)[pl_bin]
    return dpt * dpl


def extract(data_reco, bkg, migration, eff_num, eff_denom,
            pot_data, pot_mc, flux_integral_m2, n_nucleons, n_iter=10,
            data_var=None):
    """Run the full chain; return (dsigma[224], dsigma_err[224]) in
    cm²/(GeV/c)²/nucleon."""
    bkg_sub = data_reco - (pot_data / pot_mc) * bkg
    if data_var is None:
        data_var = np.abs(data_reco).astype(np.float64)   # Poisson
    unfolded, var, _ = dagostini_unfold(bkg_sub, migration, prior=eff_num,
                                        n_iter=n_iter, data_var=data_var)
    eff = np.divide(eff_num, eff_denom, out=np.zeros_like(eff_num), where=eff_denom > 0)
    eff_corr = np.divide(unfolded, eff, out=np.zeros_like(unfolded), where=eff > 0)
    eff_corr_var = np.divide(var, eff ** 2, out=np.zeros_like(var), where=eff > 0)

    norm = 1e4 / flux_integral_m2 / (n_nucleons * pot_data)
    xsec_cell = binning.to_measurement(eff_corr * norm)
    xsec_var = binning.to_measurement(eff_corr_var * norm ** 2)
    areas = cell_areas()
    return xsec_cell / areas, np.sqrt(xsec_var) / areas


def main():
    parser = make_parser("Extract the absolute 2D cross section from an "
                         "ingredients.npz (bkg-sub -> unfold -> eff -> normalize).")
    parser.add_argument("--ingredients", required=True)
    parser.add_argument("--playlist", default="minervame1A")
    parser.add_argument("--n-iter", type=int, default=10)
    parser.add_argument("--outdir", default=None)
    add_label(parser)
    args = parser.parse_args()

    outdir = Path(args.outdir) if args.outdir else default_outdir(__file__)
    outdir.mkdir(parents=True, exist_ok=True)

    ing = np.load(args.ingredients)
    pot_data, pot_mc = float(ing["pot_data"]), float(ing["pot_mc"])
    flux_m2 = FluxCV(args.playlist).integral(0.0, 100.0)        # ν/m²/POT
    n_nuc, _ = targets.tracker_n_nucleons()

    with RunLog(__file__, f"extract d2sigma {args.playlist}",
                inputs={**args_to_inputs(args), "flux_integral_m2": flux_m2,
                        "n_nucleons": n_nuc, "pot_data": pot_data}) as log:
        kw = dict(migration=ing["migration"], eff_num=ing["eff_num"],
                  eff_denom=ing["eff_denom"], pot_mc=pot_mc,
                  flux_integral_m2=flux_m2, n_nucleons=n_nuc, n_iter=args.n_iter)

        dsigma, dsigma_err = extract(ing["data_reco"], ing["bkg"],
                                     pot_data=pot_data, **kw)

        # --- MC self-closure: MC reco as pseudo-data at MC POT == truth xsec ---
        mc_reco = ing["migration"].sum(axis=1) + ing["bkg"]
        mc_xsec, _ = extract(mc_reco, ing["bkg"], pot_data=pot_mc, **kw)
        # direct MC truth xsec = normalized effDenom
        norm_mc = 1e4 / flux_m2 / (n_nuc * pot_mc)
        truth_xsec = binning.to_measurement(ing["eff_denom"] * norm_mc) / cell_areas()
        closure = np.abs(mc_xsec - truth_xsec)[truth_xsec > 0]
        closure_max = float((closure / truth_xsec[truth_xsec > 0]).max())

        np.savez(outdir / "xsec.npz", dsigma=dsigma, dsigma_err=dsigma_err,
                 truth_xsec=truth_xsec, pt_edges=binning.PT_EDGES_GEV,
                 pl_edges=binning.PL_EDGES_GEV, pot_data=pot_data,
                 flux_integral_m2=flux_m2, n_nucleons=n_nuc)

        summary = {
            "playlist": args.playlist, "n_iter": args.n_iter,
            "flux_integral_nu_per_cm2_per_pot": flux_m2 * 1e-4,
            "n_nucleons": n_nuc, "pot_data": pot_data,
            "integrated_xsec_cm2_per_nucleon": float((dsigma * cell_areas()).sum()),
            "n_cells_positive": int((dsigma > 0).sum()),
            "mc_closure_max_rel": closure_max,
            "mc_closure_pass": bool(closure_max < 1e-3),
        }
        (outdir / "summary.json").write_text(json.dumps(summary, indent=2))
        log.out("outdir", str(outdir))
        log.out("xsec", str(outdir / "xsec.npz"))
        log.out("summary", summary)
        print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
