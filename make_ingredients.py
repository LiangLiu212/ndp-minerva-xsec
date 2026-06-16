"""Stage-3 ingredients producer: the six cross-section ingredients in one npz.

Streams a playlist (xrootd, no local copies) and produces, in count-conserving
flat slot space (xsec/binning.py, N_SLOTS=288 incl. under/overflow):

  - data_reco : selected data, reco slots (unweighted)
  - bkg       : selected MC background (reco-selected, NOT signal), reco slots
  - migration : M[reco_slot, true_slot] for reco-selected SIGNAL
  - eff_num   : efficiency numerator = migration.sum(reco axis) per true slot
                (reco-selected signal, true-binned — free, since no event is
                 dropped in slot space)
  - eff_denom : efficiency DENOMINATOR = signal AND phase-space, true-binned,
                from the MC Truth tree (the genuinely new streaming loop)

MC fills carry the MnvTune v1 CV weight (--weights cv): the RECO side
(bkg, migration/eff_num) gets the full reco weight incl. MINOS efficiency; the
TRUTH side (eff_denom) gets the truth-only weight (no MINOS eff) — matching
Model::GetWeight semantics. Data is never weighted.

Efficiency asymmetry (tutorial Cutter semantics, verified): numerator =
reco-selected signal with NO truth phase-space cut; denominator = signal AND
phase space. eff = num/denom corrects acceptance.

Output: <outdir>/ingredients.npz + summary.json. Downstream extract_xsec.py
projects to the 224 measurement cells via binning.to_measurement.
"""
import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import uproot

from runlog_tools import (RunLog, add_label, args_to_inputs, default_outdir,
                          make_parser)
from xsec import binning, cuts, signal, weights
from xsec.flux import FluxCV
from xsec.kinematics import reco_pt_pz_gev, true_theta_p

DATA_BRANCHES = list(cuts.RECO_SELECTION_BRANCHES) + ["MasterAnaDev_leptonE"]
MC_RECO_BRANCHES = list(dict.fromkeys(
    DATA_BRANCHES + list(signal.SIGNAL_BRANCHES) + list(weights.RECO_WEIGHT_BRANCHES)))
TRUTH_BRANCHES = list(dict.fromkeys(
    list(signal.SIGNAL_BRANCHES) + list(signal.PHASE_SPACE_BRANCHES)
    + list(weights.TRUTH_WEIGHT_BRANCHES)))


def true_pt_pl(mc_primFSLepton):
    lep = np.asarray(mc_primFSLepton, dtype=np.float64)
    theta_t, p_t = true_theta_p(lep[:, 0], lep[:, 1], lep[:, 2])
    return p_t * np.sin(theta_t) / 1000.0, p_t * np.cos(theta_t) / 1000.0


def read_data(url):
    out = {"url": url, "role": "data"}
    with uproot.open(url) as f:
        out["pot_used"] = float(f["Meta"]["POT_Used"].array(library="np").sum())
        arrs = f["MasterAnaDev"].arrays(DATA_BRANCHES, library="np")
    sel = cuts.reco_selection(arrs)
    pt, pl = reco_pt_pz_gev(arrs["MasterAnaDev_leptonE"],
                            arrs["muon_thetaX"], arrs["muon_thetaY"])
    out["data_reco"] = binning.hist_slots(pt[sel], pl[sel])
    out["n_selected"] = int(sel.sum())
    return out


def read_mc(url, weighters):
    """One file open, both trees: reco -> bkg+migration, Truth -> eff_denom."""
    out = {"url": url, "role": "mc"}
    with uproot.open(url) as f:
        out["pot_used"] = float(f["Meta"]["POT_Used"].array(library="np").sum())
        reco = f["MasterAnaDev"].arrays(MC_RECO_BRANCHES, library="np")
        truth = f["Truth"].arrays(TRUTH_BRANCHES, library="np")

    # ---- reco tree: background + migration (signal) ----
    sel = cuts.reco_selection(reco)
    sig = signal.is_signal(reco["mc_incoming"], reco["mc_current"])
    sel_sig, sel_bkg = sel & sig, sel & ~sig
    pt_r, pl_r = reco_pt_pz_gev(reco["MasterAnaDev_leptonE"],
                                reco["muon_thetaX"], reco["muon_thetaY"])
    pt_t, pl_t = true_pt_pl(reco["mc_primFSLepton"])

    if weighters is None:
        w_bkg = w_sig = None
    else:
        w = weights.reco_cv_weight(reco, weighters["flux"], weighters["w2p2h"],
                                   weighters["rpa"])
        w_bkg, w_sig = w[sel_bkg], w[sel_sig]

    out["bkg"] = binning.hist_slots(pt_r[sel_bkg], pl_r[sel_bkg], weights=w_bkg)
    out["migration"] = binning.migration_slots(
        binning.slot_ids(pt_r[sel_sig], pl_r[sel_sig]),
        binning.slot_ids(pt_t[sel_sig], pl_t[sel_sig]), weights=w_sig)
    out["n_signal_sel"] = int(sel_sig.sum())
    out["n_bkg_sel"] = int(sel_bkg.sum())

    # ---- Truth tree: efficiency denominator ----
    denom = signal.is_efficiency_denominator(
        truth["mc_incoming"], truth["mc_current"],
        truth["mc_vtx"], truth["mc_primFSLepton"])
    ptd, pld = true_pt_pl(truth["mc_primFSLepton"])
    if weighters is None:
        wd = None
    else:
        wd = weights.truth_cv_weight(truth, weighters["flux"],
                                     weighters["w2p2h"], weighters["rpa"])
    out["eff_denom"] = binning.hist_slots(ptd[denom], pld[denom],
                                          weights=None if wd is None else wd[denom])
    out["n_denom_unwgt"] = int(denom.sum())
    out["n_truth_entries"] = int(truth["mc_incoming"].size)
    return out


def with_retry(fn, *a, retries=1):
    for attempt in range(retries + 1):
        try:
            return fn(*a)
        except Exception as err:
            if attempt == retries:
                url = a[0]
                return {"url": url, "error": f"{type(err).__name__}: {err}"}
            time.sleep(2.0)


def main():
    parser = make_parser("Stage-3 ingredients (data/bkg/migration/eff num+denom) "
                         "in count-conserving slot space, streaming-only.")
    parser.add_argument("--data-list",
                        default="config/playlists/MediumEnergy_FHC_Data_Playlist1A.txt")
    parser.add_argument("--mc-list",
                        default="config/playlists/MediumEnergy_FHC_StandardMC_Playlist1A.txt")
    parser.add_argument("--max-data-files", type=int, default=None)
    parser.add_argument("--max-mc-files", type=int, default=None)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--weights", choices=["none", "cv"], default="cv")
    parser.add_argument("--playlist", default="minervame1A")
    parser.add_argument("--outdir", default=None)
    add_label(parser)
    args = parser.parse_args()

    data_urls = [u.strip() for u in Path(args.data_list).read_text().splitlines() if u.strip()]
    mc_urls = [u.strip() for u in Path(args.mc_list).read_text().splitlines() if u.strip()]
    if args.max_data_files:
        data_urls = data_urls[:args.max_data_files]
    if args.max_mc_files:
        mc_urls = mc_urls[:args.max_mc_files]

    outdir = Path(args.outdir) if args.outdir else default_outdir(__file__)
    outdir.mkdir(parents=True, exist_ok=True)

    weighters = None
    if args.weights == "cv":
        weighters = {"flux": FluxCV(args.playlist),
                     "w2p2h": weights.TwoP2HWeight(), "rpa": weights.RPAWeight()}
        weighters["rpa"].weight(np.array([0.1]), np.array([0.4]),
                                np.array([1]), np.array([6]))   # pre-warm lru_cache

    with RunLog(__file__, f"ingredients {args.playlist} ({args.weights}): "
                          f"{len(data_urls)} data + {len(mc_urls)} MC",
                inputs={**args_to_inputs(args), "n_data_files": len(data_urls),
                        "n_mc_files": len(mc_urls)}) as log:
        t0 = time.time()
        results, failures = [], []
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futs = {}
            for u in data_urls:
                futs[pool.submit(with_retry, read_data, u)] = "data"
            for u in mc_urls:
                futs[pool.submit(with_retry, read_mc, u, weighters)] = "mc"
            for i, fut in enumerate(as_completed(futs), 1):
                r = fut.result()
                (failures if "error" in r else results).append(r)
                if i % 25 == 0 or i == len(futs):
                    print(f"  [{i}/{len(futs)}] {len(failures)} failures, "
                          f"{time.time()-t0:.0f}s", flush=True)

        def acc(key):
            tot = None
            for r in results:
                if key in r:
                    tot = r[key].copy() if tot is None else tot + r[key]
            return tot

        data_reco = acc("data_reco")
        bkg = acc("bkg")
        migration = acc("migration")
        eff_denom = acc("eff_denom")
        eff_num = migration.sum(axis=0)            # per true slot (free)
        pot_data = sum(r["pot_used"] for r in results if r["role"] == "data")
        pot_mc = sum(r["pot_used"] for r in results if r["role"] == "mc")

        np.savez(outdir / "ingredients.npz",
                 data_reco=data_reco, bkg=bkg, migration=migration,
                 eff_num=eff_num, eff_denom=eff_denom,
                 pt_edges=binning.PT_EDGES_GEV, pl_edges=binning.PL_EDGES_GEV,
                 meas_slots=binning.MEAS_SLOTS, n_slots=binning.N_SLOTS,
                 pot_data=pot_data, pot_mc=pot_mc, weight_mode=np.array(args.weights))

        # gates / summary (projected to the 224 measurement cells)
        num_m = binning.to_measurement(eff_num)
        den_m = binning.to_measurement(eff_denom)
        n_sig_sel = sum(r.get("n_signal_sel", 0) for r in results)
        summary = {
            "playlist": args.playlist, "weight_mode": args.weights,
            "files_ok": {"data": sum(r["role"] == "data" for r in results),
                         "mc": sum(r["role"] == "mc" for r in results)},
            "files_failed": [f["url"] for f in failures],
            "pot": {"data": pot_data, "mc": pot_mc},
            "data_selected_in_grid": float(binning.to_measurement(data_reco).sum()),
            "mc_signal_selected_total": int(n_sig_sel),
            "migration_total": float(migration.sum()),
            "count_conserved": bool(abs(migration.sum()
                                        - sum(r.get("migration", np.zeros(1)).sum()
                                              for r in results if "migration" in r)) < 1e-6),
            "eff_denom_unweighted_total": int(sum(r.get("n_denom_unwgt", 0) for r in results)),
            "average_efficiency_meas": float(num_m.sum() / den_m.sum()) if den_m.sum() else None,
            "wall_s": round(time.time() - t0, 1),
        }
        (outdir / "summary.json").write_text(json.dumps(summary, indent=2))
        log.out("outdir", str(outdir))
        log.out("ingredients", str(outdir / "ingredients.npz"))
        log.out("summary", summary)
        print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
