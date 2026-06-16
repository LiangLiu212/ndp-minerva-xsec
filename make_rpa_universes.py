"""RPA systematic universes (vertical weight-matrix streaming pass).

The RPA tune band is two ±1σ pairs (GetRPASystematicsMap, MnvTuneSystematics.cxx
:119-127): "HighQ2" (enhancement, variations 1/2) and "LowQ2" (suppression,
variations 3/4) — 4 universes. RPA is vertical (weight-only), gated on true-QE
(mc_intType==1) on Z>=6; each universe swaps the CV RPA weight for its band
variation (weights.RPAWeight.variation_ratio). One streaming pass (reco + Truth)
fills the 4 per-universe ingredient sets.

cov_RPA = pair_covariance(HighQ2⁺,HighQ2⁻) + pair_covariance(LowQ2⁺,LowQ2⁻)
(each 2-universe band -> CalcCovMx pair), assembled by assemble_total.py --rpa.
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

# (band, sign) order; labels stored in the npz for the assembler.
RPA_UNIVERSES = [("HighQ2", 1.0), ("HighQ2", -1.0), ("LowQ2", 1.0), ("LowQ2", -1.0)]
RPA_LABELS = ["HighQ2_p", "HighQ2_m", "LowQ2_p", "LowQ2_m"]

RECO_BR = list(dict.fromkeys(list(cuts.RECO_SELECTION_BRANCHES)
               + ["MasterAnaDev_leptonE"] + list(signal.SIGNAL_BRANCHES)
               + list(weights.RECO_WEIGHT_BRANCHES)))
TRUTH_BR = list(dict.fromkeys(list(signal.SIGNAL_BRANCHES)
                + list(signal.PHASE_SPACE_BRANCHES)
                + list(weights.TRUTH_WEIGHT_BRANCHES)))


def true_pt_pl(lep):
    lep = np.asarray(lep, np.float64)
    th, p = true_theta_p(lep[:, 0], lep[:, 1], lep[:, 2])
    return p * np.sin(th) / 1000.0, p * np.cos(th) / 1000.0


def variation_weight_columns(arrs, base_w, rpa):
    """(4, n) per-universe weights = base CV weight × (RPA variation / CV) for
    HighQ2±, LowQ2± in RPA_UNIVERSES order."""
    q0, q3 = weights.truth_q0q3_gev(arrs["mc_incomingE"], arrs["mc_primFSLepton"],
                                    arrs["mc_Q2"])
    cols = np.empty((len(RPA_UNIVERSES), base_w.size))
    for i, (band, sign) in enumerate(RPA_UNIVERSES):
        cols[i] = base_w * rpa.variation_ratio(q0, q3, arrs["mc_intType"],
                                               arrs["mc_targetZ"], band, sign)
    return cols


def read_mc(url, wts):
    out = {"url": url}
    with uproot.open(url) as f:
        out["pot_used"] = float(f["Meta"]["POT_Used"].array(library="np").sum())
        a = f["MasterAnaDev"].arrays(RECO_BR, library="np")
        t = f["Truth"].arrays(TRUTH_BR, library="np")
    nU = len(RPA_UNIVERSES)
    sel = cuts.reco_selection(a)
    sig = signal.is_signal(a["mc_incoming"], a["mc_current"])
    ss, sb = sel & sig, sel & ~sig
    pt, pl = reco_pt_pz_gev(a["MasterAnaDev_leptonE"], a["muon_thetaX"], a["muon_thetaY"])
    ptt, plt = true_pt_pl(a["mc_primFSLepton"])
    base = weights.reco_cv_weight(a, wts["flux"], wts["w2p2h"], wts["rpa"])
    reco_s = binning.slot_ids(pt[ss], pl[ss]); true_s = binning.slot_ids(ptt[ss], plt[ss])
    wU = variation_weight_columns({k: a[k] for k in RECO_BR}, base, wts["rpa"])
    mig = np.zeros((nU, binning.N_SLOTS, binning.N_SLOTS))
    bkg = np.zeros((nU, binning.N_SLOTS))
    for u in range(nU):
        mig[u] = binning.migration_slots(reco_s, true_s, weights=wU[u][ss])
        bkg[u] = binning.hist_slots(pt[sb], pl[sb], weights=wU[u][sb])
    out["migration"], out["bkg"] = mig, bkg
    denom = signal.is_efficiency_denominator(t["mc_incoming"], t["mc_current"],
                                             t["mc_vtx"], t["mc_primFSLepton"])
    ptd, pld = true_pt_pl(t["mc_primFSLepton"])
    baset = weights.truth_cv_weight(t, wts["flux"], wts["w2p2h"], wts["rpa"])
    wUt = variation_weight_columns({k: t[k] for k in TRUTH_BR}, baset, wts["rpa"])
    den = np.zeros((nU, binning.N_SLOTS))
    dsl = binning.slot_ids(ptd[denom], pld[denom])
    for u in range(nU):
        den[u], _ = np.histogram(dsl, bins=np.arange(binning.N_SLOTS + 1) - 0.5,
                                 weights=wUt[u][denom])
    out["eff_denom"] = den
    return out


def with_retry(url, wts, retries=1):
    for attempt in range(retries + 1):
        try:
            return read_mc(url, wts)
        except Exception as err:
            if attempt == retries:
                return {"url": url, "error": f"{type(err).__name__}: {err}"}
            time.sleep(2.0)


def main():
    parser = make_parser("RPA systematic universes (4-universe vertical pass).")
    parser.add_argument("--mc-list",
                        default="config/playlists/MediumEnergy_FHC_StandardMC_Playlist1A.txt")
    parser.add_argument("--playlist", default="minervame1A")
    parser.add_argument("--max-mc-files", type=int, default=None)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--outdir", default=None)
    add_label(parser)
    args = parser.parse_args()

    mc = [u.strip() for u in Path(args.mc_list).read_text().splitlines() if u.strip()]
    if args.max_mc_files:
        mc = mc[:args.max_mc_files]
    outdir = Path(args.outdir) if args.outdir else default_outdir(__file__)
    outdir.mkdir(parents=True, exist_ok=True)
    wts = {"flux": FluxCV(args.playlist), "w2p2h": weights.TwoP2HWeight(), "rpa": weights.RPAWeight()}
    wts["rpa"].weight(np.array([0.1]), np.array([0.4]), np.array([1]), np.array([6]))  # prime cache

    with RunLog(__file__, f"RPA universes ({len(RPA_UNIVERSES)}), {args.playlist}",
                inputs={**args_to_inputs(args), "n_mc_files": len(mc), "universes": RPA_LABELS}) as log:
        t0 = time.time(); fail = []; n_ok = 0; pot = 0.0
        tot = {"migration": None, "bkg": None, "eff_denom": None}
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futs = {pool.submit(with_retry, u, wts): u for u in mc}
            for i, fut in enumerate(as_completed(futs), 1):
                r = fut.result()
                if "error" in r:
                    fail.append(r)
                else:                       # incremental accumulation (low memory)
                    n_ok += 1; pot += r["pot_used"]
                    for k in tot:
                        tot[k] = r[k].copy() if tot[k] is None else tot[k] + r[k]
                    r.clear()
                if i % 5 == 0 or i == len(mc):
                    print(f"  [{i}/{len(mc)}] {len(fail)} fail, {time.time()-t0:.0f}s", flush=True)
        np.savez(outdir / "rpa_universes.npz", migration=tot["migration"],
                 bkg=tot["bkg"], eff_denom=tot["eff_denom"],
                 universes=np.array(RPA_LABELS), pot_mc=pot)
        summary = {"playlist": args.playlist, "n_universes": len(RPA_UNIVERSES),
                   "files_ok": n_ok, "files_failed": [f["url"] for f in fail],
                   "wall_s": round(time.time() - t0, 1)}
        (outdir / "summary.json").write_text(json.dumps(summary, indent=2))
        log.out("outdir", str(outdir)); log.out("summary", summary)
        print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
