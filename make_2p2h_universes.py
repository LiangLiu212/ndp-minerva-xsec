"""S5 — 2p2h low-recoil tune systematic (vertical weight-matrix streaming pass).

The MnvTune 2p2h band has THREE universes (Get2p2hSystematics,
MnvTuneSystematics.cxx:64-73), each at nsigma=1 with a different MEC-pair fit:
  mode 1 = nn/pp-only, mode 2 = np-only, mode 3 = QE->2p2h.
2p2h is vertical (weight-only): the reco/true slots are the CV ones; each
universe swaps the CV 2p2h weight for its variation (weights.twop2h_variation_
ratio). One streaming pass (reco + Truth) fills the 3 per-universe ingredient
sets (migration, bkg, eff_denom); eff_num is each migration's column sum.

The covariance is the sample covariance of the three universe cross sections
about their mean (MnvVertErrorBand::CalcCovMx, fUseSpreadError=false for >1
universe) — assembled downstream by assemble_total.py.
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

MODES = [1, 2, 3]                       # nn/pp, np, QE->2p2h (cxx:20-24)
NUCLEON_BR = ["mc_targetNucleon"]       # struck-pair gating (cxx:48-50)

RECO_BR = list(dict.fromkeys(list(cuts.RECO_SELECTION_BRANCHES)
               + ["MasterAnaDev_leptonE"] + list(signal.SIGNAL_BRANCHES)
               + list(weights.RECO_WEIGHT_BRANCHES) + NUCLEON_BR))
TRUTH_BR = list(dict.fromkeys(list(signal.SIGNAL_BRANCHES)
                + list(signal.PHASE_SPACE_BRANCHES)
                + list(weights.TRUTH_WEIGHT_BRANCHES) + NUCLEON_BR))


def true_pt_pl(lep):
    lep = np.asarray(lep, np.float64)
    th, p = true_theta_p(lep[:, 0], lep[:, 1], lep[:, 2])
    return p * np.sin(th) / 1000.0, p * np.cos(th) / 1000.0


def variation_weight_columns(arrs, base_w, w2p2h):
    """(3, n) per-universe weights = base CV weight × (variation / CV) for the
    three 2p2h modes, in MODES order."""
    q0, q3 = weights.truth_q0q3_gev(arrs["mc_incomingE"], arrs["mc_primFSLepton"],
                                    arrs["mc_Q2"])
    cols = np.empty((len(MODES), base_w.size))
    for i, m in enumerate(MODES):
        cols[i] = base_w * weights.twop2h_variation_ratio(
            q0, q3, arrs["mc_intType"], arrs["mc_targetZ"],
            arrs["mc_targetNucleon"], w2p2h, m)
    return cols


def read_mc(url, wts):
    out = {"url": url}
    with uproot.open(url) as f:
        out["pot_used"] = float(f["Meta"]["POT_Used"].array(library="np").sum())
        a = f["MasterAnaDev"].arrays(RECO_BR, library="np")
        t = f["Truth"].arrays(TRUTH_BR, library="np")
    nU = len(MODES)
    # reco: migration + bkg per universe (slots fixed; weights vary)
    sel = cuts.reco_selection(a)
    sig = signal.is_signal(a["mc_incoming"], a["mc_current"])
    ss, sb = sel & sig, sel & ~sig
    pt, pl = reco_pt_pz_gev(a["MasterAnaDev_leptonE"], a["muon_thetaX"], a["muon_thetaY"])
    ptt, plt = true_pt_pl(a["mc_primFSLepton"])
    base = weights.reco_cv_weight(a, wts["flux"], wts["w2p2h"], wts["rpa"])
    reco_s = binning.slot_ids(pt[ss], pl[ss]); true_s = binning.slot_ids(ptt[ss], plt[ss])
    wU = variation_weight_columns({k: a[k] for k in RECO_BR}, base, wts["w2p2h"])
    mig = np.zeros((nU, binning.N_SLOTS, binning.N_SLOTS))
    bkg = np.zeros((nU, binning.N_SLOTS))
    for u in range(nU):
        mig[u] = binning.migration_slots(reco_s, true_s, weights=wU[u][ss])
        bkg[u] = binning.hist_slots(pt[sb], pl[sb], weights=wU[u][sb])
    out["migration"], out["bkg"] = mig, bkg
    # truth: eff_denom per universe
    denom = signal.is_efficiency_denominator(t["mc_incoming"], t["mc_current"],
                                             t["mc_vtx"], t["mc_primFSLepton"])
    ptd, pld = true_pt_pl(t["mc_primFSLepton"])
    baset = weights.truth_cv_weight(t, wts["flux"], wts["w2p2h"], wts["rpa"])
    wUt = variation_weight_columns({k: t[k] for k in TRUTH_BR}, baset, wts["w2p2h"])
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
    parser = make_parser("2p2h systematic universes (3-universe vertical pass).")
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
    wts["rpa"].weight(np.array([0.1]), np.array([0.4]), np.array([1]), np.array([6]))

    with RunLog(__file__, f"2p2h universes ({len(MODES)}), {args.playlist}",
                inputs={**args_to_inputs(args), "n_mc_files": len(mc), "modes": MODES}) as log:
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
        np.savez(outdir / "twop2h_universes.npz", migration=tot["migration"],
                 bkg=tot["bkg"], eff_denom=tot["eff_denom"],
                 modes=np.array(MODES), pot_mc=pot)
        summary = {"playlist": args.playlist, "n_universes": len(MODES),
                   "files_ok": n_ok, "files_failed": [f["url"] for f in fail],
                   "wall_s": round(time.time() - t0, 1)}
        (outdir / "summary.json").write_text(json.dumps(summary, indent=2))
        log.out("outdir", str(outdir)); log.out("summary", summary)
        print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
