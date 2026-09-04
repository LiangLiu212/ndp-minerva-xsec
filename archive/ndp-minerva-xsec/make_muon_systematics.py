"""Muon-reconstruction systematics — all bands in one reco-only streaming pass.

The dominant systematic category (Fig 8 "Muon Reconstruction"). Every band is a
reco-side shift, so the TRUTH tree (efficiency denominator) is untouched — only
the reco migration + background move (and, for the vertical MINOS-efficiency
band, get reweighted). The CV eff_denom from ingredients.npz is reused in the
assembly. So this streams only the MasterAnaDev reco tree (no Truth, no stall).

Bands (each a ±1σ pair → its own pair covariance; source in xsec.muon_syst):
  momentum-scale (re-bin only, selection is angle-only):
    Muon_Energy_MINERvA, Muon_Energy_MINOS, Muon_Energy_Resolution
  track-angle (re-bin AND re-apply the 20° muon-angle cut):
    BeamAngleX, BeamAngleY, MuonAngleXResolution, MuonAngleYResolution
  vertical weight (slots = CV, reweight):
    MINOS_Efficiency

cov_energyscale = MINERvA ⊕ MINOS (validated vs the anc file); the full muon-
reco group = all eight bands, assembled by assemble_muon.py.
"""
import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import uproot

from runlog_tools import (RunLog, add_label, args_to_inputs, default_outdir,
                          make_parser)
from xsec import binning, cuts, kinematics, muon_syst, signal, weights
from xsec.flux import FluxCV
from xsec.kinematics import reco_pt_pz_gev, true_theta_p

MOM_BANDS = ["Muon_Energy_MINERvA", "Muon_Energy_MINOS", "Muon_Energy_Resolution"]
ANGLE_BANDS = ["BeamAngleX", "BeamAngleY", "MuonAngleXResolution", "MuonAngleYResolution"]
VERT_BANDS = ["MINOS_Efficiency"]
BANDS = MOM_BANDS + ANGLE_BANDS + VERT_BANDS

RECO_BR = list(dict.fromkeys(
    list(cuts.RECO_SELECTION_BRANCHES) + ["MasterAnaDev_leptonE",
    "MasterAnaDev_minos_trk_p", "MasterAnaDev_minos_used_curvature"]
    + list(signal.SIGNAL_BRANCHES) + list(weights.RECO_WEIGHT_BRANCHES)))


def true_pt_pl(lep):
    lep = np.asarray(lep, np.float64)
    th, p = true_theta_p(lep[:, 0], lep[:, 1], lep[:, 2])
    return p * np.sin(th) / 1000.0, p * np.cos(th) / 1000.0


def read_mc(url, wts):
    out = {"url": url}
    with uproot.open(url) as f:
        out["pot_used"] = float(f["Meta"]["POT_Used"].array(library="np").sum())
        a = f["MasterAnaDev"].arrays(RECO_BR, library="np")

    sel_cv = cuts.reco_selection(a)
    sig = signal.is_signal(a["mc_incoming"], a["mc_current"])
    lep = a["MasterAnaDev_leptonE"]
    tx0, ty0 = a["muon_thetaX"], a["muon_thetaY"]
    p_total = np.linalg.norm(np.asarray(lep, np.float64)[:, :3], axis=1)
    p_minos = a["MasterAnaDev_minos_trk_p"]
    used_curv = a["MasterAnaDev_minos_used_curvature"]
    p_true = np.linalg.norm(np.asarray(a["mc_primFSLepton"], np.float64)[:, :3], axis=1)
    true_pt, true_pl = true_pt_pl(a["mc_primFSLepton"])
    truex, truey = muon_syst.true_theta_xy(a["mc_primFSLepton"])
    pt_cv, pl_cv = reco_pt_pz_gev(lep, tx0, ty0)
    base = weights.reco_cv_weight(a, wts["flux"], wts["w2p2h"], wts["rpa"])
    theta_deg = np.degrees(kinematics.theta3d(tx0, ty0))
    eff_args = (a["MasterAnaDev_minos_trk_p"], a["numi_pot"],
                a["batch_structure"], a["reco_vertex_batch"])

    def fill(sel, pt, pl, w):
        ss, sb = sel & sig, sel & ~sig
        mig = binning.migration_slots(binning.slot_ids(pt[ss], pl[ss]),
                                      binning.slot_ids(true_pt[ss], true_pl[ss]),
                                      weights=w[ss])
        bkg = binning.hist_slots(pt[sb], pl[sb], weights=w[sb])
        return mig, bkg

    def reselect(tx, ty):
        a2 = dict(a); a2["muon_thetaX"] = tx; a2["muon_thetaY"] = ty
        return cuts.reco_selection(a2)

    res = {}
    for band in BANDS:
        for tag, ns in (("plus", 1.0), ("minus", -1.0)):
            if band == "Muon_Energy_MINERvA":
                sc = muon_syst.pmu_minerva_scale(p_total, ns)
                sel, pt, pl, w = sel_cv, pt_cv * sc, pl_cv * sc, base
            elif band == "Muon_Energy_MINOS":
                sc = muon_syst.pmu_minos_scale(p_total, p_minos, used_curv, ns)
                sel, pt, pl, w = sel_cv, pt_cv * sc, pl_cv * sc, base
            elif band == "Muon_Energy_Resolution":
                sc = muon_syst.pmu_resolution_scale(p_total, p_true, ns)
                sel, pt, pl, w = sel_cv, pt_cv * sc, pl_cv * sc, base
            elif band in ("BeamAngleX", "BeamAngleY"):
                ax = "x" if band.endswith("X") else "y"
                tx, ty = muon_syst.beam_angle_shift(tx0, ty0, ax, ns)
                pt, pl = reco_pt_pz_gev(lep, tx, ty); sel, w = reselect(tx, ty), base
            elif band in ("MuonAngleXResolution", "MuonAngleYResolution"):
                ax = "x" if "X" in band else "y"
                tx, ty = muon_syst.angle_resolution_shift(tx0, ty0, truex, truey, ax, ns)
                pt, pl = reco_pt_pz_gev(lep, tx, ty); sel, w = reselect(tx, ty), base
            else:                                       # MINOS_Efficiency (vertical)
                ratio = weights.minos_efficiency_ratio(*eff_args, theta_deg, ns)
                sel, pt, pl, w = sel_cv, pt_cv, pl_cv, base * ratio
            mig, bkg = fill(sel, pt, pl, w)
            res[f"{band}__mig__{tag}"] = mig
            res[f"{band}__bkg__{tag}"] = bkg
    out["res"] = res
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
    parser = make_parser("Muon-reconstruction systematics (all bands, reco-only stream).")
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

    with RunLog(__file__, f"muon-reco systematics ({len(BANDS)} bands), {args.playlist}",
                inputs={**args_to_inputs(args), "n_mc_files": len(mc), "bands": BANDS}) as log:
        t0 = time.time(); fail = []; n_ok = 0; pot = 0.0; tot = None
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futs = {pool.submit(with_retry, u, wts): u for u in mc}
            for i, fut in enumerate(as_completed(futs), 1):
                r = fut.result()
                if "error" in r:
                    fail.append(r)
                else:                                   # incremental accumulation
                    n_ok += 1; pot += r["pot_used"]
                    if tot is None:
                        tot = {k: v.copy() for k, v in r["res"].items()}
                    else:
                        for k, v in r["res"].items():
                            tot[k] += v
                    r.clear()
                if i % 5 == 0 or i == len(mc):
                    print(f"  [{i}/{len(mc)}] {len(fail)} fail, {time.time()-t0:.0f}s", flush=True)
        np.savez(outdir / "muon_universes.npz", bands=np.array(BANDS),
                 pot_mc=pot, **tot)
        summary = {"playlist": args.playlist, "n_bands": len(BANDS),
                   "files_ok": n_ok, "files_failed": [f["url"] for f in fail],
                   "wall_s": round(time.time() - t0, 1)}
        (outdir / "summary.json").write_text(json.dumps(summary, indent=2))
        log.out("outdir", str(outdir)); log.out("summary", summary)
        print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
