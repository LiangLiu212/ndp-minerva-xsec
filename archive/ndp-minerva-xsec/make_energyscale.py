"""S4 — muon energy-scale (lateral) shifted ingredients.

A reco muon-momentum scale shift p -> p(1±δ) scales reco p_T and p_∥ by (1±δ),
which re-bins events (the 6-cut selection is |p|-independent, so the selected
SET, the truth, and the efficiency num/denom are unchanged — only the reco
distribution of the migration and background moves). So this streams the MC
reco tree once and produces the +δ and -δ shifted migration + background; the
CV eff_num/eff_denom/data come from the existing ingredients.npz.

δ defaults to the MINOS-range fractional uncertainty (NSFDefaults
MinosMuonPRange_Err = 0.984 %), the dominant component for these MINOS-matched
forward muons. The full per-event model (MINERvA absolute dE/dx + material, and
MINOS curvature terms) is a momentum-dependent refinement; this flat-fractional
model captures the leading ~1 % scale.

Cov_energyscale = pair_covariance(σ⁺, σ⁻) is assembled downstream.
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

MINOS_MUON_SCALE = 0.00984   # NSFDefaults::MinosMuonPRange_Err
MC_BRANCHES = list(dict.fromkeys(
    list(cuts.RECO_SELECTION_BRANCHES) + ["MasterAnaDev_leptonE"]
    + list(signal.SIGNAL_BRANCHES) + list(weights.RECO_WEIGHT_BRANCHES)))


def true_pt_pl(lep):
    lep = np.asarray(lep, dtype=np.float64)
    th, p = true_theta_p(lep[:, 0], lep[:, 1], lep[:, 2])
    return p * np.sin(th) / 1000.0, p * np.cos(th) / 1000.0


def read_mc(url, weighters, delta):
    out = {"url": url}
    with uproot.open(url) as f:
        out["pot_used"] = float(f["Meta"]["POT_Used"].array(library="np").sum())
        a = f["MasterAnaDev"].arrays(MC_BRANCHES, library="np")
    sel = cuts.reco_selection(a)
    sig = signal.is_signal(a["mc_incoming"], a["mc_current"])
    sel_sig, sel_bkg = sel & sig, sel & ~sig
    pt, pl = reco_pt_pz_gev(a["MasterAnaDev_leptonE"], a["muon_thetaX"], a["muon_thetaY"])
    ptt, plt = true_pt_pl(a["mc_primFSLepton"])
    w = weights.reco_cv_weight(a, weighters["flux"], weighters["w2p2h"], weighters["rpa"])
    true_sig = binning.slot_ids(ptt[sel_sig], plt[sel_sig])
    for tag, s in (("plus", 1 + delta), ("minus", 1 - delta)):
        out[f"migration_{tag}"] = binning.migration_slots(
            binning.slot_ids(pt[sel_sig] * s, pl[sel_sig] * s), true_sig, weights=w[sel_sig])
        out[f"bkg_{tag}"] = binning.hist_slots(pt[sel_bkg] * s, pl[sel_bkg] * s, weights=w[sel_bkg])
    return out


def with_retry(url, weighters, delta, retries=1):
    for attempt in range(retries + 1):
        try:
            return read_mc(url, weighters, delta)
        except Exception as err:
            if attempt == retries:
                return {"url": url, "error": f"{type(err).__name__}: {err}"}
            time.sleep(2.0)


def main():
    parser = make_parser("Muon energy-scale shifted migration+bkg (reco-only stream).")
    parser.add_argument("--mc-list",
                        default="config/playlists/MediumEnergy_FHC_StandardMC_Playlist1A.txt")
    parser.add_argument("--playlist", default="minervame1A")
    parser.add_argument("--delta", type=float, default=MINOS_MUON_SCALE)
    parser.add_argument("--max-mc-files", type=int, default=None)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--outdir", default=None)
    add_label(parser)
    args = parser.parse_args()

    mc_urls = [u.strip() for u in Path(args.mc_list).read_text().splitlines() if u.strip()]
    if args.max_mc_files:
        mc_urls = mc_urls[:args.max_mc_files]
    outdir = Path(args.outdir) if args.outdir else default_outdir(__file__)
    outdir.mkdir(parents=True, exist_ok=True)

    weighters = {"flux": FluxCV(args.playlist), "w2p2h": weights.TwoP2HWeight(),
                 "rpa": weights.RPAWeight()}
    weighters["rpa"].weight(np.array([0.1]), np.array([0.4]), np.array([1]), np.array([6]))

    with RunLog(__file__, f"energy-scale shift δ={args.delta}, {args.playlist}",
                inputs={**args_to_inputs(args), "n_mc_files": len(mc_urls)}) as log:
        t0 = time.time()
        results, failures = [], []
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futs = {pool.submit(with_retry, u, weighters, args.delta): u for u in mc_urls}
            for i, fut in enumerate(as_completed(futs), 1):
                r = fut.result()
                (failures if "error" in r else results).append(r)
                if i % 10 == 0 or i == len(mc_urls):
                    print(f"  [{i}/{len(mc_urls)}] {len(failures)} fail, {time.time()-t0:.0f}s", flush=True)

        acc = lambda k: sum(r[k] for r in results if k in r)
        np.savez(outdir / "energyscale.npz",
                 migration_plus=acc("migration_plus"), migration_minus=acc("migration_minus"),
                 bkg_plus=acc("bkg_plus"), bkg_minus=acc("bkg_minus"),
                 delta=args.delta, pot_mc=sum(r["pot_used"] for r in results))
        summary = {"playlist": args.playlist, "delta": args.delta,
                   "files_ok": len(results), "files_failed": [f["url"] for f in failures],
                   "wall_s": round(time.time() - t0, 1)}
        (outdir / "summary.json").write_text(json.dumps(summary, indent=2))
        log.out("outdir", str(outdir))
        log.out("energyscale", str(outdir / "energyscale.npz"))
        log.out("summary", summary)
        print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
