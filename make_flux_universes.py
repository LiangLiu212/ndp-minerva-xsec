"""Flux PPFX systematic universes (vertical weight pass) — the flux *shape* term.

The flux uncertainty is the per-universe spread of the PPFX-throw fluxes,
ν-e-constraint-weighted (MnvHistoConstrainer::CorrectFluxUniv -> SetUnivWgt ->
MnvVertErrorBand::CalcCovMx weighted covariance). This upgrades the S2
normalization model (a flat 3.23 %) to the per-cell, off-diagonal-resolved
covariance.

Each universe u reweights EVERY MC fill by U_u(Enu)/Φ_constrained(Enu)
(FluxCV.universe_ratio) — vertical, so reco/true slots are the CV ones — and
additionally uses its own integrated flux Φ_u in the denominator (applied
downstream in assemble_flux.py, from FluxCV.universe_integrals).

N = 100 universes by default (the first 100 PPFX throws; matches the CCQENu
`fluxUniverses: 100` config). Memory-light: the per-universe ratio is computed
and discarded inside the fill loop, so the full (N, n_events) weight matrix is
never built (the Truth tree has ~544k events).
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


def read_mc(url, wts, n_univ):
    out = {"url": url}
    with uproot.open(url) as f:
        out["pot_used"] = float(f["Meta"]["POT_Used"].array(library="np").sum())
        a = f["MasterAnaDev"].arrays(RECO_BR, library="np")
        t = f["Truth"].arrays(TRUTH_BR, library="np")
    flux = wts["flux"]
    # reco: migration + bkg per universe (slots fixed; only the flux weight varies)
    sel = cuts.reco_selection(a)
    sig = signal.is_signal(a["mc_incoming"], a["mc_current"])
    ss, sb = sel & sig, sel & ~sig
    pt, pl = reco_pt_pz_gev(a["MasterAnaDev_leptonE"], a["muon_thetaX"], a["muon_thetaY"])
    ptt, plt = true_pt_pl(a["mc_primFSLepton"])
    base = weights.reco_cv_weight(a, flux, wts["w2p2h"], wts["rpa"])
    enu = np.asarray(a["mc_incomingE"], np.float64) / 1000.0
    reco_s = binning.slot_ids(pt[ss], pl[ss]); true_s = binning.slot_ids(ptt[ss], plt[ss])
    pt_b, pl_b = pt[sb], pl[sb]
    mig = np.zeros((n_univ, binning.N_SLOTS, binning.N_SLOTS))
    bkg = np.zeros((n_univ, binning.N_SLOTS))
    for u in range(n_univ):
        wu = base * flux.universe_ratio(u, enu)
        mig[u] = binning.migration_slots(reco_s, true_s, weights=wu[ss])
        bkg[u] = binning.hist_slots(pt_b, pl_b, weights=wu[sb])
    out["migration"], out["bkg"] = mig, bkg
    # truth: eff_denom per universe (flux weight also reshapes the denominator)
    denom = signal.is_efficiency_denominator(t["mc_incoming"], t["mc_current"],
                                             t["mc_vtx"], t["mc_primFSLepton"])
    ptd, pld = true_pt_pl(t["mc_primFSLepton"])
    baset = weights.truth_cv_weight(t, flux, wts["w2p2h"], wts["rpa"])
    enut = np.asarray(t["mc_incomingE"], np.float64) / 1000.0
    dsl = binning.slot_ids(ptd[denom], pld[denom])
    edges = np.arange(binning.N_SLOTS + 1) - 0.5
    den = np.zeros((n_univ, binning.N_SLOTS))
    for u in range(n_univ):
        wut = baset * flux.universe_ratio(u, enut)
        den[u], _ = np.histogram(dsl, bins=edges, weights=wut[denom])
    out["eff_denom"] = den
    return out


def with_retry(url, wts, n_univ, retries=1):
    for attempt in range(retries + 1):
        try:
            return read_mc(url, wts, n_univ)
        except Exception as err:
            if attempt == retries:
                return {"url": url, "error": f"{type(err).__name__}: {err}"}
            time.sleep(2.0)


def main():
    parser = make_parser("Flux PPFX systematic universes (vertical weight pass).")
    parser.add_argument("--mc-list",
                        default="config/playlists/MediumEnergy_FHC_StandardMC_Playlist1A.txt")
    parser.add_argument("--playlist", default="minervame1A")
    parser.add_argument("--n-universes", type=int, default=100,
                        help="first N PPFX universes (CCQENu uses 100)")
    parser.add_argument("--max-mc-files", type=int, default=None)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--outdir", default=None)
    add_label(parser)
    args = parser.parse_args()

    mc = [u.strip() for u in Path(args.mc_list).read_text().splitlines() if u.strip()]
    if args.max_mc_files:
        mc = mc[:args.max_mc_files]
    outdir = Path(args.outdir) if args.outdir else default_outdir(__file__)
    outdir.mkdir(parents=True, exist_ok=True)
    flux = FluxCV(args.playlist)
    n_univ = min(args.n_universes, flux.n_universes)
    wts = {"flux": flux, "w2p2h": weights.TwoP2HWeight(), "rpa": weights.RPAWeight()}
    wts["rpa"].weight(np.array([0.1]), np.array([0.4]), np.array([1]), np.array([6]))

    with RunLog(__file__, f"flux PPFX universes ({n_univ}), {args.playlist}",
                inputs={**args_to_inputs(args), "n_mc_files": len(mc),
                        "n_universes": n_univ}) as log:
        t0 = time.time(); fail = []; n_ok = 0; pot = 0.0
        tot = {"migration": None, "bkg": None, "eff_denom": None}
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futs = {pool.submit(with_retry, u, wts, n_univ): u for u in mc}
            for i, fut in enumerate(as_completed(futs), 1):
                r = fut.result()
                if "error" in r:
                    fail.append(r)
                else:
                    n_ok += 1; pot += r["pot_used"]
                    for k in tot:
                        tot[k] = r[k].copy() if tot[k] is None else tot[k] + r[k]
                    r.clear()
                if i % 5 == 0 or i == len(mc):
                    print(f"  [{i}/{len(mc)}] {len(fail)} fail, {time.time()-t0:.0f}s", flush=True)
        np.savez(outdir / "flux_universes.npz", migration=tot["migration"],
                 bkg=tot["bkg"], eff_denom=tot["eff_denom"],
                 n_universes=n_univ, pot_mc=pot)
        summary = {"playlist": args.playlist, "n_universes": n_univ,
                   "files_ok": n_ok, "files_failed": [f["url"] for f in fail],
                   "wall_s": round(time.time() - t0, 1)}
        (outdir / "summary.json").write_text(json.dumps(summary, indent=2))
        log.out("outdir", str(outdir)); log.out("summary", summary)
        print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
