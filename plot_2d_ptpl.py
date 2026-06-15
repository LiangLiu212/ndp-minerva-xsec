"""Full-playlist 2D (p_T x p_parallel) distributions + migration matrix.

Streams every AnaTuple in the given playlist file lists (xrootd, no local
copies), applies the certified reco selection to data and MC, splits MC by
the truth signal definition, and fills on the paper's 224-cell grid
(xsec/binning.py == anc/bin_mapping.txt):

  - data reco 2D (selected)
  - MC reco 2D (selected: all / signal / background)
  - migration matrix M[true_gid, reco_gid] for selected signal events
    (true kinematics from mc_primFSLepton with the NuMI beam rotation)

Migration follows the tutorial semantics: filled for reco-selected events
passing the SIGNAL DEFINITION only (no truth phase-space requirement);
pairs with either side outside the grid are counted separately, not filled.

Outputs: <outdir>/hists.npz + plots (data_2d.png, mc_2d_signal.png,
migration.png) + RunLog with POT and count summaries.
"""
import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot
from matplotlib.colors import LogNorm

from runlog_tools import (RunLog, add_label, args_to_inputs, default_outdir,
                          make_parser)
from xsec import binning, cuts, signal, weights
from xsec.flux import FluxCV
from xsec.kinematics import reco_pt_pz_gev, true_theta_p

DATA_BRANCHES = list(cuts.RECO_SELECTION_BRANCHES) + ["MasterAnaDev_leptonE"]
MC_BRANCHES_NOWGT = DATA_BRANCHES + list(signal.SIGNAL_BRANCHES) + ["mc_primFSLepton"]
# weighted runs additionally need the MnvTunev1 CV-weight inputs
MC_BRANCHES_CV = list(dict.fromkeys(MC_BRANCHES_NOWGT + list(weights.RECO_WEIGHT_BRANCHES)))


def read_file(url, role, weighters=None):
    """Stream one AnaTuple; return fills + counters for accumulation.

    weighters: None (unweighted) or {flux, w2p2h, rpa}. Data is never weighted;
    when weighting MC, fills carry the per-event reco CV weight and the run
    also tracks sum(w^2) for the selected-signal fills (Poisson-error scaling).
    """
    out = {"url": url, "role": role}
    mc_branches = MC_BRANCHES_CV if weighters else MC_BRANCHES_NOWGT
    with uproot.open(url) as f:
        meta = f["Meta"].arrays(["POT_Used"], library="np")
        out["pot_used"] = float(meta["POT_Used"].sum())
        arrs = f["MasterAnaDev"].arrays(DATA_BRANCHES if role == "data" else mc_branches,
                                        library="np")
    sel = cuts.reco_selection(arrs)
    pt, pl = reco_pt_pz_gev(arrs["MasterAnaDev_leptonE"],
                            arrs["muon_thetaX"], arrs["muon_thetaY"])
    out["n_total"] = int(sel.size)
    out["n_selected"] = int(sel.sum())

    if role == "data":
        out["h2_data"] = binning.hist2d(pt[sel], pl[sel])
        out["n_sel_in_grid"] = int(out["h2_data"].sum())
        return out

    sig = signal.is_signal(arrs["mc_incoming"], arrs["mc_current"])
    sel_sig, sel_bkg = sel & sig, sel & ~sig
    out["n_signal"] = int(sel_sig.sum())
    out["n_background"] = int(sel_bkg.sum())

    if weighters is None:
        w_sel = w_sig = w_bkg = None
        out["sumw_signal"] = float(sel_sig.sum())
    else:
        w = weights.reco_cv_weight(arrs, weighters["flux"], weighters["w2p2h"],
                                   weighters["rpa"])
        w_sel, w_sig, w_bkg = w[sel], w[sel_sig], w[sel_bkg]
        out["sumw_signal"] = float(w_sig.sum())

    out["h2_mc_selected"] = binning.hist2d(pt[sel], pl[sel], weights=w_sel)
    out["h2_mc_signal"] = binning.hist2d(pt[sel_sig], pl[sel_sig], weights=w_sig)
    out["h2_mc_bkg"] = binning.hist2d(pt[sel_bkg], pl[sel_bkg], weights=w_bkg)

    lep = np.asarray(arrs["mc_primFSLepton"], dtype=np.float64)[sel_sig]
    theta_t, p_t = true_theta_p(lep[:, 0], lep[:, 1], lep[:, 2])
    pt_true = p_t * np.sin(theta_t) / 1000.0
    pl_true = p_t * np.cos(theta_t) / 1000.0
    reco_ids = binning.cell_ids(pt[sel_sig], pl[sel_sig])
    true_ids = binning.cell_ids(pt_true, pl_true)
    out["migration"] = binning.migration_matrix(true_ids, reco_ids, weights=w_sig)
    out["mig_both_in"] = int(((true_ids >= 0) & (reco_ids >= 0)).sum())
    out["mig_true_out"] = int(((true_ids < 0) & (reco_ids >= 0)).sum())
    out["mig_reco_out"] = int(((true_ids >= 0) & (reco_ids < 0)).sum())
    out["mig_both_out"] = int(((true_ids < 0) & (reco_ids < 0)).sum())
    return out


def read_with_retry(url, role, weighters=None, retries=1):
    for attempt in range(retries + 1):
        try:
            return read_file(url, role, weighters)
        except Exception as err:
            if attempt == retries:
                return {"url": url, "role": role, "error": f"{type(err).__name__}: {err}"}
            time.sleep(2.0)


def accumulate(results, keys_2d, counter_keys):
    total = {k: None for k in keys_2d}
    counts = dict.fromkeys(counter_keys, 0)
    pot = 0.0
    for r in results:
        pot += r["pot_used"]
        for k in keys_2d:
            if k in r:
                total[k] = r[k] if total[k] is None else total[k] + r[k]
        for k in counter_keys:
            counts[k] += r.get(k, 0)
    return total, counts, pot


def plot_grid_hist(h, title, path):
    fig, ax = plt.subplots(figsize=(8, 5))
    masked = np.ma.masked_equal(h, 0.0)
    pcm = ax.pcolormesh(binning.PL_EDGES_GEV, binning.PT_EDGES_GEV, masked,
                        norm=LogNorm(), cmap="viridis")
    ax.set_xscale("log")
    ax.set_xlabel(r"$p_{\parallel,\mu}$ [GeV/c]")
    ax.set_ylabel(r"$p_{T,\mu}$ [GeV/c]")
    ax.set_title(title)
    fig.colorbar(pcm, ax=ax, label="selected events / cell")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_migration(m, path):
    row_sums = m.sum(axis=1, keepdims=True)
    norm = np.divide(m, row_sums, out=np.zeros_like(m), where=row_sums > 0)
    fig, ax = plt.subplots(figsize=(8.5, 7))
    pcm = ax.imshow(np.ma.masked_equal(norm, 0.0), origin="lower",
                    norm=LogNorm(vmin=1e-4, vmax=1.0), cmap="viridis",
                    extent=(-0.5, binning.N_CELLS - 0.5, -0.5, binning.N_CELLS - 0.5))
    for g in range(0, binning.N_CELLS + 1, binning.N_PL):
        ax.axhline(g - 0.5, color="w", lw=0.25, alpha=0.5)
        ax.axvline(g - 0.5, color="w", lw=0.25, alpha=0.5)
    ax.set_xlabel("reco cell GlobalID  (blocks of 16 = one $p_T$ bin)")
    ax.set_ylabel("true cell GlobalID")
    ax.set_title(r"Migration $P(\mathrm{reco}\,|\,\mathrm{true})$ — selected signal")
    fig.colorbar(pcm, ax=ax, label="row-normalized probability")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main():
    parser = make_parser("Full-playlist 2D pT x p|| distributions + migration "
                         "matrix on the paper grid (streaming-only).")
    parser.add_argument("--data-list",
                        default="config/playlists/MediumEnergy_FHC_Data_Playlist1A.txt")
    parser.add_argument("--mc-list",
                        default="config/playlists/MediumEnergy_FHC_StandardMC_Playlist1A.txt")
    parser.add_argument("--max-data-files", type=int, default=None,
                        help="process only the first N data files (testing)")
    parser.add_argument("--max-mc-files", type=int, default=None,
                        help="process only the first N MC files (testing)")
    parser.add_argument("--workers", type=int, default=6,
                        help="concurrent xrootd streams")
    parser.add_argument("--weights", choices=["none", "cv"], default="none",
                        help="MC weighting: 'none' (unweighted CV) or 'cv' (MnvTunev1)")
    parser.add_argument("--playlist", default="minervame1A",
                        help="playlist name (selects the flux group for --weights cv)")
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
                     "w2p2h": weights.TwoP2HWeight(),
                     "rpa": weights.RPAWeight()}
        # pre-warm the RPA carbon table (lru_cache) before threads dispatch
        weighters["rpa"].weight(np.array([0.1]), np.array([0.4]),
                                np.array([1]), np.array([6]))
    wtag = f"weighted ({args.weights})" if args.weights != "none" else "unweighted CV"

    with RunLog(__file__, f"2D pT x p|| + migration, {args.playlist}, {wtag} "
                          f"({len(data_urls)} data + {len(mc_urls)} MC files)",
                inputs={**args_to_inputs(args),
                        "n_data_files": len(data_urls),
                        "n_mc_files": len(mc_urls)}) as log:
        t0 = time.time()
        jobs = [(u, "data") for u in data_urls] + [(u, "mc") for u in mc_urls]
        results, failures = [], []
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futs = {pool.submit(read_with_retry, u, role,
                                None if role == "data" else weighters): (u, role)
                    for u, role in jobs}
            for i, fut in enumerate(as_completed(futs), 1):
                r = fut.result()
                (failures if "error" in r else results).append(r)
                if i % 25 == 0 or i == len(jobs):
                    print(f"  [{i}/{len(jobs)}] streamed, {len(failures)} failures, "
                          f"{time.time()-t0:.0f}s elapsed", flush=True)

        data_res = [r for r in results if r["role"] == "data"]
        mc_res = [r for r in results if r["role"] == "mc"]

        d_tot, d_counts, pot_data = accumulate(
            data_res, ["h2_data"], ["n_total", "n_selected", "n_sel_in_grid"])
        m_tot, m_counts, pot_mc = accumulate(
            mc_res, ["h2_mc_selected", "h2_mc_signal", "h2_mc_bkg", "migration"],
            ["n_total", "n_selected", "n_signal", "n_background", "sumw_signal",
             "mig_both_in", "mig_true_out", "mig_reco_out", "mig_both_out"])

        np.savez(outdir / "hists.npz",
                 data2d=d_tot["h2_data"], mc2d_selected=m_tot["h2_mc_selected"],
                 mc2d_signal=m_tot["h2_mc_signal"], mc2d_bkg=m_tot["h2_mc_bkg"],
                 migration=m_tot["migration"],
                 pt_edges=binning.PT_EDGES_GEV, pl_edges=binning.PL_EDGES_GEV,
                 pot_data=pot_data, pot_mc=pot_mc,
                 weight_mode=np.array(args.weights))

        plot_grid_hist(d_tot["h2_data"],
                       f"Data, {args.playlist}: {int(d_tot['h2_data'].sum())} selected "
                       f"in grid ({pot_data:.3e} POT)", outdir / "data_2d.png")
        plot_grid_hist(m_tot["h2_mc_signal"],
                       f"MC signal ({wtag}), {args.playlist}: "
                       f"{m_tot['h2_mc_signal'].sum():.0f} selected in grid "
                       f"({pot_mc:.3e} POT)", outdir / "mc_2d_signal.png")
        plot_migration(m_tot["migration"], outdir / "migration.png")

        summary = {
            "playlist": args.playlist,
            "weight_mode": args.weights,
            "files_ok": {"data": len(data_res), "mc": len(mc_res)},
            "files_failed": [f["url"] for f in failures],
            "pot": {"data": pot_data, "mc": pot_mc,
                    "ratio_data_over_mc": pot_data / pot_mc if pot_mc else None},
            "data": d_counts,
            "mc": m_counts,
            "mc_signal_weighted_in_grid": float(m_tot["h2_mc_signal"].sum()),
            "data_over_mc_signal_potscaled":
                float(d_tot["h2_data"].sum()
                      / (m_tot["h2_mc_signal"].sum() * pot_data / pot_mc))
                if pot_mc and m_tot["h2_mc_signal"].sum() else None,
            "migration_in_grid_fraction":
                m_counts["mig_both_in"] / max(1, m_counts["n_signal"]),
            "wall_s": round(time.time() - t0, 1),
        }
        (outdir / "summary.json").write_text(json.dumps(summary, indent=2))

        log.out("outdir", str(outdir))
        log.out("hists", str(outdir / "hists.npz"))
        log.out("plots", [str(outdir / p) for p in
                          ("data_2d.png", "mc_2d_signal.png", "migration.png")])
        log.out("summary", summary)

        print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
