"""Diagnostics runs that are not model tests: apply a channel's truth-level signal definition
to the cached official MC and record what it selects, in a manifest-backed run directory.

    runs/<date>_signal_<channel>[_N]/
      manifest.json     inputs (cache fingerprint), channel snapshot hash, git state, versions
      channel.json      the channel manifest snapshot
      cutflow.json      cumulative counts per step, for all cached truth events and inside the fiducial vertex
      summary.json      composition by interaction type, proton multiplicities, veto inventory, parity vs MAT
      report.md         the same, readable
      figs/             muon and leading-proton kinematics of the selected signal
"""
from __future__ import annotations

from pathlib import Path
import json

import numpy as np

from .config import SiteConfig, load_site_config
from .channels import load_channel, ChannelSpec
from .channels import observables as obs
from .channels import signal as sig
from .events import TruthTable, INT_TYPES
from .io import dump_json, cheap_fingerprint, git_state, versions, timestamp, unique_run_dir, ensure_dir, sha256_text
from . import __version__

#: MAT-MINERvA `IsQELike` (universes/CCQE3DFitsSystematics.cxx) enumerated PDG lists, for parity.
_MAT_MESONS = (211, 321, 323, 111, 130, 310, 311, 313)            # abs() for the first three
_MAT_HEAVY = (3112, 3122, 3212, 3222, 4112, 4122, 4212, 4222, 411, 421, 111)


def mat_isqelike(t: TruthTable, photon_E_max_gev: float = 0.010) -> np.ndarray:
    """Transcription of MAT's IsQELike on a TruthTable (one |mu| in the FS list, no listed mesons /
    heavy baryons, no photons above threshold). No kinematic window — vetoes only."""
    pdg, E = t["fs_pdg"], t["fs_E"]
    apdg = np.abs(pdg)
    n_mu = sig.fs_count(t, apdg == 13)
    n_ph = sig.fs_count(t, (pdg == 22) & (E > photon_E_max_gev))
    mes = np.isin(apdg, _MAT_MESONS[:3]) | np.isin(pdg, _MAT_MESONS[3:])
    n_mes = sig.fs_count(t, mes & ~((pdg == 22) & (E > photon_E_max_gev)))
    n_heavy = sig.fs_count(t, np.isin(pdg, _MAT_HEAVY) & ~mes)
    return (n_mu == 1) & (n_mes == 0) & (n_heavy == 0) & (n_ph == 0)


def _hist(ax, x, edges, label, xlabel):
    ax.hist(x[np.isfinite(x)], bins=edges, histtype="step", lw=1.4, label=label)
    ax.set_xlabel(xlabel); ax.set_ylabel("events (unweighted)")


def run_signal_diagnostics(channel_name: str, cfg: SiteConfig | None = None, *, cache: str | Path | None = None,
                           out_root: str | Path | None = None, slug: str | None = None) -> Path:
    cfg = cfg or load_site_config()
    ch: ChannelSpec = load_channel(channel_name)
    if cache is None:
        from .adapters.minerva_anatuple import cache_tag
        mc = ch.data.get("reco_mc_files", [])
        if not mc or cfg.data_dir is None:
            raise FileNotFoundError("channel names no reco_mc_files or site has no data_dir; pass cache=")
        cache = cfg.data_dir / "cache" / f"truth_{cache_tag(mc[0])}.npz"
    cache = Path(cache)
    t = TruthTable.load(cache)
    frame = ch.frame

    # ---- cutflow, with and without the fiducial vertex --------------------------------------
    steps = ch.signal_cutflow(t)
    fid = ch.in_phase_space(t)                     # for this kind of channel: the fiducial vertex only
    w = t["weight"]
    rows = [{"step": "all_cached_truth_events", "n": int(t.n), "sumw": float(w.sum()),
             "n_fiducial": int(fid.sum()), "sumw_fiducial": float(w[fid].sum())}]
    for name, m in steps:
        rows.append({"step": name, "n": int(m.sum()), "sumw": float(w[m].sum()),
                     "n_fiducial": int((m & fid).sum()), "sumw_fiducial": float(w[m & fid].sum())})
    signal = steps[-1][1]
    sel = signal & fid

    # ---- composition + leading proton ----------------------------------------------------------
    it = t["int_type"]
    comp = {INT_TYPES.get(int(k), str(k)): {"n": int((sel & (it == k)).sum()),
                                           "fraction": float((sel & (it == k)).sum() / max(sel.sum(), 1))}
            for k in np.unique(it[sel])}
    lp = sig.leading_proton(t, frame, ch.signal.get("proton", {}))
    cls = sig.fs_classes(t)
    n_all_protons = sig.fs_count(t, cls["proton"])
    mult_window = {int(k): int(v) for k, v in zip(*np.unique(lp["n_in_window"][sel], return_counts=True))}
    mult_all = {int(k): int(v) for k, v in zip(*np.unique(n_all_protons[sel], return_counts=True))}

    # ---- veto inventory on the kinematically accepted CC events ------------------------------
    kin = dict(steps)["proton_in_window"] & fid
    pdg = t["fs_pdg"]
    ev_kin = kin[t.fs_event_index()]
    inventory = {}
    for cname in ("meson", "heavy_baryon", "photon", "charged_lepton", "neutrino", "nucleus", "pseudo"):
        m = cls[cname] & ev_kin
        codes, counts = np.unique(pdg[m], return_counts=True)
        inventory[cname] = {"n_events_with_any": int((sig.fs_count(t, cls[cname]) > 0)[kin].sum()),
                            "pdg_counts": {int(c): int(n) for c, n in zip(codes, counts)}}
    hard = cls["photon"] & (t["fs_E"] > float(ch.signal.get("veto", {}).get("photon_E_max_gev", 0.010)))
    inventory["photon"]["n_events_with_photon_above_threshold"] = int((sig.fs_count(t, hard) > 0)[kin].sum())
    anti = np.isin(pdg, [-2212, -2112]) & ev_kin
    inventory["antinucleons_not_vetoed"] = {"n_particles": int(anti.sum()),
                                            "n_events": int((sig.fs_count(t, np.isin(pdg, [-2212, -2112])) > 0)[kin].sum())}
    n_fs_mu = sig.fs_count(t, np.abs(pdg) == 13)
    extra_mu = {"n_events_with_fs_muon_count_ne_1": int((n_fs_mu[kin] != 1).sum())}

    # ---- parity of the veto classes vs MAT's enumerated IsQELike -----------------------------
    ours_veto = (sig.fs_count(t, cls["meson"]) == 0) & (sig.fs_count(t, cls["heavy_baryon"]) == 0) & (sig.fs_count(t, hard) == 0)
    mat = mat_isqelike(t)
    numu_cc = steps[0][1]
    parity = {"scope": "numu CC events with a mu- primary lepton, vetoes only (no kinematic window)",
              "n_scope": int(numu_cc.sum()),
              "n_pass_ours": int((ours_veto & numu_cc).sum()), "n_pass_mat": int((mat & numu_cc).sum()),
              "n_disagree": int(((ours_veto != mat) & numu_cc).sum())}
    if parity["n_disagree"]:
        dis = (ours_veto != mat) & numu_cc
        codes, counts = np.unique(pdg[dis[t.fs_event_index()]], return_counts=True)
        parity["pdg_in_disagreeing_events"] = {int(c): int(n) for c, n in zip(codes, counts)}

    # ---- kinematics of the selected signal ----------------------------------------------------
    mu_p, mu_th = obs.lep_p(t, frame)[sel], obs.lep_theta_deg(t, frame)[sel]
    pr_p, pr_th, pr_pt = lp["p"][sel], np.rad2deg(lp["theta"][sel]), lp["pT"][sel]
    kin_summary = {
        "muon_p_gev": {"median": float(np.median(mu_p)), "p16": float(np.percentile(mu_p, 16)), "p84": float(np.percentile(mu_p, 84))},
        "muon_theta_deg": {"median": float(np.median(mu_th)), "p16": float(np.percentile(mu_th, 16)), "p84": float(np.percentile(mu_th, 84))},
        "proton_p_gev": {"median": float(np.median(pr_p)), "p16": float(np.percentile(pr_p, 16)), "p84": float(np.percentile(pr_p, 84))},
        "proton_theta_deg": {"median": float(np.median(pr_th)), "p16": float(np.percentile(pr_th, 16)), "p84": float(np.percentile(pr_th, 84))},
        "proton_pT_gev": {"median": float(np.median(pr_pt)), "p16": float(np.percentile(pr_pt, 16)), "p84": float(np.percentile(pr_pt, 84))},
    }

    # ---- run directory ------------------------------------------------------------------------
    run_dir = unique_run_dir(out_root or cfg.runs, slug or f"signal_{ch.name}")
    figs = ensure_dir(run_dir / "figs")
    dump_json(ch.to_dict(), run_dir / "channel.json")
    dump_json({"channel": ch.name, "frame": frame, "cache": str(cache), "rows": rows}, run_dir / "cutflow.json")
    summary = {"channel": ch.name, "signal": ch.signal, "frame": frame, "cache": str(cache),
               "n_signal": int(signal.sum()), "n_signal_fiducial": int(sel.sum()),
               "signal_fraction_of_numu_cc_fiducial": float(sel.sum() / max((numu_cc & fid).sum(), 1)),
               "composition_int_type_fiducial": comp,
               "protons_in_window_multiplicity_fiducial": mult_window, "all_protons_multiplicity_fiducial": mult_all,
               "veto_inventory_on_kinematically_accepted": inventory, "fs_muon_count_check": extra_mu,
               "mat_isqelike_parity": parity, "kinematics_fiducial": kin_summary}
    dump_json(summary, run_dir / "summary.json")

    # ---- TKI observables of the selected signal (only if the channel defines the leading proton) ---
    tki_summary, tki_vals = {}, {}
    if "proton" in ch.signal:
        tki_spec = [("dpT", "δp_T [GeV/c]", np.linspace(0, 2.0, 41)), ("dpTx", "δp_Tx [GeV/c]", np.linspace(-1.5, 1.5, 41)),
                    ("dpTy", "δp_Ty [GeV/c]", np.linspace(-2.5, 1.0, 41)), ("dalphaT_deg", "δα_T [deg]", np.linspace(0, 180, 37)),
                    ("dphiT_deg", "φ_T [deg]", np.linspace(0, 180, 37)), ("dpL", "δp_L [GeV/c]", np.linspace(-1.0, 1.0, 41)),
                    ("pn", "p_n [GeV/c]", np.linspace(0, 2.0, 41))]
        for name, _, _ in tki_spec:
            try:
                v = ch.evaluate(name, t)[sel]
            except KeyError as e:          # dpL / pn without observable_params.tki
                tki_summary[name] = {"error": str(e)}
                continue
            tki_vals[name] = v
            f = v[np.isfinite(v)]
            tki_summary[name] = {"n_finite": int(f.size), "median": float(np.median(f)), "p16": float(np.percentile(f, 16)),
                                 "p84": float(np.percentile(f, 84)), "mean": float(f.mean())}
        summary["tki_fiducial_signal"] = tki_summary
        dump_json(summary, run_dir / "summary.json")

    fig_paths = []
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        if tki_vals:
            fig, axs = plt.subplots(2, 4, figsize=(16, 7.5))
            for ax, (name, xlabel, edges) in zip(axs.flat, tki_spec):
                if name in tki_vals:
                    _hist(ax, tki_vals[name], edges, "signal", xlabel)
                    ax.grid(alpha=0.3)
            axs.flat[-1].axis("off")
            fig.suptitle(f"{ch.name}: TKI of the truth-level signal, {cache.name} ({frame} frame), fiducial vertex")
            fig.tight_layout()
            p = figs / "signal_tki.png"
            fig.savefig(p, dpi=130); plt.close(fig)
            fig_paths.append(p)
        fig, axs = plt.subplots(2, 3, figsize=(13, 7.5))
        _hist(axs[0, 0], mu_p, np.linspace(0, 25, 51), "signal", "muon p [GeV/c]")
        _hist(axs[0, 1], mu_th, np.linspace(0, 20, 41), "signal", "muon theta [deg]")
        _hist(axs[0, 2], n_all_protons[sel].astype(float), np.arange(-0.5, 10.5, 1), "all protons", "final-state protons per event")
        _hist(axs[0, 2], lp["n_in_window"][sel].astype(float), np.arange(-0.5, 10.5, 1), "in window", "final-state protons per event")
        axs[0, 2].legend()
        _hist(axs[1, 0], pr_p, np.linspace(0.4, 1.2, 41), "leading proton", "leading proton p [GeV/c]")
        _hist(axs[1, 1], pr_th, np.linspace(0, 80, 41), "leading proton", "leading proton theta [deg]")
        _hist(axs[1, 2], pr_pt, np.linspace(0, 1.2, 41), "leading proton", "leading proton pT [GeV/c]")
        for ax in axs.flat:
            ax.grid(alpha=0.3)
        fig.suptitle(f"{ch.name}: truth-level signal on {cache.name} ({frame} frame), fiducial vertex")
        fig.tight_layout()
        p = figs / "signal_kinematics.png"
        fig.savefig(p, dpi=130); plt.close(fig)
        fig_paths.append(p)
    except ImportError:
        pass

    lines = [f"# Truth-level signal diagnostics — `{ch.name}`", "",
             f"Cache `{cache}` ({t.n} truth entries), frame `{frame}`, signal type `{sig.signal_type(ch.signal)}`.", "",
             "## Cutflow (cumulative)", "", "| step | all | fiducial vertex |", "|---|---|---|"]
    lines += [f"| {r['step']} | {r['n']} | {r['n_fiducial']} |" for r in rows]
    lines += ["", f"Signal fraction of fiducial νμ CC (μ⁻) events: {summary['signal_fraction_of_numu_cc_fiducial']:.4f}", "",
              "## Composition of the fiducial signal by interaction type", "", "| type | n | fraction |", "|---|---|---|"]
    lines += [f"| {k} | {v['n']} | {v['fraction']:.4f} |" for k, v in comp.items()]
    lines += ["", "## Proton multiplicity (fiducial signal)", "",
              f"in window: {mult_window}", "", f"all final-state protons: {mult_all}", "",
              "## Veto inventory on kinematically accepted events (μ window + ≥1 proton in window, fiducial)", ""]
    for k, v in inventory.items():
        lines.append(f"- {k}: {json.dumps(v)}")
    lines += [f"- {json.dumps(extra_mu)}", "", "## Parity of the veto classes vs MAT `IsQELike`", "", f"{json.dumps(parity)}", "",
              "## Kinematics of the fiducial signal (median, 16th–84th percentile)", ""]
    lines += [f"- {k}: {v['median']:.3f} ({v['p16']:.3f}–{v['p84']:.3f})" for k, v in kin_summary.items()]
    if tki_summary:
        lines += ["", "## TKI observables of the fiducial signal (median, 16th–84th percentile; n finite)", ""]
        lines += [f"- {k}: {v['median']:.3f} ({v['p16']:.3f}–{v['p84']:.3f}); n = {v['n_finite']}" if "median" in v
                  else f"- {k}: {v['error']}" for k, v in tki_summary.items()]
    if fig_paths:
        lines += ["", "## Figures", ""] + [f"![{p.stem}]({p.relative_to(run_dir)})" for p in fig_paths]
    (run_dir / "report.md").write_text("\n".join(lines) + "\n")

    manifest = {"run_id": run_dir.name, "kind": "signal_diagnostics", "timestamp": timestamp(),
                "platform_version": __version__, "channel": ch.name,
                "channel_file": str(ch.path), "channel_sha256": sha256_text(Path(ch.path).read_text()) if ch.path else None,
                "inputs": {"truth_cache": cheap_fingerprint(cache), "truth_meta": t.meta},
                "git": git_state(cfg.repo_root), "versions": versions(),
                "outputs": ["channel.json", "cutflow.json", "summary.json", "report.md"] + [str(p.relative_to(run_dir)) for p in fig_paths],
                "results": {"n_signal": summary["n_signal"], "n_signal_fiducial": summary["n_signal_fiducial"],
                            "mat_parity_disagreements": parity["n_disagree"]}}
    dump_json(manifest, run_dir / "manifest.json")
    return run_dir


# =============================================================================================
# Reconstruction-level selection: cutflow, purity / efficiency, data vs official MC
# =============================================================================================
_CATEGORIES = ("signal QE", "signal 2p2h", "signal RES", "signal DIS/other",
               "bkg 1 pi+-", "bkg 1 pi0", "bkg multi-pi", "bkg other (no pion)")
_CAT_COLORS = {"signal QE": "#1f77b4", "signal 2p2h": "#17becf", "signal RES": "#2ca02c", "signal DIS/other": "#98df8a",
               "bkg 1 pi+-": "#ff7f0e", "bkg 1 pi0": "#ffbb78", "bkg multi-pi": "#d62728", "bkg other (no pion)": "#9e9e9e"}


def mc_categories(ch: ChannelSpec, rt: TruthTable) -> np.ndarray:
    """Per-candidate category label (see _CATEGORIES) from the reco rows' truth: the channel's signal
    definition + fiducial phase space split by GENIE process; background split by final-state pions."""
    sig = ch.is_signal(rt) & ch.in_phase_space(rt)
    it = rt["int_type"]
    cat = np.full(rt.n, "bkg other (no pion)", dtype=object)
    cat[sig & (it == 1)] = "signal QE"; cat[sig & (it == 5)] = "signal 2p2h"; cat[sig & (it == 2)] = "signal RES"
    cat[sig & ~np.isin(it, [1, 2, 5])] = "signal DIS/other"
    pdg = rt["fs_pdg"]
    n_pic = sig_mod_count(rt, np.abs(pdg) == 211); n_pi0 = sig_mod_count(rt, pdg == 111)
    cat[~sig & (n_pic == 1) & (n_pi0 == 0)] = "bkg 1 pi+-"
    cat[~sig & (n_pic == 0) & (n_pi0 == 1)] = "bkg 1 pi0"
    cat[~sig & (n_pic + n_pi0 >= 2)] = "bkg multi-pi"
    return cat


def sig_mod_count(rt: TruthTable, particle_mask: np.ndarray) -> np.ndarray:
    return sig.fs_count(rt, particle_mask)


def _load_reco_side(ch: ChannelSpec, cfg: SiteConfig):
    from .adapters.minerva_anatuple import cache_tag, load_reco_cache, read_pot
    data_dir = cfg.require("data_dir"); cache = data_dir / "cache"
    def concat(files, is_mc):
        recos, truths, pot, srcs = [], [], 0.0, []
        for fn in files:
            tag = cache_tag(fn)
            r = load_reco_cache(cache / f"reco_{tag}.npz")
            recos.append(r); srcs.append(str(cache / f"reco_{tag}.npz"))
            pot += read_pot(data_dir / fn)["pot_used"]
            if is_mc:
                t = TruthTable.load(cache / f"reco_{tag}_truthcols.npz")
                if not t.has_fs:
                    raise ValueError(f"{cache / f'reco_{tag}_truthcols.npz'} has no final-state particles: rebuild the cache "
                                     f"(`python -m ndp data cache --channel {ch.name} --reco-only`)")
                truths.append(t)
        keys = [k for k in recos[0] if k != "__meta__"]
        r = {k: np.concatenate([x[k] for x in recos]) for k in keys}
        r["__meta__"] = recos[0].get("__meta__", {})
        return r, (TruthTable.concatenate(truths) if truths else None), pot, srcs
    rd, _, pot_d, src_d = concat(ch.data["reco_data_files"], False)
    rm, rt, pot_m, src_m = concat(ch.data["reco_mc_files"], True)
    truth = TruthTable.concatenate([TruthTable.load(cache / f"truth_{cache_tag(fn)}.npz") for fn in ch.data["reco_mc_files"]])
    return rd, rm, rt, truth, pot_d, pot_m, src_d + src_m


def run_selection_comparison(channel_name: str, cfg: SiteConfig | None = None, *, out_root: str | Path | None = None,
                             slug: str | None = None) -> Path:
    from .channels import list_measurements, load_measurement
    from .channels.selections import cutflow as sel_cutflow, select
    from .compare.plots import DATA_COLOR
    cfg = cfg or load_site_config()
    ch = load_channel(channel_name)
    rd, rm, rt, truth, pot_d, pot_m, sources = _load_reco_side(ch, cfg)
    scale = pot_d / pot_m
    params = ch.observable_params

    # ---- cutflow with purity and efficiency ---------------------------------------------------
    steps_d, steps_m = sel_cutflow(ch, rd), sel_cutflow(ch, rm)
    sig_rt = ch.is_signal(rt) & ch.in_phase_space(rt)
    den = int((ch.is_signal(truth) & ch.in_phase_space(truth)).sum())
    rows = []
    for (lab, md), (_, mm) in zip(steps_d, steps_m):
        n_sig = int((mm & sig_rt).sum())
        rows.append({"step": lab, "n_data": int(md.sum()), "n_mc": int(mm.sum()), "n_mc_scaled": float(mm.sum() * scale),
                     "mc_purity": float(n_sig / max(mm.sum(), 1)), "mc_efficiency": float(n_sig / max(den, 1)),
                     "data_over_mc": float(md.sum() / max(mm.sum() * scale, 1e-12))})
    sel_d, sel_m = steps_d[-1][1], steps_m[-1][1]
    cats = mc_categories(ch, rt)
    comp = {c: {"n_mc": int((sel_m & (cats == c)).sum()), "fraction": float((sel_m & (cats == c)).sum() / max(sel_m.sum(), 1))}
            for c in _CATEGORIES}

    # ---- proton-score scan (all other cuts as in the manifest) --------------------------------
    scan = []
    p0 = ch.selection.get("params", {})
    if "proton_score1_min" in p0:
        from .adapters.minerva_anatuple import ccqelike_1mu1p_cutflow
        for thr in (0.0, 0.2, 0.35, 0.5, 0.6, 0.7, 0.8):
            pp = dict(p0, proton_score1_min=thr)
            pd_, _ = ccqelike_1mu1p_cutflow(rd, pp); pm_, _ = ccqelike_1mu1p_cutflow(rm, pp)
            ns = int((pm_ & sig_rt).sum())
            scan.append({"proton_score1_min": thr, "n_data": int(pd_.sum()), "n_mc_scaled": float(pm_.sum() * scale),
                         "mc_purity": float(ns / max(pm_.sum(), 1)), "mc_efficiency": float(ns / max(den, 1))})

    # ---- run directory + figures ---------------------------------------------------------------
    run_dir = unique_run_dir(out_root or cfg.runs, slug or f"selection_{ch.name}")
    figs = ensure_dir(run_dir / "figs")
    dump_json(ch.to_dict(), run_dir / "channel.json")
    dump_json({"channel": ch.name, "selection": ch.selection, "pot_data": pot_d, "pot_mc": pot_m, "scale": scale,
               "efficiency_denominator": den, "rows": rows, "score_scan": scan}, run_dir / "cutflow.json")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    per_meas = {}
    fig_paths = []
    names = [n for n in list_measurements(ch) if n != "published"]
    for name in names:
        try:
            m = load_measurement(ch, name)
            xd = m.reco_observables(rd, params=params)[0][sel_d]
            xm = m.reco_observables(rm, params=params)[0][sel_m]
        except KeyError as e:
            per_meas[name] = {"error": str(e)}
            continue
        edges = np.asarray(m.x.edges, float)
        hd, _ = np.histogram(xd[np.isfinite(xd)], bins=edges)
        hm = {c: np.histogram(xm[np.isfinite(xm) & (cats[sel_m] == c)], bins=edges)[0] for c in _CATEGORIES}
        tot = sum(hm.values()).astype(float)
        per_meas[name] = {"edges": edges.tolist(), "data": hd.tolist(), "mc_scaled": (tot * scale).tolist(),
                          "mc_by_category_scaled": {c: (v * scale).tolist() for c, v in hm.items()},
                          "n_data_in_range": int(hd.sum()), "n_mc_scaled_in_range": float(tot.sum() * scale),
                          "n_data_nan": int((~np.isfinite(xd)).sum()), "n_mc_nan": int((~np.isfinite(xm)).sum())}
        fig, (ax, axr) = plt.subplots(2, 1, figsize=(7, 6.5), gridspec_kw={"height_ratios": [3, 1]}, sharex=True)
        bottom = np.zeros(len(edges) - 1)
        for c in _CATEGORIES:
            v = hm[c] * scale
            if v.sum() > 0:
                ax.bar(edges[:-1], v, width=np.diff(edges), bottom=bottom, align="edge", color=_CAT_COLORS[c], label=c, lw=0)
                bottom += v
        centres = 0.5 * (edges[:-1] + edges[1:])
        ax.errorbar(centres, hd, xerr=np.diff(edges) / 2, yerr=np.sqrt(hd), fmt="o", ms=4, color=DATA_COLOR, capsize=2, label="data")
        ax.set_ylabel("events per bin (MC scaled to data POT)")
        ax.set_title(f"{ch.name}: {m.x.label or m.x.reco}  ({int(hd.sum())} data, {tot.sum() * scale:.1f} MC)", fontsize=10)
        ax.legend(fontsize=7, ncol=2)
        mc_tot = tot * scale
        with np.errstate(invalid="ignore", divide="ignore"):
            ratio = np.where(mc_tot > 0, hd / mc_tot, np.nan); rerr = np.where(mc_tot > 0, np.sqrt(hd) / mc_tot, np.nan)
            band = np.where(mc_tot > 0, np.sqrt(tot) * scale / mc_tot, np.nan)
        axr.fill_between(edges, 1 - np.append(band, band[-1]), 1 + np.append(band, band[-1]), step="post", color="#9e9e9e", alpha=0.3, lw=0, label="MC stat")
        axr.errorbar(centres, ratio, xerr=np.diff(edges) / 2, yerr=rerr, fmt="o", ms=4, color=DATA_COLOR, capsize=2)
        axr.axhline(1, color="k", lw=0.8, ls="--"); axr.set_ylim(0, 2.5); axr.set_ylabel("data / MC")
        axr.set_xlabel(f"{m.x.label or m.x.reco} [{m.x.units}]" if m.x.units else (m.x.label or m.x.reco))
        fig.tight_layout()
        p = figs / f"data_vs_mc_{name}.png"
        fig.savefig(p, dpi=120); plt.close(fig); fig_paths.append(p)

    summary = {"channel": ch.name, "selection": ch.selection.get("name"), "pot_data": pot_d, "pot_mc": pot_m, "scale": scale,
               "n_data_selected": int(sel_d.sum()), "n_mc_selected": int(sel_m.sum()), "n_mc_selected_scaled": float(sel_m.sum() * scale),
               "mc_purity": rows[-1]["mc_purity"], "mc_efficiency": rows[-1]["mc_efficiency"], "efficiency_denominator": den,
               "paper_reference": {"tracker_efficiency": 0.28, "tracker_purity": 0.60, "source": "arXiv:2503.15047 Table I / Sec. Efficiency (paper_2503.15047.md Sec. 4)"},
               "mc_composition_selected": comp, "score_scan": scan, "measurements": per_meas}
    dump_json(summary, run_dir / "summary.json")

    lines = [f"# Reco selection `{ch.selection.get('name')}` on `{ch.name}`: data vs official MC", "",
             f"POT data {pot_d:.4g}, MC {pot_m:.4g} (scale {scale:.5g}). Efficiency denominator = truth signal in the fiducial volume: {den} events.", "",
             "## Cutflow (cumulative; MC purity = signal fraction, efficiency = selected signal / denominator)", "",
             "| step | data | MC (scaled) | data/MC | MC purity | MC efficiency |", "|---|---|---|---|---|---|"]
    lines += [f"| {r['step']} | {r['n_data']} | {r['n_mc_scaled']:.1f} | {r['data_over_mc']:.3f} | {r['mc_purity']:.3f} | {r['mc_efficiency']:.3f} |" for r in rows]
    lines += ["", f"Paper (CH tracker): efficiency 28%, purity 60% (arXiv:2503.15047). This run: efficiency {rows[-1]['mc_efficiency']:.3f}, purity {rows[-1]['mc_purity']:.3f}.", "",
              "## Composition of the selected MC (categories from the reco rows' truth)", "", "| category | n MC | fraction |", "|---|---|---|"]
    lines += [f"| {c} | {v['n_mc']} | {v['fraction']:.3f} |" for c, v in comp.items()]
    if scan:
        lines += ["", "## Proton-score threshold scan (all other cuts as in the manifest)", "",
                  "| proton_score1 > | data | MC (scaled) | MC purity | MC efficiency |", "|---|---|---|---|---|"]
        lines += [f"| {s['proton_score1_min']} | {s['n_data']} | {s['n_mc_scaled']:.1f} | {s['mc_purity']:.3f} | {s['mc_efficiency']:.3f} |" for s in scan]
    lines += ["", "## Data vs MC per released grid (selected sample)", "", "| measurement | data in range | MC scaled in range | data NaN | MC NaN |", "|---|---|---|---|---|"]
    lines += [f"| {k} | {v['n_data_in_range']} | {v['n_mc_scaled_in_range']:.1f} | {v['n_data_nan']} | {v['n_mc_nan']} |" if "data" in v else f"| {k} | ERROR {v['error']} | | | |"
              for k, v in per_meas.items()]
    lines += ["", "## Figures", ""] + [f"![{p.stem}]({p.relative_to(run_dir)})" for p in fig_paths]
    (run_dir / "report.md").write_text("\n".join(lines) + "\n")

    data_dir = cfg.require("data_dir")
    manifest = {"run_id": run_dir.name, "kind": "selection_comparison", "timestamp": timestamp(), "platform_version": __version__,
                "channel": ch.name, "channel_file": str(ch.path), "channel_sha256": sha256_text(Path(ch.path).read_text()) if ch.path else None,
                "selection": ch.selection, "inputs": {"caches": [cheap_fingerprint(s) for s in sources],
                                                      "data_files": [cheap_fingerprint(data_dir / fn) for fn in ch.data["reco_data_files"]],
                                                      "mc_files": [cheap_fingerprint(data_dir / fn) for fn in ch.data["reco_mc_files"]]},
                "pot": {"data": pot_d, "mc": pot_m, "scale": scale}, "git": git_state(cfg.repo_root), "versions": versions(),
                "outputs": ["channel.json", "cutflow.json", "summary.json", "report.md"] + [str(p.relative_to(run_dir)) for p in fig_paths],
                "results": {"n_data_selected": int(sel_d.sum()), "n_mc_selected": int(sel_m.sum()), "mc_purity": rows[-1]["mc_purity"],
                            "mc_efficiency": rows[-1]["mc_efficiency"]}}
    dump_json(manifest, run_dir / "manifest.json")
    return run_dir
