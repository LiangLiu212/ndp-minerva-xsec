#!/usr/bin/env python3
"""Render the GiBUU-overlay report: efficiency maps, ansatz closure, background composition and the
per-grid data vs (GiBUU signal + MC background) comparison, from an `ndp efficiency run` directory and
an `ndp run ... --measurement all` directory.

    python report/make_overlay_report.py --efficiency runs/<date>_efficiency_minerva_me_ccqelike_1mu1p \
        --overlay runs/<date>_gibuu_2025_me_fhc_c12__minerva_me_ccqelike_1mu1p__all18 --out report/GiBUU_overlay_FHC.md [--tag gibuu]

Every number is read from the run directories; figures are copied to <out dir>/figs/<tag>_*.png.
"""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def _splines_path() -> str | None:
    """The GENIE spline XML from the site config, for the sigma(E) cross-check (optional)."""
    try:
        import yaml
        p = str(yaml.safe_load((ROOT / "ndp.yaml").read_text()).get("genie_splines") or "")
        return p if p and Path(p).exists() else None
    except Exception:  # noqa: BLE001
        return None


GRID_ORDER = ["muon_p", "muon_theta", "muon_pt", "proton_p", "proton_theta", "proton_pt", "dpt", "dpt_fine", "dptx", "dpty",
              "alpha", "phi", "pl", "pn", "muon_costheta", "proton_costheta", "muon_p_costheta", "proton_p_costheta"]
LABELS = {"muon_p": "muon p", "muon_theta": "muon θ", "muon_pt": "muon p_T", "proton_p": "leading proton p", "proton_theta": "leading proton θ",
          "proton_pt": "leading proton p_T", "dpt": "δp_T", "dpt_fine": "δp_T (fine)", "dptx": "δp_Tx", "dpty": "δp_Ty", "alpha": "δα_T",
          "phi": "φ_T", "pl": "δp_L", "pn": "p_n", "muon_costheta": "cos θ_μ", "proton_costheta": "cos θ_p",
          "muon_p_costheta": "muon p × cos θ_μ", "proton_p_costheta": "proton p × cos θ_p"}


def n(x) -> str:
    return f"{int(round(float(x))):,}".replace(",", " ")


def sci(x: float, digits: int = 4) -> str:
    sup = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")
    m, e = f"{x:.{digits - 1}e}".split("e")
    return f"{m} × 10{str(int(e)).translate(sup)}"


def copy_fig(src: Path, dst: Path) -> str | None:
    if not src.exists():
        return None
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dst)
    return dst.name


def sigma_of_energy_figure(run_json: dict, out: Path, cfg_repo_root: Path, splines: str | None) -> Path | None:
    """The generated sample's own sigma_CC(E) per energy point, against the GENIE spline on the same target."""
    pts = run_json.get("energy_points")
    if not pts:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    E = np.array([q["energy"] for q in pts]); sig = np.array([q["sigma_cc_1e-38cm2"] for q in pts])
    nev = np.array([q["n_events"] for q in pts]); frac = np.array([q["flux_fraction"] for q in pts])
    fig, (ax, axf) = plt.subplots(2, 1, figsize=(7.5, 7), gridspec_kw={"height_ratios": [3, 1]}, sharex=True)
    ax.plot(E, sig, "o", ms=4, color="#08519c", label="GiBUU σ$_{CC}$/nucleon (per energy point)")
    if splines:
        try:
            from ndp.theory import splines as spl
            s = spl.total_xsec_splines(splines, 14, [1000060120], proc_contains="Weak[CC]", cache_dir=str(cfg_repo_root / "runs/_generator_cache/spline_cache"))
            Eg, sg = spl.per_nucleon_total(s, {1000060120: 1.0})
            m = (Eg >= E.min() * 0.8) & (Eg <= E.max() * 1.2)
            ax.plot(Eg[m], sg[m] / 1e-38, "-", lw=1.4, color="#a63603", label="GENIE G18_02a spline (¹²C, per nucleon)")
        except Exception:  # noqa: BLE001
            pass
    ax.set_ylabel("σ$_{CC}$ per nucleon [10⁻³⁸ cm²]"); ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax.set_title("GiBUU energy scan: cross section per point and the flux weight", fontsize=10)
    axf.bar(E, frac, width=np.diff(E).min() * 0.9, color="#9e9e9e", label="flux fraction of the point")
    axf2 = axf.twinx(); axf2.plot(E, nev, ".", ms=4, color="#2ca25f"); axf2.set_ylabel("events", color="#2ca25f", fontsize=8)
    axf.set_xlabel("E$_ν$ [GeV]"); axf.set_ylabel("flux fraction", fontsize=8); axf.set_yscale("log")
    fig.tight_layout(); fig.savefig(out, dpi=120); plt.close(fig)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--efficiency", required=True); ap.add_argument("--overlay", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--tag", default="gibuu")
    a = ap.parse_args()
    E, O = Path(a.efficiency).resolve(), Path(a.overlay).resolve()
    out = Path(a.out).resolve(); figs = out.parent / "figs"; rel = figs.name; tag = a.tag
    es = json.loads((E / "summary.json").read_text()); clo = json.loads((E / "ansatz_closure.json").read_text()) if (E / "ansatz_closure.json").exists() else None
    sc = json.loads((O / "scorecard.json").read_text()); om = json.loads((O / "manifest.json").read_text())
    ts = json.loads((O / "truth_summary.json").read_text())
    model = sc["model"]; prov = sc.get("provenance", {}); meta = ts.get("meta", {})
    grids = [g for g in GRID_ORDER if g in sc["measurements"]] + [g for g in sc["measurements"] if g not in GRID_ORDER]
    first = next(iter(sc["measurements"].values()))
    pot_data = first["pot_data"]; exp = first["expected"]

    L = [f"# GiBUU 1μ1p signal on the MINERvA ME FHC selection: efficiency, background and data overlay", "",
         f"**Status:** rendered from `{E.relative_to(ROOT)}` (efficiency maps, background, ansatz closure) and `{O.relative_to(ROOT)}` "
         f"(data vs GiBUU signal + MC background on {len(grids)} grids), manifest git `{(om.get('platform_git') or {}).get('sha', '')[:7]}`. "
         "Numbers are quoted from those directories; this file was rendered by `report/make_overlay_report.py`.", "",
         "**Model:** " + model.get("description", "").strip(), "", "---", "",
         "## 1. Inputs", "", "| item | value |", "|---|---|",
         f"| data | ME FHC playlists 1A–1P, {sci(pot_data)} POT_Used, {n(first['n_data_selected'])} selected events |",
         f"| official MC (efficiency and background) | {sci(es['pot_mc'])} POT, {n(es['n_signal_fiducial'])} truth signal events in the fiducial volume, {n(es['n_selected_signal'])} selected (⟨ε⟩ = {es['mean_efficiency']:.4f}) |",
         f"| GiBUU sample | {meta.get('generator', 'GiBUU')}, {n(meta.get('n_generated', ts['summary'].get('n_events', 0)))} events from {meta.get('n_jobs_merged', '?')} job(s), fingerprint `{meta.get('fingerprint', '')}` |",
         f"| GiBUU σ_CC (flux-averaged, per nucleon) | {meta.get('sigma_flux_avg_per_nucleon_cm2', float('nan')):.4e} cm² |",
         f"| normalisation | {exp['how']} |", ""]

    # the generated sample itself
    run_json = None
    src = meta.get("source")
    if src and (Path(src) / "gibuu_run.json").exists():
        run_json = json.loads((Path(src) / "gibuu_run.json").read_text())
    if run_json and run_json.get("energy_points"):
        pts = run_json["energy_points"]
        f = sigma_of_energy_figure(run_json, figs / f"{tag}_sigma_of_energy.png", ROOT, _splines_path())
        L += ["## 1b. The GiBUU sample", "",
              f"{len(pts)} energy points {pts[0]['energy']:.2f}–{pts[-1]['energy']:.2f} GeV covering "
              f"{run_json['flux_fraction_covered']:.4f} of the 0–100 GeV flux, {sum(q['n_jobs'] for q in pts)} grid jobs, "
              f"{run_json['n_generated']:,} events; flux-averaged σ$_{{CC}}$ = {run_json['sigma_flux_avg_per_nucleon_cm2']:.4e} cm²/nucleon "
              f"(over the covered flux). Missing points: {run_json.get('energy_points_missing') or 'none'}.".replace(",", " "), ""]
        if f:
            L += [f"![sigma of energy]({rel}/{f})", "*σ$_{CC}$(E) of every point against the GENIE spline, with each point's flux weight and event count.*", ""]

    # efficiency maps
    L += ["## 2. Selection efficiency (official MC)", "",
          "ε = selected truth signal / truth signal in the fiducial volume, per true cell (the paper's definition), binomial errors. "
          "The two maps are what a truth sample can be weighted with; the 1D grids are their projections and the released grids' efficiencies.", ""]
    for name in ("muon_p_costheta", "proton_p_costheta"):
        p = E / f"eff_{name}.json"
        if not p.exists():
            continue
        t = json.loads(p.read_text())
        f1 = copy_fig(E / "figs" / f"eff_map_{name}.png", figs / f"{tag}_eff_map_{name}.png")
        f2 = copy_fig(E / "figs" / f"eff_proj_{name}.png", figs / f"{tag}_eff_proj_{name}.png")
        L += [f"### ε({t['x']['label']}, {t['y']['label']})", "",
              "| " + t["x"]["label"] + " \\ " + t["y"]["label"] + " | " + " | ".join(f"[{lo:.4g}, {hi:.4g})" for lo, hi in zip(t["y"]["edges"][:-1], t["y"]["edges"][1:])) + " |",
              "|---|" + "---|" * (len(t["y"]["edges"]) - 1)]
        for i, (lo, hi) in enumerate(zip(t["x"]["edges"][:-1], t["x"]["edges"][1:])):
            L.append(f"| [{lo:g}, {hi:g}) | " + " | ".join(f"{t['eff'][i][j]:.3f} ± {t['err'][i][j]:.3f}" for j in range(len(t["y"]["edges"]) - 1)) + " |")
        L += ["", f"![{name} map]({rel}/{f1})" if f1 else "", f"![{name} projections]({rel}/{f2})" if f2 else "", ""]
    L += ["| grid | ⟨ε⟩ | background at data POT |", "|---|---|---|"]
    for g in grids:
        m = es["measurements"].get(g)
        if m:
            L.append(f"| {LABELS.get(g, g)} | {m['eff']:.4f} | {m['n_background_at_data_pot']:.0f} |")
    L.append("")
    if clo:
        an = es["ansatz"]
        L += ["## 3. The factorised ansatz w = ε_μ(p_μ, cos θ_μ) · ε_p(p_p, cos θ_p) / ⟨ε⟩ — closure on the MC", "",
              f"Weighting the MC's own truth signal with the maps and comparing with the actually selected signal per true bin: "
              f"total ratio {an['total_ratio_ansatz_over_actual']:.4f}. Per grid (the error the shortcut makes on that grid):", "",
              "| grid | ansatz / actual | max |rel. dev.| (bins > 50 events) |", "|---|---|---|"]
        for g in grids:
            r = an["per_grid_total_ratio"].get(g); d = an["per_grid_max_abs_rel_dev"].get(g)
            if r is not None:
                L.append(f"| {LABELS.get(g, g)} | {r:.4f} | {d:.3f} |" if d is not None else f"| {LABELS.get(g, g)} | {r:.4f} | – |")
        L.append("")
        for g in ("dpt", "alpha", "pn", "muon_theta", "proton_theta"):
            f = copy_fig(E / "figs" / f"ansatz_{g}.png", figs / f"{tag}_ansatz_{g}.png")
            if f:
                L.append(f"![ansatz {g}]({rel}/{f})")
        L.append("")
    # background
    L += ["## 4. Background of the selected sample (official MC, at the data POT)", "",
          "Selected reco events that are not truth signal inside the fiducial volume, by the category of the reco rows' truth, plus feed-in "
          "(signal whose true observable lies outside the grid). The MC is the unweighted central value.", "",
          "| grid | total | single π± | single π⁰ | multi-π | no pion | feed-in |", "|---|---|---|---|---|---|---|"]
    import numpy as np
    for g in grids:
        p = E / f"background_{g}.json"
        if not p.exists():
            continue
        b = json.loads(p.read_text()); c = b["by_category_at_data_pot"]
        L.append(f"| {LABELS.get(g, g)} | {b['n_at_data_pot']:.0f} | " + " | ".join(f"{np.asarray(c[k]).sum():.0f}" if k in c else "–" for k in
                 ("bkg 1 pi+-", "bkg 1 pi0", "bkg multi-pi", "bkg other (no pion)")) + f" | {np.asarray(b['feedin_at_data_pot']).sum():.0f} |")
    L.append("")
    # overlay
    L += ["## 5. Data versus GiBUU signal + MC background", "",
          "Three predictions per grid: **full** = GiBUU truth cells × efficiency × migration (the binned response learned on the official MC) + background; "
          "**eff-only** = GiBUU truth cells × efficiency, no migration; **ansatz** = GiBUU truth events × w from the two maps, histogrammed in truth bins. "
          "−2lnL is the Baker–Cousins Poisson likelihood ratio of the data against the full prediction.", "",
          "| grid | data | full | eff-only | ansatz | background | data/full | −2lnL/ndf (full) |", "|---|---|---|---|---|---|---|---|"]
    for g in grids:
        r = sc["measurements"][g]; tf, te, ta = r["totals"]["full"], r["totals"]["eff_only"], r["totals"]["ansatz"]; gf = r["gof"]["full"]
        L.append(f"| {LABELS.get(g, g)} | {tf['data']:.0f} | {tf['pred']:.0f} | {te['pred']:.0f} | {ta['pred']:.0f} | {tf['bkg']:.0f} | {tf['ratio_data_over_pred']:.3f} | {gf['minus2lnL']:.1f}/{gf['ndf']} |"
                 if ta else f"| {LABELS.get(g, g)} | {tf['data']:.0f} | {tf['pred']:.0f} | {te['pred']:.0f} | – | {tf['bkg']:.0f} | {tf['ratio_data_over_pred']:.3f} | {gf['minus2lnL']:.1f}/{gf['ndf']} |")
    L.append("")
    k = 1
    for g in grids:
        f = copy_fig(O / "figs" / f"overlay_{g}.png", figs / f"{tag}_overlay_{g}.png") or copy_fig(O / "figs" / f"folded_{g}.png", figs / f"{tag}_folded_{g}.png")
        if f:
            L += [f"![{g}]({rel}/{f})", f"*Figure {k} — {LABELS.get(g, g)}.*", ""]; k += 1
    L += ["## 6. Caveats", "",
          "- The background and the efficiency come from the unweighted central-value official MC (no MINERvA tune, flux or detector weights); the data/MC "
          "ratio of the selected sample is 0.84 with that MC, so absolute agreement or disagreement of the overlay carries that uncertainty.",
          "- The fiducial volume and n_nucleons are the inclusive channel's (the paper's own are unstated); the flux integral is the inclusive channel's Φ.",
          "- GiBUU: carbon only (hydrogen gives no signal), the channel flux rebinned to 0.5 GeV, card physics from GiBUU's MINERvA-ME card, statistics and "
          "weights as recorded in the sample's `gibuu_run.json`.",
          "- The ansatz prediction has no detector migration; its closure table (Sec. 3) is the size of that approximation on each grid.", ""]
    out.write_text("\n".join(x for x in L if x is not None) + "\n")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
