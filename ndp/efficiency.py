"""Selection-efficiency maps, the background of the selected sample, and the factorised-efficiency ansatz.

Everything here is read from binned surrogates built on the official MC (`ndp surrogate build`): for a
measurement grid, `BinnedResponse.eff` is the selected truth signal over the truth signal in the
fiducial volume per true cell (the paper's efficiency definition), `den_counts`/`num_counts` give its
binomial error, and `bkg_per_pot` / `bkg_by_category_counts` the selected non-signal MC per reco cell.

The 2D maps eps_mu(p_mu, cos theta_mu) and eps_p(p_p, cos theta_p) can be applied to any truth sample
(e.g. GiBUU) as a per-event weight w = eps_mu * eps_p / <eps>. That factorised ansatz is an
approximation; `ansatz_closure` measures it on the MC itself: the ansatz-weighted truth signal versus
the actually selected signal, per released grid.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from .channels import load_channel, load_measurement, list_measurements
from .config import SiteConfig, load_site_config
from .io import cheap_fingerprint, dump_json, ensure_dir, git_state, sha256_text, timestamp, unique_run_dir, versions
from .pipeline import resolve_surrogate
from .surrogate.chunked import BKG_CATEGORIES
from . import __version__

DEFAULT_MAPS = ("muon_p_costheta", "proton_p_costheta", "muon_costheta", "proton_costheta")


# ---- efficiency tables -----------------------------------------------------------------------------
def efficiency_table(measurement, surrogate, sur_path: Path) -> dict:
    """Edges, eff, binomial error, den, num as JSON-ready grids ([n_x][n_y]; 1D grids have n_y = 1)."""
    b = measurement.binning
    err = surrogate.eff_uncertainty()
    den, num = surrogate.den_counts, surrogate.num_counts
    total = {"den": float(den.sum()), "num": float(num.sum()), "eff": float(num.sum() / den.sum()) if den.sum() else None}
    return {"measurement": measurement.name, "channel": measurement.channel, "is_1d": measurement.is_1d,
            "x": {"observable": measurement.x.truth, "label": measurement.x.label or measurement.x.truth, "units": measurement.x.units, "edges": list(b.x_edges)},
            "y": {"observable": measurement.y.truth, "label": measurement.y.label or measurement.y.truth, "units": measurement.y.units, "edges": list(b.y_edges)},
            "eff": b.to_grid(surrogate.eff).tolist(), "err": b.to_grid(err).tolist(),
            "den": b.to_grid(den).tolist(), "num": b.to_grid(num).tolist(),
            "projection_x": _projection(b, den, num, "x"), "projection_y": _projection(b, den, num, "y") if not measurement.is_1d else None,
            "total": total, "pot_mc": surrogate.meta.get("pot_mc"), "surrogate": str(sur_path),
            "definition": "eff = selected truth signal / truth signal in the fiducial volume, per true cell (unweighted MC counts); "
                          "err = sqrt(eff (1 - eff) / den)"}


def _projection(b, den, num, axis: str) -> dict:
    d = b.project(den, axis, per_width=False); n = b.project(num, axis, per_width=False)
    with np.errstate(invalid="ignore", divide="ignore"):
        eff = np.where(d > 0, n / d, 0.0); err = np.where(d > 0, np.sqrt(eff * (1 - eff) / np.where(d > 0, d, 1)), 0.0)
    edges = b.x_edges if axis == "x" else b.y_edges
    return {"edges": list(edges), "den": d.tolist(), "num": n.tolist(), "eff": eff.tolist(), "err": err.tolist()}


def background_table(measurement, surrogate, pot_data: float) -> dict:
    """The selected non-signal MC per reco cell at the data POT: total (feed-in included), by category, feed-in."""
    b = measurement.binning
    pot_mc = float(surrogate.meta["pot_mc"]); scale = pot_data / pot_mc
    total_counts = surrogate.bkg_per_pot * pot_mc
    by_cat = surrogate.bkg_by_category_counts or {}
    feed = surrogate.feedin_counts if surrogate.feedin_counts is not None else np.zeros(b.n_cells)
    check = None
    if by_cat:
        s = sum(by_cat.values()) + feed
        check = {"sum_categories_plus_feedin_equals_total": bool(np.allclose(s, total_counts, atol=1e-6)),
                 "max_abs_dev": float(np.abs(s - total_counts).max())}
    return {"measurement": measurement.name, "pot_data": pot_data, "pot_mc": pot_mc, "scale": scale,
            "x_edges": list(b.x_edges), "y_edges": list(b.y_edges), "is_1d": measurement.is_1d,
            "total_at_data_pot": b.to_grid(total_counts * scale).tolist(),
            "mc_stat_err_at_data_pot": b.to_grid(np.sqrt(total_counts) * scale).tolist(),
            "by_category_at_data_pot": {c: b.to_grid(v * scale).tolist() for c, v in by_cat.items()},
            "feedin_at_data_pot": b.to_grid(feed * scale).tolist(),
            "mc_counts_total": float(total_counts.sum()), "n_at_data_pot": float(total_counts.sum() * scale),
            "consistency": check,
            "definition": "selected reco events that are not truth signal inside the fiducial volume (by category, from the reco rows' truth) "
                          "plus feed-in (selected signal whose true observable is outside the grid); MC counts scaled by POT_data / POT_mc"}


# ---- factorised ansatz --------------------------------------------------------------------------------
class EfficiencyMaps:
    """eps_mu(p, cos theta) and eps_p(p, cos theta) from two efficiency tables, applied per event."""

    def __init__(self, mu: dict, pr: dict, mean_eff: float):
        self.mu, self.pr, self.mean_eff = mu, pr, mean_eff

    @staticmethod
    def _lookup(tab: dict, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        xe, ye = np.asarray(tab["x"]["edges"]), np.asarray(tab["y"]["edges"])
        grid = np.asarray(tab["eff"])
        ix = np.searchsorted(xe, x, side="right") - 1
        iy = np.searchsorted(ye, y, side="right") - 1
        inside = (x >= xe[0]) & (x < xe[-1]) & (y >= ye[0]) & (y < ye[-1])
        out = np.zeros(len(x))
        out[inside] = grid[ix[inside], iy[inside]]
        return out

    def weights(self, channel, t) -> np.ndarray:
        """w = eps_mu(p_mu, cos th_mu) x eps_p(p_p, cos th_p) / <eps>; 0 outside either map (NaN proton -> 0)."""
        pm = channel.evaluate(self.mu["x"]["observable"], t); cm = channel.evaluate(self.mu["y"]["observable"], t)
        pp = channel.evaluate(self.pr["x"]["observable"], t); cp = channel.evaluate(self.pr["y"]["observable"], t)
        ok = np.isfinite(pm) & np.isfinite(cm) & np.isfinite(pp) & np.isfinite(cp)
        w = np.zeros(t.n)
        w[ok] = self._lookup(self.mu, pm[ok], cm[ok]) * self._lookup(self.pr, pp[ok], cp[ok]) / self.mean_eff
        return w

    @staticmethod
    def load(run_dir: str | Path) -> "EfficiencyMaps":
        d = Path(run_dir)
        mu = json.loads((d / "eff_muon_p_costheta.json").read_text()); pr = json.loads((d / "eff_proton_p_costheta.json").read_text())
        s = json.loads((d / "summary.json").read_text())
        return EfficiencyMaps(mu, pr, float(s["mean_efficiency"]))


def ansatz_closure(channel, measurements: list, maps: EfficiencyMaps, surrogates: dict, cfg, log=print) -> dict:
    """Weight the MC's own truth signal with the ansatz and compare, per grid, with the selected signal per true cell."""
    from .products import iter_mc_chunks
    acc = {m.name: np.zeros(m.binning.n_cells) for m in measurements}
    n_sig = 0; sum_w = 0.0
    for label, rm, rt, truth, pot, srcs in iter_mc_chunks(cfg, channel):
        del rm, rt
        sig = channel.is_signal(truth) & channel.in_phase_space(truth)
        w = maps.weights(channel, truth)[sig]
        n_sig += int(sig.sum()); sum_w += float(w.sum())
        for m in measurements:
            x, y = m.truth_observables(channel, truth)
            acc[m.name] += m.binning.histogram(x[sig], y[sig], w)[0]
        log(f"ansatz closure {label}: {int(sig.sum())} signal events, sum of ansatz weights {w.sum():.1f}")
        del truth
    out = {"n_signal_fiducial": n_sig, "ansatz_selected_total": sum_w, "grids": {}}
    for m in measurements:
        num = surrogates[m.name].num_counts
        pred = acc[m.name]
        with np.errstate(invalid="ignore", divide="ignore"):
            ratio = np.where(num > 0, pred / num, np.nan)
        big = num > 50
        out["grids"][m.name] = {"edges_x": list(m.binning.x_edges), "edges_y": list(m.binning.y_edges),
                                "actual_selected_signal": m.binning.to_grid(num).tolist(), "ansatz": m.binning.to_grid(pred).tolist(),
                                "ratio": m.binning.to_grid(np.nan_to_num(ratio, nan=0.0)).tolist(),
                                "total_ratio": float(pred.sum() / num.sum()) if num.sum() else None,
                                "max_abs_rel_dev_cells_gt_50": float(np.nanmax(np.abs(ratio[big] - 1))) if big.any() else None,
                                "rms_rel_dev_cells_gt_50": float(np.sqrt(np.nanmean((ratio[big] - 1) ** 2))) if big.any() else None}
    out["actual_selected_total"] = float(surrogates[measurements[0].name].num_counts.sum()) if measurements else None
    return out


# ---- figures -----------------------------------------------------------------------------------------------
def _fig_map(tab: dict, out: Path, title: str):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    eff, err, den = np.asarray(tab["eff"]), np.asarray(tab["err"]), np.asarray(tab["den"])
    xe, ye = np.asarray(tab["x"]["edges"]), np.asarray(tab["y"]["edges"])
    fig, ax = plt.subplots(figsize=(9.5, 7))
    im = ax.pcolormesh(np.arange(len(xe)), np.arange(len(ye)), eff.T, vmin=0, vmax=max(0.6, float(np.nanmax(eff))), cmap="viridis")
    for i in range(eff.shape[0]):
        for j in range(eff.shape[1]):
            if den[i, j] > 0:
                ax.text(i + 0.5, j + 0.5, f"{eff[i, j]:.3f}\n±{err[i, j]:.3f}\n({int(den[i, j])})", ha="center", va="center", fontsize=6,
                        color="w" if eff[i, j] < 0.35 else "k")
    ax.set_xticks(np.arange(len(xe))); ax.set_xticklabels([f"{v:g}" for v in xe], fontsize=8)
    ax.set_yticks(np.arange(len(ye))); ax.set_yticklabels([f"{v:.4g}" for v in ye], fontsize=8)
    ax.set_xlabel(f"{tab['x']['label']} [{tab['x']['units']}]" if tab["x"]["units"] else tab["x"]["label"])
    ax.set_ylabel(f"{tab['y']['label']} [{tab['y']['units']}]" if tab["y"]["units"] else tab["y"]["label"])
    fig.colorbar(im, ax=ax, label="efficiency (selected signal / fiducial signal)")
    ax.set_title(title, fontsize=10); fig.tight_layout(); fig.savefig(out, dpi=130); plt.close(fig)
    return out


def _fig_proj(tab: dict, out: Path, title: str):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    projs = [("x", tab["projection_x"], tab["x"])] + ([("y", tab["projection_y"], tab["y"])] if tab.get("projection_y") else [])
    fig, axes = plt.subplots(1, len(projs), figsize=(6.5 * len(projs), 4.8), squeeze=False)
    for ax, (_, p, axis) in zip(axes[0], projs):
        e = np.asarray(p["edges"]); c = 0.5 * (e[:-1] + e[1:])
        ax.errorbar(c, p["eff"], xerr=np.diff(e) / 2, yerr=p["err"], fmt="o", ms=4, capsize=2, color="#08519c")
        ax.set_ylim(0, max(0.6, 1.15 * max(p["eff"]))); ax.set_ylabel("efficiency"); ax.grid(alpha=0.3)
        ax.set_xlabel(f"{axis['label']} [{axis['units']}]" if axis["units"] else axis["label"])
        for xc, v, n in zip(c, p["eff"], p["den"]):
            ax.annotate(f"{v:.3f}", (xc, v), textcoords="offset points", xytext=(0, 7), ha="center", fontsize=7)
    fig.suptitle(title, fontsize=10); fig.tight_layout(); fig.savefig(out, dpi=130); plt.close(fig)
    return out


def _fig_ansatz(name: str, g: dict, out: Path, label: str):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    act = np.asarray(g["actual_selected_signal"]).ravel(); ans = np.asarray(g["ansatz"]).ravel()
    e = np.asarray(g["edges_x"]) if len(act) == len(g["edges_x"]) - 1 else np.arange(len(act) + 1)
    c = 0.5 * (e[:-1] + e[1:])
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(6.5, 6), gridspec_kw={"height_ratios": [3, 1]}, sharex=True)
    ax.step(e, np.append(act, act[-1]), where="post", color="#08519c", lw=1.6, label="selected truth signal (actual, per true bin)")
    ax.step(e, np.append(ans, ans[-1]), where="post", color="#d62728", lw=1.4, ls="--", label="ansatz: signal × ε_μ ε_p / ⟨ε⟩")
    ax.set_ylabel("MC events"); ax.legend(fontsize=8); ax.set_title(f"factorised-efficiency closure: {name}", fontsize=10)
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where(act > 0, ans / act, np.nan)
    axr.step(e, np.append(r, r[-1]), where="post", color="#d62728", lw=1.4); axr.axhline(1, color="k", lw=0.8, ls="--")
    axr.set_ylim(0.6, 1.4); axr.set_ylabel("ansatz / actual"); axr.set_xlabel(label)
    fig.tight_layout(); fig.savefig(out, dpi=120); plt.close(fig)
    return out


# ---- the run ------------------------------------------------------------------------------------------------
def run_efficiency(channel_name: str, cfg: SiteConfig | None = None, *, maps=DEFAULT_MAPS, grids: list | None = None,
                   out_root: str | Path | None = None, slug: str | None = None, closure: bool = True, log=print) -> Path:
    from .products import sources_fingerprints, total_pot
    cfg = cfg or load_site_config()
    ch = load_channel(channel_name)
    names = [n for n in list_measurements(ch) if n != "published"] if grids is None else list(grids)
    names = list(dict.fromkeys([*maps, *names]))
    meas = {n: load_measurement(ch, n) for n in names}
    surs, paths = {}, {}
    for n, m in meas.items():
        s, p = resolve_surrogate(ch, m, cfg)
        if s is None or s.kind != "binned_response" or s.den_counts is None:
            raise FileNotFoundError(f"no binned surrogate for measurement {n!r}: run `ndp surrogate build --channel {ch.name} --measurement all --kind binned`")
        surs[n], paths[n] = s, p
    pot_data = total_pot(cfg, ch, "data")
    run_dir = unique_run_dir(out_root or cfg.runs, slug or f"efficiency_{ch.name}")
    figs = ensure_dir(run_dir / "figs")
    dump_json(ch.to_dict(), run_dir / "channel.json")
    tables, bkgs, fig_paths = {}, {}, []
    for n, m in meas.items():
        tab = efficiency_table(m, surs[n], paths[n]); tables[n] = tab
        dump_json(tab, run_dir / f"eff_{n}.json")
        np.savez_compressed(run_dir / f"eff_{n}.npz", x_edges=np.asarray(tab["x"]["edges"]), y_edges=np.asarray(tab["y"]["edges"]),
                            eff=np.asarray(tab["eff"]), err=np.asarray(tab["err"]), den=np.asarray(tab["den"]), num=np.asarray(tab["num"]))
        title = f"{ch.name}: efficiency vs {tab['x']['label']}" + (f" × {tab['y']['label']}" if not m.is_1d else "")
        if not m.is_1d:
            fig_paths.append(_fig_map(tab, figs / f"eff_map_{n}.png", title))
        fig_paths.append(_fig_proj(tab, figs / f"eff_proj_{n}.png", title))
        bk = background_table(m, surs[n], pot_data); bkgs[n] = bk
        dump_json(bk, run_dir / f"background_{n}.json")
        log(f"{n}: <eff> {tab['total']['eff']:.4f} over {int(tab['total']['den'])} signal events; background {bk['n_at_data_pot']:.0f} at data POT")
    ref = surs[names[0]]
    mean_eff = float(ref.num_counts.sum() / ref.den_counts.sum())
    mu_tab, pr_tab = tables.get("muon_p_costheta"), tables.get("proton_p_costheta")
    summary = {"channel": ch.name, "pot_data": pot_data, "pot_mc": ref.meta.get("pot_mc"), "mean_efficiency": mean_eff,
               "n_signal_fiducial": float(ref.den_counts.sum()), "n_selected_signal": float(ref.num_counts.sum()),
               "measurements": {n: {"eff": tables[n]["total"]["eff"], "n_background_at_data_pot": bkgs[n]["n_at_data_pot"],
                                    "surrogate": str(paths[n])} for n in names},
               "maps": list(maps), "ansatz": None}
    if closure and mu_tab and pr_tab:
        em = EfficiencyMaps(mu_tab, pr_tab, mean_eff)
        grids_for_closure = [meas[n] for n in names]
        clo = ansatz_closure(ch, grids_for_closure, em, surs, cfg, log=log)
        dump_json(clo, run_dir / "ansatz_closure.json")
        for n, g in clo["grids"].items():
            if meas[n].is_1d:
                fig_paths.append(_fig_ansatz(n, g, figs / f"ansatz_{n}.png", meas[n].x.axis_label("truth")))
        summary["ansatz"] = {"definition": "w = eps_mu(p_mu, cos th_mu) * eps_p(p_p, cos th_p) / <eps>, <eps> = selected signal / fiducial signal",
                             "mean_efficiency": mean_eff, "total_ratio_ansatz_over_actual": clo["ansatz_selected_total"] / clo["actual_selected_total"],
                             "per_grid_total_ratio": {n: g["total_ratio"] for n, g in clo["grids"].items()},
                             "per_grid_max_abs_rel_dev": {n: g["max_abs_rel_dev_cells_gt_50"] for n, g in clo["grids"].items()}}
    dump_json(summary, run_dir / "summary.json")
    _write_report(run_dir, ch, meas, tables, bkgs, summary, fig_paths)
    manifest = {"run_id": run_dir.name, "kind": "efficiency", "timestamp": timestamp(), "platform_version": __version__,
                "channel": ch.name, "channel_file": str(ch.path), "channel_sha256": sha256_text(Path(ch.path).read_text()) if ch.path else None,
                "surrogates": {n: {"path": str(paths[n]), **cheap_fingerprint(paths[n] / "arrays.npz")} for n in names},
                "inputs": {"fingerprints": sources_fingerprints(cfg, ch)}, "pot_data": pot_data, "git": git_state(cfg.repo_root), "versions": versions(),
                "outputs": ["summary.json", "report.md", "ansatz_closure.json"] + [f"eff_{n}.json" for n in names] + [f"background_{n}.json" for n in names]
                           + [str(p.relative_to(run_dir)) for p in fig_paths]}
    dump_json(manifest, run_dir / "manifest.json")
    return run_dir


def _write_report(run_dir: Path, ch, meas: dict, tables: dict, bkgs: dict, summary: dict, fig_paths: list):
    L = [f"# Selection efficiency and background of `{ch.name}` (selection `{ch.selection.get('name')}`)", "",
         f"Official MC {summary['pot_mc']:.4g} POT: {summary['n_signal_fiducial']:.0f} truth signal events in the fiducial volume, "
         f"{summary['n_selected_signal']:.0f} of them selected: mean efficiency {summary['mean_efficiency']:.4f}. "
         f"Backgrounds are scaled to the data POT {summary['pot_data']:.4g}. Efficiency = selected truth signal / truth signal in the "
         "fiducial volume per true cell (the paper's definition); error = binomial.", "",
         "## Efficiency per grid", "", "| measurement | grid | <eff> | den | background at data POT |", "|---|---|---|---|---|"]
    for n, t in tables.items():
        m = meas[n]; g = f"{len(t['x']['edges']) - 1}" + ("" if m.is_1d else f" × {len(t['y']['edges']) - 1}")
        L.append(f"| {n} | {g} | {t['total']['eff']:.4f} | {int(t['total']['den'])} | {bkgs[n]['n_at_data_pot']:.0f} |")
    for n in ("muon_p_costheta", "proton_p_costheta"):
        if n in tables:
            t = tables[n]
            L += ["", f"### {n}: ε({t['x']['label']}, {t['y']['label']})", "", "| " + t["x"]["label"] + " \\ " + t["y"]["label"] + " | "
                  + " | ".join(f"[{lo:.4g}, {hi:.4g})" for lo, hi in zip(t["y"]["edges"][:-1], t["y"]["edges"][1:])) + " |",
                  "|---|" + "---|" * (len(t["y"]["edges"]) - 1)]
            for i, (lo, hi) in enumerate(zip(t["x"]["edges"][:-1], t["x"]["edges"][1:])):
                L.append(f"| [{lo:g}, {hi:g}) | " + " | ".join(f"{t['eff'][i][j]:.3f} ± {t['err'][i][j]:.3f}" for j in range(len(t["y"]["edges"]) - 1)) + " |")
    if summary.get("ansatz"):
        a = summary["ansatz"]
        L += ["", "## Factorised ansatz w = ε_μ(p_μ, cos θ_μ) · ε_p(p_p, cos θ_p) / ⟨ε⟩ — closure on the MC", "",
              f"⟨ε⟩ = {a['mean_efficiency']:.4f}; total ansatz-selected / actually selected = {a['total_ratio_ansatz_over_actual']:.4f}.", "",
              "| grid | ansatz / actual (total) | max |rel. dev.| over bins with > 50 events |", "|---|---|---|"]
        for n, r in a["per_grid_total_ratio"].items():
            d = a["per_grid_max_abs_rel_dev"].get(n)
            L.append(f"| {n} | {r:.4f} | {d:.3f} |" if d is not None else f"| {n} | {r:.4f} | – |")
        L += ["", "The ansatz is exact by construction for the two map grids' totals; the per-grid deviations are the error made by "
              "weighting a truth sample with the maps instead of folding it through the full response."]
    L += ["", "## Background composition (selected non-signal MC at the data POT)", "", "| grid | total | " + " | ".join(BKG_CATEGORIES) + " | feed-in |",
          "|---|---|" + "---|" * (len(BKG_CATEGORIES) + 1)]
    for n, b in bkgs.items():
        cats = b["by_category_at_data_pot"]
        L.append(f"| {n} | {b['n_at_data_pot']:.0f} | " + " | ".join(f"{np.asarray(cats[c]).sum():.0f}" if c in cats else "–" for c in BKG_CATEGORIES)
                 + f" | {np.asarray(b['feedin_at_data_pot']).sum():.0f} |")
    L += ["", "## Figures", ""] + [f"![{p.stem}]({p.relative_to(run_dir)})" for p in fig_paths]
    (run_dir / "report.md").write_text("\n".join(L) + "\n")
