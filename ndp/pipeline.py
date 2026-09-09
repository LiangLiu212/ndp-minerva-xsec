"""Orchestration: one theorist model + one channel (+ one measurement) -> a manifest-backed run directory.

    run_dir/
      model.yaml|json     the model spec as run
      channel.json        the channel manifest snapshot
      measurement.json    the observable pair + binning the comparison was made on
      manifest.json       inputs (fingerprints), versions, git state, timings, results summary
      scorecard.json      every number produced
      report.md           human summary
      figs/               comparison figures
      truth_summary.json  the realised sample (counts, interaction mix, normalisation)

The primary comparison is *folded*: the model's truth is pushed through the detector
surrogate learned for this measurement and compared with the selected data counts. The
*unfolded* comparison against the experiment's published cross section is added when the
measurement is one the experiment published (the channel's `published` grid).
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from . import __version__
from .config import load_site_config, SiteConfig
from .io import (dump_json, dump_yaml_or_json, cheap_fingerprint, git_state, versions, timestamp, unique_run_dir,
                 ensure_dir)
from .channels import load_channel, ChannelSpec, Measurement, load_measurement
from .theory.models import ModelSpec, realize, RealizeContext, Prediction
from .surrogate.base import load_surrogate
from .compare import unfolded as unf, folded as fld, plots, report


def resolve_surrogate(channel: ChannelSpec, measurement: Measurement, cfg: SiteConfig, explicit: str | None = None):
    """Explicit path, else the measurement's `surrogate.default`/`fallback`, else the conventional
    surrogates/<channel>[/<measurement>]/binned_<tag> and parametric_<tag> directories."""
    cands = [explicit] if explicit else [measurement.surrogate.get("default"), measurement.surrogate.get("fallback")]
    if not explicit:
        root = measurement.surrogate_root(cfg)
        if root.exists():
            cands += [str(p) for p in sorted(root.glob("binned_*"))] + [str(p) for p in sorted(root.glob("parametric_*"))]
    for c in cands:
        if not c:
            continue
        p = Path(c)
        if not p.is_absolute():
            p = cfg.repo_root / p
        if (p / "surrogate.json").exists():
            return load_surrogate(p), p
    return None, None


def run_model(model: str | Path | ModelSpec, channel_name: str, *, measurement: str | Path | None = None,
              cfg: SiteConfig | None = None, out_root: str | Path | None = None, surrogate_path: str | None = None,
              modes: tuple = ("folded", "unfolded"), fold_events: bool = False, slug: str | None = None) -> Path:
    t_start = time.time()
    cfg = cfg or load_site_config()
    channel = load_channel(channel_name)
    meas = load_measurement(channel, measurement)
    spec = model if isinstance(model, ModelSpec) else ModelSpec.load(model)
    default_slug = f"{spec.name}__{channel.name}" + (f"__{meas.name}" if meas.name != "published" else "")
    run_dir = unique_run_dir(out_root or cfg.runs, slug or default_slug)
    figs = ensure_dir(run_dir / "figs")
    warnings: list[str] = []
    timings: dict = {}
    dump_yaml_or_json(spec.to_dict(), run_dir / ("model.yaml" if spec.path and spec.path.suffix in (".yaml", ".yml") else "model.json"))
    dump_json(channel.to_dict(), run_dir / "channel.json")
    dump_json(meas.to_dict(), run_dir / "measurement.json")
    warnings += [f"measurement note: {n}" for n in meas.notes]

    # ---- 1. realise the model -------------------------------------------------------------
    t0 = time.time()
    ctx = RealizeContext(cfg, workdir=cfg.runs / "_generator_cache")
    pred: Prediction = realize(spec, channel, ctx)
    timings["realize_s"] = round(time.time() - t0, 1)
    pred_summary = {}
    if pred.truth is not None:
        pred_summary = pred.truth.summary()
        dump_json({"summary": pred_summary, "meta": pred.truth.meta}, run_dir / "truth_summary.json")
    warnings += pred.notes
    if pred.truth is not None and pred.truth.meta.get("has_geometry") is False and channel.phase_space.get("vertex"):
        warnings.append("model sample has no detector geometry: the fiducial-vertex phase-space cut was not applied "
                        "(the sample is taken as generated on the fiducial target; normalisation uses n_nucleons)")

    ctx_out = {"platform_version": __version__, "model": spec.to_dict(),
               "channel": {"name": channel.name, "description": channel.description},
               "measurement": meas.to_dict(),
               "prediction_summary": {k: v for k, v in pred_summary.items() if k != "norm"} | ({"normalisation": pred_summary.get("norm")} if pred_summary else {}),
               "provenance": pred.provenance, "figures": [], "warnings": warnings}
    fig_paths = []
    title_meas = f"{meas.x.axis_label('reco')} × {meas.y.axis_label('reco')}" if not meas.is_1d else meas.x.axis_label("reco")

    # ---- 2. folded-space comparison (forward folding: the primary mode) --------------------------
    if "folded" in modes and pred.truth is not None:
        t0 = time.time()
        try:
            sur, sur_path = resolve_surrogate(channel, meas, cfg, surrogate_path)
            if sur is None:
                warnings.append(f"no detector surrogate for measurement {meas.name!r} "
                                f"(run `ndp surrogate build --channel {channel.name} --measurement {meas.name}`)")
            else:
                data = fld.data_reco_cells(channel, meas, cfg)
                res = fld.compare_folded(channel, meas, pred.truth, sur, data, use_events=fold_events)
                res_out = {k: v for k, v in res.items() if not k.endswith("_cells")}
                res_out.update({"surrogate": str(sur_path.relative_to(cfg.repo_root)) if str(sur_path).startswith(str(cfg.repo_root)) else str(sur_path),
                                "surrogate_kind": sur.kind, "surrogate_meta": sur.meta, "data_sources": data["sources"]})
                ctx_out["folded"] = res_out
                np.savez(run_dir / "folded_cells.npz", **{k: v for k, v in res.items() if k.endswith("_cells")})
                fig_paths.append(plots.folded_projections(res, figs / "folded_projections.png",
                                                          f"{spec.name} → {sur.kind} → {channel.experiment} data ({title_meas}, folded space)"))
                if not meas.is_1d:
                    fig_paths.append(plots.cell_ratio_map(meas.binning, res["data_cells"], res["pred_cells"], figs / "folded_cell_ratio.png",
                                                          f"data / ({spec.name} → surrogate), per reco cell",
                                                          meas.x.axis_label("reco"), meas.y.axis_label("reco")))
        except Exception as e:  # keep the run alive; report the failure loudly
            warnings.append(f"folded comparison failed: {type(e).__name__}: {e}")
        timings["folded_s"] = round(time.time() - t0, 1)
    elif "folded" in modes:
        warnings.append("folded comparison skipped: the model has no truth events to smear (shipped curve)")

    # ---- 3. unfolded-space comparison (only where the experiment published this grid) -------------
    rel = None
    if "unfolded" in modes and meas.release:
        t0 = time.time()
        try:
            from .compare.minerva_bridge import PaperRelease
            rel = PaperRelease(cfg.require("minerva_repo"), meas.release)
            if pred.xsec_vector is not None:
                vec, var, how = pred.xsec_vector, None, "shipped curve (already d2sigma per nucleon)"
            else:
                x = unf.xsec_vector_from_truth(channel, meas, pred.truth)
                vec, var, how = x["vec"], x["var"], x["normalisation"]
                ctx_out["prediction_summary"]["sigma_total_phase_space_cm2_per_nucleon"] = x["sigma_total_phase_space_cm2"]
                ctx_out["prediction_summary"]["n_signal_in_phase_space"] = x["n_signal_in_ps"]
            sc = unf.score_unfolded(meas, rel, vec, var, spec.name)
            sc["ranking"] = unf.shipped_ranking(rel)
            sc["normalisation"] = how
            ctx_out["unfolded"] = sc
            np.save(run_dir / "unfolded_model_vector.npy", vec)
            cv_tune = next((m["name"] for m in rel.manifest.get("models", []) if m.get("cv_tune")), None)
            overlay = {spec.name: vec}
            if cv_tune:
                overlay[cv_tune] = rel.shipped_curve(cv_tune)
            fig_paths.append(plots.unfolded_projections(rel, overlay, figs / "unfolded_projections.png",
                                                        f"{spec.name} vs arXiv:{rel.arxiv} (unfolded space)"))
            fig_paths.append(plots.cell_ratio_map(meas.binning, vec, np.where(rel.mask, rel.data, 0.0), figs / "unfolded_cell_ratio.png",
                                                  f"{spec.name} / data, per cell (unfolded)",
                                                  meas.x.axis_label("truth"), meas.y.axis_label("truth")))
        except Exception as e:
            warnings.append(f"unfolded comparison failed: {type(e).__name__}: {e}")
        timings["unfolded_s"] = round(time.time() - t0, 1)
    elif "unfolded" in modes:
        warnings.append(f"unfolded comparison not available: measurement {meas.name!r} has no published release "
                        "(user-defined observables are compared in folded space only)")

    # ---- 4. report + manifest -----------------------------------------------------------------
    ctx_out["figures"] = [str(p.relative_to(run_dir)) for p in fig_paths]
    ctx_out["warnings"] = warnings
    report.write_report(run_dir, ctx_out)
    inputs = []
    for key in ("reco_data_files", "reco_mc_files"):
        for fn in channel.data.get(key, []):
            p = Path(cfg.data_dir or "") / fn
            if p.exists():
                inputs.append({"role": key, **cheap_fingerprint(p)})
    if pred.truth is not None and pred.truth.meta.get("source"):
        inputs.append({"role": "model_truth_source", "path": str(pred.truth.meta["source"])})
    manifest = {
        "run_id": run_dir.name, "timestamp": timestamp(), "platform_version": __version__,
        "platform_git": git_state(cfg.repo_root), "minerva_repo_git": git_state(cfg.minerva_repo) if cfg.minerva_repo else None,
        "model": spec.to_dict(), "model_fingerprint": spec.fingerprint(), "channel": channel.name,
        "channel_file": str(channel.path), "measurement": meas.name, "measurement_file": str(meas.path) if meas.path else None,
        "modes": list(modes), "inputs": inputs, "versions": versions(),
        "site_config": cfg.as_dict(), "timings_s": timings | {"total_s": round(time.time() - t_start, 1)},
        "outputs": ["scorecard.json", "report.md", "measurement.json", *ctx_out["figures"]], "warnings": warnings,
        "results_summary": _summary(ctx_out),
    }
    dump_json(manifest, run_dir / "manifest.json")
    return run_dir


def _summary(ctx: dict) -> dict:
    out = {"measurement": ctx.get("measurement", {}).get("name")}
    if ctx.get("folded"):
        g = ctx["folded"]["gof"]; t = ctx["folded"]["totals"]
        out["folded"] = {"minus2lnL_per_ndf": g["minus2lnL_per_ndf"], "pearson_per_ndf": g["pearson_per_ndf"],
                         "data_over_pred": t["ratio_data_over_pred"], "n_data": t["data"], "n_pred": t["pred"],
                         "surrogate": ctx["folded"].get("surrogate")}
    if ctx.get("unfolded"):
        r0 = ctx["unfolded"]["rows"][0]
        out["unfolded"] = {k: r0[k] for k in ("chi2_total_per_ndf", "chi2_shape_per_ndf", "alpha_shape", "norm_offset_pct", "n_populated")}
    return out
