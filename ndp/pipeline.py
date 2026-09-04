"""Orchestration: one theorist model + one channel -> a manifest-backed run directory.

    run_dir/
      model.yaml|json     the model spec as run
      channel.json        the channel manifest snapshot
      manifest.json       inputs (fingerprints), versions, git state, timings, results summary
      scorecard.json      every number produced
      report.md           human summary
      figs/               comparison figures
      truth_summary.json  the realised sample (counts, interaction mix, normalisation)
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
from .channels import load_channel, ChannelSpec
from .theory.models import ModelSpec, realize, RealizeContext, Prediction
from .surrogate.base import load_surrogate
from .compare import unfolded as unf, folded as fld, plots, report


def resolve_surrogate(channel: ChannelSpec, cfg: SiteConfig, explicit: str | None = None):
    cands = [explicit] if explicit else [channel.surrogate.get("default"), channel.surrogate.get("fallback")]
    for c in cands:
        if not c:
            continue
        p = Path(c)
        if not p.is_absolute():
            p = cfg.repo_root / p
        if (p / "surrogate.json").exists():
            return load_surrogate(p), p
    return None, None


def run_model(model: str | Path | ModelSpec, channel_name: str, *, cfg: SiteConfig | None = None,
              out_root: str | Path | None = None, surrogate_path: str | None = None,
              modes: tuple = ("unfolded", "folded"), fold_events: bool = False, slug: str | None = None) -> Path:
    t_start = time.time()
    cfg = cfg or load_site_config()
    channel = load_channel(channel_name)
    spec = model if isinstance(model, ModelSpec) else ModelSpec.load(model)
    run_dir = unique_run_dir(out_root or cfg.runs, slug or f"{spec.name}__{channel.name}")
    figs = ensure_dir(run_dir / "figs")
    warnings: list[str] = []
    timings: dict = {}
    dump_yaml_or_json(spec.to_dict(), run_dir / ("model.yaml" if spec.path and spec.path.suffix in (".yaml", ".yml") else "model.json"))
    dump_json(channel.to_dict(), run_dir / "channel.json")

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

    ctx_out = {"platform_version": __version__, "model": spec.to_dict(), "channel": {"name": channel.name, "description": channel.description},
               "prediction_summary": {k: v for k, v in pred_summary.items() if k != "norm"} | ({"normalisation": pred_summary.get("norm")} if pred_summary else {}),
               "provenance": pred.provenance, "figures": [], "warnings": warnings}
    fig_paths = []

    # ---- 2. unfolded-space comparison --------------------------------------------------------
    rel = None
    if "unfolded" in modes and channel.data.get("paper_manifest"):
        t0 = time.time()
        try:
            from .compare.minerva_bridge import PaperRelease
            rel = PaperRelease(cfg.require("minerva_repo"), channel.data["paper_manifest"])
            if pred.xsec_vector is not None:
                vec, var, how = pred.xsec_vector, None, "shipped curve (already d2sigma per nucleon)"
            else:
                x = unf.xsec_vector_from_truth(channel, pred.truth)
                vec, var, how = x["vec"], x["var"], x["normalisation"]
                ctx_out["prediction_summary"]["sigma_total_phase_space_cm2_per_nucleon"] = x["sigma_total_phase_space_cm2"]
                ctx_out["prediction_summary"]["n_signal_in_phase_space"] = x["n_signal_in_ps"]
            sc = unf.score_unfolded(channel, rel, vec, var, spec.name)
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
            fig_paths.append(plots.cell_ratio_map(channel.binning, vec, np.where(rel.mask, rel.data, 0.0), figs / "unfolded_cell_ratio.png",
                                                  f"{spec.name} / data, per cell (unfolded)"))
        except Exception as e:  # keep the run alive; report the failure loudly
            warnings.append(f"unfolded comparison failed: {type(e).__name__}: {e}")
        timings["unfolded_s"] = round(time.time() - t0, 1)

    # ---- 3. folded-space comparison ------------------------------------------------------------
    if "folded" in modes and pred.truth is not None:
        t0 = time.time()
        try:
            sur, sur_path = resolve_surrogate(channel, cfg, surrogate_path)
            if sur is None:
                warnings.append("no detector surrogate available for this channel (run `ndp surrogate build`)")
            else:
                data = fld.data_reco_cells(channel, cfg)
                res = fld.compare_folded(channel, pred.truth, sur, data, use_events=fold_events)
                res_out = {k: v for k, v in res.items() if not k.endswith("_cells")}
                res_out.update({"surrogate": str(sur_path.relative_to(cfg.repo_root)) if str(sur_path).startswith(str(cfg.repo_root)) else str(sur_path),
                                "surrogate_kind": sur.kind, "surrogate_meta": sur.meta})
                ctx_out["folded"] = res_out
                np.savez(run_dir / "folded_cells.npz", **{k: v for k, v in res.items() if k.endswith("_cells")})
                fig_paths.append(plots.folded_projections(res, figs / "folded_projections.png", f"{spec.name} → {sur.kind} → MINERvA data (folded space)"))
                fig_paths.append(plots.cell_ratio_map(channel.binning, res["data_cells"], res["pred_cells"], figs / "folded_cell_ratio.png",
                                                      f"data / ({spec.name} → surrogate), per reco cell"))
        except Exception as e:
            warnings.append(f"folded comparison failed: {type(e).__name__}: {e}")
        timings["folded_s"] = round(time.time() - t0, 1)
    elif "folded" in modes:
        warnings.append("folded comparison skipped: the model has no truth events to smear (shipped curve)")

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
        "channel_file": str(channel.path), "modes": list(modes), "inputs": inputs, "versions": versions(),
        "site_config": cfg.as_dict(), "timings_s": timings | {"total_s": round(time.time() - t_start, 1)},
        "outputs": ["scorecard.json", "report.md", *ctx_out["figures"]], "warnings": warnings,
        "results_summary": _summary(ctx_out),
    }
    dump_json(manifest, run_dir / "manifest.json")
    return run_dir


def _summary(ctx: dict) -> dict:
    out = {}
    if ctx.get("unfolded"):
        r0 = ctx["unfolded"]["rows"][0]
        out["unfolded"] = {k: r0[k] for k in ("chi2_total_per_ndf", "chi2_shape_per_ndf", "alpha_shape", "norm_offset_pct", "n_populated")}
    if ctx.get("folded"):
        g = ctx["folded"]["gof"]; t = ctx["folded"]["totals"]
        out["folded"] = {"minus2lnL_per_ndf": g["minus2lnL_per_ndf"], "pearson_per_ndf": g["pearson_per_ndf"],
                         "data_over_pred": t["ratio_data_over_pred"], "n_data": t["data"], "n_pred": t["pred"]}
    return out
