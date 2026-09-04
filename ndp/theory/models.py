"""Theorist-facing model specification and its realisation.

A model YAML has a `kind` and kind-specific fields:

  kind: shipped_curve     a generator curve shipped with the paper's data release
                          (unfolded comparison only)          fields: curve
  kind: reference_mc      the experiment's own truth-level MC as the model
                          (POT-normalised)                    fields: sample (default: channel MC)
  kind: reweight          a reference sample with per-event weights from a numpy expression
                          or a python function of the truth columns / observables
                          fields: base (a model spec or 'reference_mc'), weight_expr | weight_module
  kind: genie             run GENIE with the channel flux and target mix
                          fields: tune, generator_list, n_events, n_jobs, seed, gxmlpath, ...
  kind: external          a truth sample provided by the theorist
                          fields: path, format (ndp_npz | genie_gst), sigma_per_nucleon_cm2 (optional)

`realize(spec, channel, ctx)` returns a `Prediction`: a TruthTable (absolute or
POT-normalised) and/or a d2sigma vector in the paper's cell basis.
"""
from __future__ import annotations

import importlib.util
import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..events import TruthTable, Normalization, INT_CODE
from ..io import load_yaml_or_json, sha256_text
from ..channels import ChannelSpec
from ..channels import observables as obs

KINDS = ("shipped_curve", "reference_mc", "reweight", "genie", "external")


@dataclass
class ModelSpec:
    name: str
    kind: str
    params: dict = field(default_factory=dict)
    description: str = ""
    author: str = ""
    path: Path | None = None

    @staticmethod
    def load(path: str | Path) -> "ModelSpec":
        raw = load_yaml_or_json(path)
        return ModelSpec.from_dict(raw, Path(path))

    @staticmethod
    def from_dict(raw: dict, path: Path | None = None) -> "ModelSpec":
        if "kind" not in raw or raw["kind"] not in KINDS:
            raise ValueError(f"model needs kind in {KINDS}, got {raw.get('kind')!r}")
        params = {k: v for k, v in raw.items() if k not in ("name", "kind", "description", "author")}
        return ModelSpec(name=raw.get("name", path.stem if path else "model"), kind=raw["kind"], params=params,
                         description=raw.get("description", ""), author=raw.get("author", ""), path=path)

    def to_dict(self) -> dict:
        return {"name": self.name, "kind": self.kind, "description": self.description, "author": self.author, **self.params}

    def fingerprint(self) -> str:
        return sha256_text(json.dumps(self.to_dict(), sort_keys=True, default=str))[:16]

    def validate(self) -> list[str]:
        errs = []
        p = self.params
        if self.kind == "shipped_curve" and "curve" not in p:
            errs.append("shipped_curve needs `curve` (a model name in the paper release)")
        if self.kind == "reweight":
            if "weight_expr" not in p and "weight_module" not in p:
                errs.append("reweight needs `weight_expr` or `weight_module`")
            if "base" not in p:
                errs.append("reweight needs `base` (a model dict or 'reference_mc')")
        if self.kind == "external" and "path" not in p:
            errs.append("external needs `path`")
        if self.kind == "genie" and "tune" not in p:
            errs.append("genie needs `tune`")
        return errs


@dataclass
class Prediction:
    model: ModelSpec
    truth: TruthTable | None = None
    xsec_vector: np.ndarray | None = None       # d2sigma per cell in the channel basis [cm^2/GeV^2/nucleon]
    provenance: dict = field(default_factory=dict)
    notes: list = field(default_factory=list)


class RealizeContext:
    """What realisation needs from the site: config, cached reference MC, GENIE runner."""

    def __init__(self, site_cfg, workdir: Path):
        self.cfg = site_cfg
        self.workdir = Path(workdir)

    def reference_truth(self, channel: ChannelSpec) -> TruthTable:
        cache = self.cfg.require("data_dir") / "cache" / "truth_mc110040.npz"
        if cache.exists():
            return TruthTable.load(cache)
        from ..adapters.minerva_anatuple import read_truth
        files = channel.data.get("reco_mc_files", [])
        if not files:
            raise FileNotFoundError("channel lists no reco_mc_files for a reference MC")
        t = read_truth(self.cfg.require("data_dir") / files[0])
        cache.parent.mkdir(parents=True, exist_ok=True)
        t.save(cache)
        return t


# --------------------------------------------------------------------------------------
def realize(spec: ModelSpec, channel: ChannelSpec, ctx: RealizeContext) -> Prediction:
    errs = spec.validate()
    if errs:
        raise ValueError("; ".join(errs))
    fn = {"shipped_curve": _shipped, "reference_mc": _reference, "reweight": _reweight,
          "genie": _genie, "external": _external}[spec.kind]
    return fn(spec, channel, ctx)


def _shipped(spec, channel, ctx) -> Prediction:
    from ..compare.minerva_bridge import PaperRelease
    rel = PaperRelease(ctx.cfg.require("minerva_repo"), channel.data["paper_manifest"])
    vec = rel.shipped_curve(spec.params["curve"])
    return Prediction(spec, xsec_vector=vec, provenance={"release": rel.arxiv, "curve": spec.params["curve"]},
                      notes=["shipped curve: unfolded comparison only (no truth events to fold)"])


def _reference(spec, channel, ctx) -> Prediction:
    t = ctx.reference_truth(channel)
    return Prediction(spec, truth=t, provenance={"source": t.meta.get("source"), "norm": t.meta.get("norm")})


def _reweight(spec, channel, ctx) -> Prediction:
    base = spec.params["base"]
    base_spec = ModelSpec.from_dict({"name": f"{spec.name}__base", **base}) if isinstance(base, dict) \
        else ModelSpec(name=f"{spec.name}__base", kind=str(base))
    base_pred = realize(base_spec, channel, ctx)
    if base_pred.truth is None:
        raise ValueError("reweight base must produce truth events")
    t = base_pred.truth
    w = _weights(spec.params, t, channel)
    if not np.all(np.isfinite(w)) or np.any(w < 0):
        raise ValueError("weight function produced negative or non-finite weights")
    t2 = t.with_weights(t["weight"] * w, note=spec.params.get("weight_expr") or spec.params.get("weight_module"))
    prov = {"base": base_pred.provenance, "weight": spec.params.get("weight_expr") or spec.params.get("weight_module"),
            "mean_weight": float(np.mean(w)), "weight_range": [float(w.min()), float(w.max())]}
    return Prediction(spec, truth=t2, provenance=prov)


def _weights(params: dict, t: TruthTable, channel: ChannelSpec) -> np.ndarray:
    if "weight_module" in params:
        mod_path, _, func = str(params["weight_module"]).partition(":")
        p = Path(mod_path)
        if not p.is_absolute() and channel.path is not None:
            for base in (Path.cwd(), channel.path.parents[1]):
                if (base / p).exists():
                    p = base / p
                    break
        spec_ = importlib.util.spec_from_file_location("ndp_user_weight", p)
        mod = importlib.util.module_from_spec(spec_)
        spec_.loader.exec_module(mod)
        return np.asarray(getattr(mod, func or "weight")(t), float)
    expr = params["weight_expr"]
    ns = {k: t[k] for k in t.columns if k not in ("fs_offsets",) and not k.startswith("fs_")}
    frame = channel.phase_space.get("frame", "detector")
    for name in obs.OBSERVABLES:
        try:
            ns[name] = obs.evaluate(name, t, frame=frame)
        except Exception:
            pass
    ns.update({k: v for k, v in INT_CODE.items()})          # QE, RES, DIS, COH, MEC as codes
    ns.update({"np": np, "where": np.where, "exp": np.exp, "log": np.log, "sqrt": np.sqrt, "abs": np.abs,
               "clip": np.clip, "minimum": np.minimum, "maximum": np.maximum, "pi": np.pi})
    w = eval(expr, {"__builtins__": {}}, ns)  # noqa: S307 — theorist-authored expression, numpy namespace only
    return np.broadcast_to(np.asarray(w, float), (t.n,)).copy()


def _genie(spec, channel, ctx) -> Prediction:
    from .generator import GenieSpec, generate
    from . import flux as fluxmod
    p = dict(spec.params)
    mix = p.pop("target_mix", None) or channel.normalization.get("target_mix_mass_fractions")
    fl = fluxmod.load_channel_flux(channel, ctx.cfg.repo_root)
    g = GenieSpec(tune=p.get("tune", "G18_02a_00_000"), generator_list=p.get("generator_list", "CC"),
                  n_events=int(p.get("n_events", 20000)), nu_pdg=int(p.get("nu_pdg", channel.signal["nu_pdg"][0])),
                  target_mix={int(k): float(v) for k, v in mix.items()}, e_min=float(p.get("e_min", 0.0)),
                  e_max=float(p.get("e_max", 100.0)), seed=int(p.get("seed", 1)), run_number=int(p.get("run_number", 1)),
                  n_jobs=int(p.get("n_jobs", 1)), env_json=p.get("env_json"), splines=p.get("splines"),
                  gxmlpath=p.get("gxmlpath"), extra_args=list(p.get("extra_args", [])))
    t = generate(g, fl["edges"], fl["density_cm2_pot_gev"], flux_source=str(channel.normalization["flux_table"]),
                 workdir=ctx.workdir, site_cfg=ctx.cfg)
    return Prediction(spec, truth=t, provenance={"genie": t.meta.get("genie_spec"), "source": t.meta.get("source"),
                                                 "sigma_flux_avg_per_nucleon_cm2": t.meta.get("sigma_flux_avg_per_nucleon_cm2"),
                                                 "cache_hit": t.meta.get("cache_hit", False)})


def _external(spec, channel, ctx) -> Prediction:
    p = spec.params
    path = Path(p["path"])
    fmt = p.get("format", "ndp_npz")
    if fmt == "ndp_npz":
        t = TruthTable.load(path)
    elif fmt == "genie_gst":
        from ..adapters.genie_gst import read_gst
        t = read_gst(path)
    else:
        raise ValueError(f"unknown external format {fmt!r}")
    if "sigma_per_nucleon_cm2" in p:
        t.meta["norm"] = Normalization(kind="xsec_per_nucleon",
                                       xsec_per_unit_weight=float(p["sigma_per_nucleon_cm2"]) / float(t["weight"].sum()),
                                       notes="sigma_per_nucleon_cm2 supplied in the model spec").to_dict()
    return Prediction(spec, truth=t, provenance={"path": str(path), "format": fmt, "norm": t.meta.get("norm")})
