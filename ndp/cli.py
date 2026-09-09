"""Command line: `python -m ndp <command>` (or `ndp` once installed)."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from . import __version__
from .config import load_site_config
from .io import dump_json


def _cmd_channels(a):
    from .channels import list_channels, load_channel, list_measurements
    for name in list_channels():
        ch = load_channel(name)
        status = ch.raw.get("status", "ready")
        print(f"{name:40s} [{status}] {ch.experiment}: {ch.description.strip()[:90]}")
        try:
            print(f"{'':40s}   measurements: {', '.join(list_measurements(ch))}")
        except Exception as e:
            print(f"{'':40s}   measurements: ERROR {e}")
    return 0


def _cmd_measurements(a):
    from .channels import load_channel, list_measurements, load_measurement
    from .pipeline import resolve_surrogate
    cfg = load_site_config(); ch = load_channel(a.channel)
    for name in list_measurements(ch):
        try:
            m = load_measurement(ch, name)
        except Exception as e:
            print(f"{name:28s} ERROR {e}"); continue
        sur, sp = resolve_surrogate(ch, m, cfg)
        grid = f"{m.x.truth}↔{m.x.reco} [{len(m.x.edges) - 1}]" + ("" if m.is_1d else f" × {m.y.truth}↔{m.y.reco} [{len(m.y.edges) - 1}]")
        print(f"{name:28s} [{m.status}] {grid}; release={m.release or '—'}; surrogate={'none' if sur is None else sp.relative_to(cfg.repo_root)}")
        if m.description:
            print(f"{'':28s} {m.description.strip()[:110]}")
    return 0


def _cmd_models(a):
    from .theory.models import ModelSpec
    cfg = load_site_config()
    for p in sorted((cfg.repo_root / "models").glob("*.y*ml")):
        try:
            m = ModelSpec.load(p)
            errs = m.validate()
            print(f"{p.name:40s} kind={m.kind:14s} {'OK' if not errs else 'INVALID: ' + '; '.join(errs)}  {m.description.strip()[:70]}")
        except Exception as e:
            print(f"{p.name:40s} ERROR {e}")
    return 0


def _cmd_validate(a):
    from .theory.models import ModelSpec
    from .channels import load_channel, load_measurement
    m = ModelSpec.load(a.model)
    errs = m.validate()
    ch = load_channel(a.channel) if a.channel else None
    meas = load_measurement(ch, a.measurement) if ch else None
    print(json.dumps({"model": m.to_dict(), "errors": errs, "channel": ch.name if ch else None,
                      "measurement": meas.to_dict() if meas else None}, indent=2, default=str))
    return 1 if errs else 0


def _cmd_run(a):
    from .pipeline import run_model
    modes = tuple(a.modes.split(","))
    rd = run_model(a.model, a.channel, measurement=a.measurement, out_root=a.out, surrogate_path=a.surrogate, modes=modes,
                   fold_events=a.fold_events, slug=a.slug)
    man = json.loads((rd / "manifest.json").read_text())
    print(f"run dir: {rd}")
    print(json.dumps(man["results_summary"], indent=2))
    if man["warnings"]:
        print("warnings:")
        for w in man["warnings"]:
            print("  -", w)
    print((rd / "report.md").read_text())
    return 0


def _cmd_surrogate_build(a):
    from .channels import load_channel, load_measurement
    from .surrogate.binned import BinnedResponse
    from .surrogate.build import build_surrogates
    from .compare import plots
    from .io import timestamp
    cfg = load_site_config(); ch = load_channel(a.channel)
    meas = load_measurement(ch, a.measurement)
    if a.source == "artifacts":
        if meas.name != "published":
            print("artifacts source only exists for the published grid"); return 2
        run = cfg.require("minerva_repo") / "runs" / "2026-06-19_me_inclusive_ddxsec" / "bench"
        sur = BinnedResponse.from_run_artifacts(meas.binning, run / "step6_migration.npy", run / "step7_efficiency_2d.npy",
                                                meta={"channel": ch.name, "built": timestamp(),
                                                      "source": "ndp-minerva-data-release-exploration run 2026-06-19_me_inclusive_ddxsec (audited, manifest-backed)"})
        out = cfg.surrogates / ch.name / "binned_from_run_2026-06-19"
        sur.save(out); plots.response_figure(sur, out / "response.png", "binned response from run artifacts")
        print(json.dumps(sur.diagnostics(), indent=2)); print("saved", out)
        return 0
    kinds = ("binned", "parametric") if a.kind == "all" else (a.kind,)
    results = build_surrogates(ch, meas, cfg, kinds=kinds, n_samples=a.n_samples, out_root=a.out)
    for r in results:
        print(json.dumps(r, indent=2, default=str))
    bad = [r for r in results if r["kind"] == "binned" and not r["closure"]["exact"]]
    return 1 if bad else 0


def _cmd_surrogate_inspect(a):
    from .surrogate.base import load_surrogate
    s = load_surrogate(a.path)
    print(json.dumps({"kind": s.kind, "diagnostics": s.diagnostics(), "meta": s.meta}, indent=2, default=str))
    return 0


def _cmd_data_status(a):
    from .channels import load_channel
    from .adapters.minerva_anatuple import cache_tag, RECO_CACHE_VERSION
    cfg = load_site_config()
    print("site config:", json.dumps(cfg.as_dict(), indent=2))
    for name in ([a.channel] if a.channel else __import__("ndp.channels", fromlist=["list_channels"]).list_channels()):
        ch = load_channel(name)
        for key in ("reco_data_files", "reco_mc_files"):
            for fn in ch.data.get(key, []):
                p = (cfg.data_dir or Path(".")) / fn
                line = f"{name}: {key} {fn} -> {'present' if p.exists() else 'MISSING'}"
                if p.exists():
                    try:
                        tag = cache_tag(fn); cache = p.parent / "cache"
                        parts = []
                        for stem in ([f"truth_{tag}", f"reco_{tag}", f"reco_{tag}_truthcols"] if key == "reco_mc_files" else [f"reco_{tag}"]):
                            f = cache / f"{stem}.npz"
                            if not f.exists():
                                parts.append(f"{stem}: missing")
                            elif stem.startswith("reco_") and not stem.endswith("_truthcols"):
                                z = np.load(f, allow_pickle=False)
                                v = json.loads(str(z["__meta__"])).get("cache_version", 1) if "__meta__" in z.files else 1
                                parts.append(f"{stem}: v{v}" + ("" if v >= RECO_CACHE_VERSION else f" (stale, need v{RECO_CACHE_VERSION})"))
                            else:
                                parts.append(f"{stem}: ok")
                        line += "; cache " + ", ".join(parts)
                    except ValueError:
                        pass
                print(line)
    return 0


def _cmd_data_cache(a):
    from .channels import load_channel
    from .adapters.minerva_anatuple import build_cache
    cfg = load_site_config(); ch = load_channel(a.channel)
    data_dir = cfg.require("data_dir"); cache = data_dir / "cache"
    which = set(a.which.split(","))
    for key, is_mc in (("reco_data_files", False), ("reco_mc_files", True)):
        if ("mc" if is_mc else "data") not in which:
            continue
        for fn in ch.data.get(key, []):
            r = build_cache(data_dir / fn, cache, is_mc=is_mc, truth=not a.reco_only, entry_stop=a.entry_stop)
            print(json.dumps(r, indent=2))
    return 0


def _cmd_flux(a):
    from .channels import load_channel
    from .theory import flux
    cfg = load_site_config(); ch = load_channel(a.channel)
    fl = flux.load_channel_flux(ch, cfg.repo_root)
    print(json.dumps({k: v for k, v in fl.items() if k not in ("edges", "density_cm2_pot_gev")}, indent=2))
    print("channel phi_per_pot_cm2:", ch.normalization.get("phi_per_pot_cm2"))
    if a.out:
        flux.write_th1_root(a.out, fl["edges"], fl["density_cm2_pot_gev"]); print("wrote", a.out)
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="ndp", description=f"Neutrino Discovery Platform {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("channels", help="list channel manifests (and their measurements)").set_defaults(fn=_cmd_channels)
    p = sub.add_parser("measurements", help="list a channel's measurements (published grid + user-defined observables)")
    p.add_argument("--channel", required=True); p.set_defaults(fn=_cmd_measurements)
    sub.add_parser("models", help="list + validate model specs in models/").set_defaults(fn=_cmd_models)
    p = sub.add_parser("validate", help="validate a model spec"); p.add_argument("model"); p.add_argument("--channel"); p.add_argument("--measurement"); p.set_defaults(fn=_cmd_validate)
    p = sub.add_parser("run", help="test a model against a channel (forward-folded; unfolded too where published)")
    p.add_argument("model"); p.add_argument("--channel", required=True)
    p.add_argument("--measurement", help="measurement name or YAML path (default: the channel's published grid)")
    p.add_argument("--out"); p.add_argument("--surrogate"); p.add_argument("--modes", default="folded,unfolded")
    p.add_argument("--fold-events", action="store_true", help="smear truth events (parametric surrogate) instead of folding true cells")
    p.add_argument("--slug"); p.set_defaults(fn=_cmd_run)
    ps = sub.add_parser("surrogate", help="build / inspect detector surrogates").add_subparsers(dest="scmd", required=True)
    p = ps.add_parser("build"); p.add_argument("--channel", required=True); p.add_argument("--measurement")
    p.add_argument("--source", choices=("mc", "artifacts"), default="mc")
    p.add_argument("--kind", choices=("binned", "parametric", "all"), default="all"); p.add_argument("--n-samples", type=int, default=20)
    p.add_argument("--out", help="output root (default surrogates/<channel>[/<measurement>])"); p.set_defaults(fn=_cmd_surrogate_build)
    p = ps.add_parser("inspect"); p.add_argument("path"); p.set_defaults(fn=_cmd_surrogate_inspect)
    pd = sub.add_parser("data", help="data availability / caches").add_subparsers(dest="dcmd", required=True)
    p = pd.add_parser("status"); p.add_argument("--channel"); p.set_defaults(fn=_cmd_data_status)
    p = pd.add_parser("cache", help="build the truth/reco .npz caches from the channel's AnaTuples")
    p.add_argument("--channel", required=True); p.add_argument("--which", default="data,mc")
    p.add_argument("--reco-only", action="store_true", help="skip the (slow) Truth tree"); p.add_argument("--entry-stop", type=int)
    p.set_defaults(fn=_cmd_data_cache)
    p = sub.add_parser("flux", help="channel flux table summary / export"); p.add_argument("--channel", required=True); p.add_argument("--out"); p.set_defaults(fn=_cmd_flux)
    a = ap.parse_args(argv)
    return a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
