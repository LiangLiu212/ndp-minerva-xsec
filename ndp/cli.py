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
    from .channels import list_channels, load_channel
    for name in list_channels():
        ch = load_channel(name)
        status = ch.raw.get("status", "ready")
        print(f"{name:40s} [{status}] {ch.experiment}: {ch.description.strip()[:90]}")
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
    from .channels import load_channel
    m = ModelSpec.load(a.model)
    errs = m.validate()
    ch = load_channel(a.channel) if a.channel else None
    print(json.dumps({"model": m.to_dict(), "errors": errs, "channel": ch.name if ch else None}, indent=2, default=str))
    return 1 if errs else 0


def _cmd_run(a):
    from .pipeline import run_model
    modes = tuple(a.modes.split(","))
    rd = run_model(a.model, a.channel, out_root=a.out, surrogate_path=a.surrogate, modes=modes, fold_events=a.fold_events, slug=a.slug)
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
    from .channels import load_channel
    from .surrogate.binned import BinnedResponse
    from .surrogate.parametric import SmearingSurrogate
    from .events import TruthTable
    from .compare import plots
    from .io import cheap_fingerprint, timestamp
    cfg = load_site_config(); ch = load_channel(a.channel)
    out_root = cfg.surrogates / ch.name
    b = ch.binning
    if a.source == "artifacts":
        run = cfg.require("minerva_repo") / "runs" / "2026-06-19_me_inclusive_ddxsec" / "bench"
        sur = BinnedResponse.from_run_artifacts(b, run / "step6_migration.npy", run / "step7_efficiency_2d.npy",
                                                meta={"channel": ch.name, "built": timestamp(),
                                                      "source": "ndp-minerva-data-release-exploration run 2026-06-19_me_inclusive_ddxsec (audited, manifest-backed)"})
        out = out_root / "binned_from_run_2026-06-19"
        sur.save(out); plots.response_figure(sur, out / "response.png", "binned response from run artifacts")
        print(json.dumps(sur.diagnostics(), indent=2)); print("saved", out)
        return 0
    # from the MINERvA MC: cached tables written by the adapter (see README "data cache")
    cache = cfg.require("data_dir") / "cache"
    truth = TruthTable.load(cache / "truth_mc110040.npz")
    z = np.load(cache / "reco_mc110040.npz"); rt = TruthTable.load(cache / "reco_mc110040_truthcols.npz")
    passed = z["passed"]; pT_r, pz_r = z["reco_pT"], z["reco_pz"]
    sig_den = ch.is_signal(truth) & ch.in_phase_space(truth)
    xd, yd = ch.observables(truth)
    sig_num = passed & ch.is_signal(rt) & ch.in_phase_space(rt)
    xn, yn = ch.observables(rt)
    bkg = passed & ~(ch.is_signal(rt) & ch.in_phase_space(rt))     # non-signal + out-of-phase-space signal
    pot_mc = truth.norm.pot
    meta = {"channel": ch.name, "built": timestamp(), "training_mc": cheap_fingerprint(cfg.data_dir / ch.data["reco_mc_files"][0]),
            "pot_mc": pot_mc, "generator": truth.meta.get("generator"), "selection": ch.selection.get("name"),
            "phase_space": ch.phase_space, "signal": ch.signal}
    kw = dict(x_true_den=xd[sig_den], y_true_den=yd[sig_den], x_true_num=xn[sig_num], y_true_num=yn[sig_num],
              x_reco_num=pT_r[sig_num], y_reco_num=pz_r[sig_num])
    if a.kind in ("binned", "all"):
        sur = BinnedResponse.fit(b, **kw, x_reco_bkg=pT_r[bkg], y_reco_bkg=pz_r[bkg], pot_mc=pot_mc, meta=meta)
        out = out_root / "binned_mc110040"; sur.save(out); plots.response_figure(sur, out / "response.png", "binned response, MC run 110040")
        print("binned:", json.dumps(sur.diagnostics(), indent=2)); print("saved", out)
    if a.kind in ("parametric", "all"):
        sur = SmearingSurrogate.fit(b, **kw, mode="diff,ratio", n_samples=a.n_samples, meta=meta)
        out = out_root / "parametric_mc110040"; sur.save(out); plots.response_figure(sur, out / "response.png", "parametric smearing, MC run 110040")
        print("parametric:", json.dumps(sur.diagnostics(), indent=2)); print("saved", out)
    return 0


def _cmd_surrogate_inspect(a):
    from .surrogate.base import load_surrogate
    s = load_surrogate(a.path)
    print(json.dumps({"kind": s.kind, "diagnostics": s.diagnostics(), "meta": s.meta}, indent=2, default=str))
    return 0


def _cmd_data_status(a):
    from .channels import load_channel
    cfg = load_site_config()
    print("site config:", json.dumps(cfg.as_dict(), indent=2))
    for name in ([a.channel] if a.channel else __import__("ndp.channels", fromlist=["list_channels"]).list_channels()):
        ch = load_channel(name)
        for key in ("reco_data_files", "reco_mc_files"):
            for fn in ch.data.get(key, []):
                p = (cfg.data_dir or Path(".")) / fn
                print(f"{name}: {key} {fn} -> {'present' if p.exists() else 'MISSING'}")
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
    sub.add_parser("channels", help="list channel manifests").set_defaults(fn=_cmd_channels)
    sub.add_parser("models", help="list + validate model specs in models/").set_defaults(fn=_cmd_models)
    p = sub.add_parser("validate", help="validate a model spec"); p.add_argument("model"); p.add_argument("--channel"); p.set_defaults(fn=_cmd_validate)
    p = sub.add_parser("run", help="test a model against a channel"); p.add_argument("model"); p.add_argument("--channel", required=True)
    p.add_argument("--out"); p.add_argument("--surrogate"); p.add_argument("--modes", default="unfolded,folded")
    p.add_argument("--fold-events", action="store_true", help="smear truth events (parametric surrogate) instead of folding true cells")
    p.add_argument("--slug"); p.set_defaults(fn=_cmd_run)
    ps = sub.add_parser("surrogate", help="build / inspect detector surrogates").add_subparsers(dest="scmd", required=True)
    p = ps.add_parser("build"); p.add_argument("--channel", required=True); p.add_argument("--source", choices=("mc", "artifacts"), default="mc")
    p.add_argument("--kind", choices=("binned", "parametric", "all"), default="all"); p.add_argument("--n-samples", type=int, default=20); p.set_defaults(fn=_cmd_surrogate_build)
    p = ps.add_parser("inspect"); p.add_argument("path"); p.set_defaults(fn=_cmd_surrogate_inspect)
    p = sub.add_parser("data", help="data availability").add_subparsers(dest="dcmd", required=True).add_parser("status"); p.add_argument("--channel"); p.set_defaults(fn=_cmd_data_status)
    p = sub.add_parser("flux", help="channel flux table summary / export"); p.add_argument("--channel", required=True); p.add_argument("--out"); p.set_defaults(fn=_cmd_flux)
    a = ap.parse_args(argv)
    return a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
