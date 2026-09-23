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
    from .pipeline import run_model, run_model_multi
    modes = tuple(a.modes.split(","))
    want = (a.measurement or "").strip()
    if want == "all" or "," in want:
        names = None if want == "all" else [n.strip() for n in want.split(",")]
        rd = run_model_multi(a.model, a.channel, measurements=names, out_root=a.out, efficiency_run=a.efficiency_run, slug=a.slug,
                             surrogate_kind=a.surrogate_kind, vbll_opts=_vbll_opts(a))
        man = json.loads((rd / "manifest.json").read_text())
        print(f"run dir: {rd}")
        if man["warnings"]:
            print("warnings:")
            for w in man["warnings"]:
                print("  -", w)
        print((rd / "report.md").read_text())
        return 0
    rd = run_model(a.model, a.channel, measurement=a.measurement, out_root=a.out, surrogate_path=a.surrogate, modes=modes,
                   fold_events=a.fold_events, slug=a.slug, surrogate_kind=a.surrogate_kind, vbll_opts=_vbll_opts(a))
    man = json.loads((rd / "manifest.json").read_text())
    print(f"run dir: {rd}")
    print(json.dumps(man["results_summary"], indent=2))
    if man["warnings"]:
        print("warnings:")
        for w in man["warnings"]:
            print("  -", w)
    print((rd / "report.md").read_text())
    return 0


def _vbll_opts(a) -> dict:
    """`--vbll-*` options of `run` / `surrogate build` -> VBLLEventSurrogate.from_parts keyword arguments."""
    return {"model_name": a.vbll_model, "n_samples": a.vbll_samples, "seed": a.vbll_seed, "truncate_to_reco_windows": not a.no_vbll_truncate}


def _add_vbll_args(p):
    p.add_argument("--vbll-model", default="x60_het", help="ported VBLL model under surrogates/<channel>/_vbll/ (default x60_het)")
    p.add_argument("--vbll-samples", type=int, default=20, help="smeared copies per truth event (default 20)")
    p.add_argument("--vbll-seed", type=int, default=0)
    p.add_argument("--no-vbll-truncate", action="store_true", help="do not condition the smeared reco on the selection's kinematic windows")


def _cmd_surrogate_build(a):
    from .channels import load_channel, load_measurement
    from .surrogate.binned import BinnedResponse
    from .surrogate.build import build_surrogates
    from .compare import plots
    from .io import timestamp
    cfg = load_site_config(); ch = load_channel(a.channel)
    meas = load_measurement(ch, a.measurement) if (a.measurement or "").strip() != "all" and "," not in (a.measurement or "") else None
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
    if a.kind == "vbll":
        # wrap each grid's binned response with a ported VBLL model (event-level migration) + closure on the MC truth
        from .channels import list_measurements
        from .surrogate.vbll import build_vbll
        want = (a.measurement or "").strip()
        names = [n for n in list_measurements(ch) if n != "published"] if want == "all" else ([n.strip() for n in want.split(",")] if want else [None])
        ms = [load_measurement(ch, n) for n in names]
        res = build_vbll(ch, ms, cfg, out_root=a.out, closure=not a.no_closure, **_vbll_opts(a))
        slim = [{k: v for k, v in r.items() if k != "closure"} |
                ({"closure": {kk: vv for kk, vv in r["closure"].items() if kk not in ("pred_cells", "reco_cells", "binned_fold_cells", "inputs")}} if "closure" in r else {})
                for r in res]
        print(json.dumps(slim, indent=2, default=str))
        return 0
    kinds = ("binned", "parametric") if a.kind == "all" else (a.kind,)
    from .channels import list_measurements
    from .products import has_products
    want = (a.measurement or "").strip()
    several = want == "all" or "," in want
    if several or has_products(ch) or a.chunked:
        # one pass over the MC, one playlist (or legacy file) at a time; binned responses only
        from .surrogate.chunked import build_surrogates_chunked
        if "parametric" in kinds and a.kind != "all":
            print("parametric surrogates need the event arrays in memory: not available with playlist products / --chunked"); return 2
        names = [n for n in list_measurements(ch) if n != "published"] if want == "all" else [n.strip() for n in want.split(",")]
        ms = [load_measurement(ch, n) for n in names]
        res = build_surrogates_chunked(ch, ms, cfg, out_root=a.out)
        print(json.dumps(res, indent=2, default=str))
        bad = [n for n, r in res.items() if not r["closure"]["exact"]]
        return 1 if bad else 0
    results = build_surrogates(ch, meas, cfg, kinds=kinds, n_samples=a.n_samples, out_root=a.out)
    for r in results:
        print(json.dumps(r, indent=2, default=str))
    bad = [r for r in results if r["kind"] == "binned" and not r["closure"]["exact"]]
    return 1 if bad else 0


def _cmd_surrogate_train_vbll(a):
    """Train a VBLL surrogate on the channel's selected signal pairs (ml environment) -> surrogates/<channel>/_vbll/<name>/."""
    from .channels import load_channel
    from .surrogate.vbll_train import collect_pairs, save_pairs, load_pairs, train_vbll, save_model
    cfg = load_site_config(); ch = load_channel(a.channel)
    pairs_path = Path(a.pairs) if a.pairs else cfg.require("data_dir") / "cache" / f"vbll_pairs_{ch.name}_{a.frame}.npz"
    if pairs_path.exists() and not a.recollect:
        pairs = load_pairs(pairs_path); print(f"pairs: {pairs['n_pairs']} from {pairs_path} (frame {pairs['frame']})")
    else:
        print(f"collecting pairs (frame {a.frame}) ...")
        pairs = collect_pairs(cfg, ch, frame=a.frame); save_pairs(pairs, pairs_path); print(f"pairs: {pairs['n_pairs']} -> {pairs_path}")
    inputs = tuple(a.inputs.split(",")); outputs = tuple(a.outputs.split(","))
    spec, state, info = train_vbll(pairs, inputs=inputs, outputs=outputs, head_type=a.head, d_embed=a.d_embed, hidden=a.hidden, n_layers=a.n_layers,
                                   noise_prior_scale=a.noise_prior_scale, lr=a.lr, epochs=a.epochs, patience=a.patience, batch_size=a.batch_size,
                                   val_fraction=a.val_fraction, seed=a.seed, threads=a.threads)
    out = save_model(cfg.surrogates / ch.name / "_vbll" / a.name, spec, state, info, pairs_path)
    print(json.dumps({"saved": str(out), "epochs_run": info["epochs_run"], "best_val_nll": info["best_val_nll"], "train_s": info["train_s"],
                      "validation_metrics": info["validation_metrics"]}, indent=1, default=str))
    return 0


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
    if a.url:                                   # one file (local path or root:// URL), grid-style outputs
        from .grid.process_file import process
        out = a.out or str(cfg.require("data_dir") / "cache")
        r = process(a.url, a.channel, a.kind or ("mc" if "_mc_" in a.url else "data"), out, skim=a.skim, entry_stop=a.entry_stop)
        print(json.dumps({k: v for k, v in r.items() if k != "fingerprint"}, indent=2, default=str))
        return 5 if r.get("status") == "failed" else 0
    data_dir = cfg.require("data_dir"); cache = data_dir / "cache"
    which = set(a.which.split(","))
    for key, is_mc in (("reco_data_files", False), ("reco_mc_files", True)):
        if ("mc" if is_mc else "data") not in which:
            continue
        for fn in ch.data.get(key, []):
            r = build_cache(data_dir / fn, cache, is_mc=is_mc, truth=not a.reco_only, entry_stop=a.entry_stop,
                            channel=ch if (a.skim or a.sidecar) else None, skim=a.skim, sidecar=a.sidecar)
            print(json.dumps({k: v for k, v in r.items() if k != "fingerprint"}, indent=2, default=str))
    return 0


def _cmd_data_merge(a):
    from .products import merge_playlist
    cfg = load_site_config()
    root = Path(a.products_dir) if a.products_dir else cfg.require("data_dir") / "products"
    for kind in a.kind.split(","):
        merge_playlist(root, a.beam, a.playlist, kind)
    return 0


def _cmd_grid(a):
    from .grid import campaign as cp
    if a.gcmd == "harvest-pot":
        cp.harvest_pot(workers=a.workers); return 0
    if a.gcmd == "plan":
        fpp = {}
        if a.files_per_process:
            for kv in a.files_per_process.split(","):
                k, v = kv.split("="); fpp[k] = int(v)
        cp.plan(a.name, a.channel, beams=tuple(a.beams.split(",")), kinds=tuple(a.kinds.split(",")),
                playlists=a.playlists.split(",") if a.playlists else None, files_per_process=fpp, pnfs_base=a.pnfs_base)
        return 0
    if a.gcmd == "status":
        cp.status(a.name, pot_check=not a.no_pot_check); return 0
    if a.gcmd == "resubmit":
        cp.resubmit(a.name); return 0
    if a.gcmd == "harvest":
        cfg = load_site_config()
        root = Path(a.products_dir) if a.products_dir else cfg.require("data_dir") / "products"
        print(json.dumps(cp.harvest(a.name, root, playlists=a.playlists.split(",") if a.playlists else None,
                                    workers=a.workers, archive=not a.no_archive)))
        return 0
    if a.gcmd == "stage-worklists":
        cp.stage_worklists(a.name, files=a.files or None); return 0
    if a.gcmd == "submit-cmd":
        c = cp.load_campaign(a.name)
        w = c["worklists"][a.worklist]
        wl_pnfs = c.get("staged_worklists", {}).get(a.file or w["file"]) or w.get("pnfs_worklist")
        if not wl_pnfs:
            raise SystemExit(f"worklist not staged on PNFS: run `ndp grid stage-worklists {a.name}` first")
        # requests: MaxRSS 1.26 GB measured on run 110040; the held smoke job showed ~1.45 GB charged during the cold
        # CVMFS import of the environment, so data jobs get 3000 MB and MC jobs 4000 MB
        # disk: 8 files x ~300 MB outputs per MC process, and playlists with more POT per file (1N: 1.65e19) exceed 4 GB
        req = {"mc": f"--memory 4000MB --disk {a.disk or '8GB'} --expected-lifetime 3h", "data": f"--memory 3000MB --disk {a.disk or '2GB'} --expected-lifetime 2h"}[w["kind"]]
        n = min(w["n_processes"], a.max_processes) if a.max_processes else w["n_processes"]
        print(f"python3 .claude/skills/jobsub-lite/scripts/jobsub.py submit --worker grid/worker.sh -N {n} --tar-label {a.tar_label} {req} "
              f"--jobsub-arg=--onsite -f {wl_pnfs} --pnfs-out {w['pnfs_out']} --runtype ndpstream --stem {a.stem or a.worklist} "
              f"-- -R @TAR_DIR@ -O @PNFS_OUT@ -W {Path(wl_pnfs).name} -K {w['kind']} -C {c['channel']} -n {w['files_per_process']}")
        return 0
    raise SystemExit(f"unknown grid command {a.gcmd}")


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


def _cmd_signal(a):
    from .diagnostics import run_signal_diagnostics
    run_dir = run_signal_diagnostics(a.channel, load_site_config(), cache=a.cache, out_root=a.out, slug=a.slug)
    print((run_dir / "report.md").read_text())
    print(f"run directory: {run_dir}")
    return 0


def _cmd_efficiency(a):
    from .efficiency import run_efficiency, write_report_from_run, EfficiencyMaps
    if a.ecmd == "report":
        p = write_report_from_run(a.run)
        print(p.read_text()); print("rewrote", p)
        return 0
    if a.ecmd == "apply":
        from .events import TruthTable
        from .channels import load_channel
        ch = load_channel(a.channel); t = TruthTable.load(a.sample)
        w = EfficiencyMaps.load(a.run).weights(ch, t)
        np.savez_compressed(a.out, weight=w)
        print(f"wrote {a.out}: {t.n} events, sum of ansatz weights {w.sum():.4f} (signal in the maps: {(w > 0).sum()})")
        return 0
    grids = None if a.grids in (None, "all") else [g.strip() for g in a.grids.split(",")]
    run_dir = run_efficiency(a.channel, load_site_config(), maps=tuple(a.maps.split(",")), grids=grids, out_root=a.out, slug=a.slug,
                             closure=not a.no_closure)
    print((run_dir / "report.md").read_text())
    print(f"run directory: {run_dir}")
    return 0


def _cmd_gibuu(a):
    from .grid import gibuu_campaign as gc
    cfg = load_site_config()
    if a.gcmd == "smoke":
        from .channels import load_channel
        from .theory.gibuu import GibuuSpec, read_job, run_local
        from .theory.models import ModelSpec
        from .theory.gibuu import stratum_specs
        ch = load_channel(a.channel); strata = stratum_specs(ModelSpec.load(a.model).params)
        g = strata[a.stratum] if a.stratum else next(iter(strata.values()))
        job = run_local(g, ch, cfg, num_ensembles=a.ensembles, seed=a.seed, job_name=a.job_name)
        t, info = read_job(job, g)
        sig = ch.is_signal(t)
        print(json.dumps({"job": str(job), **info, "n_signal": int(sig.sum()), "sigma_signal_1e-38cm2": float(t["weight"][sig].sum()),
                          "manifest": json.loads((job / "manifest_job.json").read_text())}, indent=2))
        return 0
    if a.gcmd == "plan":
        gc.plan(a.name, a.model, a.channel, cfg, stratum=a.stratum, n_jobs=a.n_jobs, num_ensembles=a.ensembles); return 0
    if a.gcmd == "submit-cmd":
        print(gc.submit_cmd(a.name, a.tar_label, memory=a.memory, disk=a.disk, lifetime=a.lifetime, n=a.n, stem=a.stem,
                            processes=[p.strip() for p in a.processes.split(",")] if a.processes else None)); return 0
    if a.gcmd == "record":
        gc.record_submission(a.name, a.jobid, a.cluster, a.tar_label, [p.strip() for p in a.processes.split(",")] if a.processes else None); return 0
    if a.gcmd == "status":
        gc.status(a.name); return 0
    if a.gcmd == "harvest":
        print(json.dumps(gc.harvest(a.name, workers=a.workers), indent=2)); return 0
    if a.gcmd == "merge":
        t = gc.merge(a.name, a.channel, cfg)
        print(json.dumps({"n_events": int(t.n), "sigma_flux_avg_per_nucleon_cm2": t.meta.get("sigma_flux_avg_per_nucleon_cm2"),
                          "flux_fraction_covered": t.meta.get("flux_fraction_covered"), "energy_points_missing": t.meta.get("energy_points_missing"),
                          "sigma_per_job_1e-38cm2": t.meta.get("sigma_per_job_1e-38cm2")}, indent=2))
        return 0
    return 2


def _cmd_selection(a):
    from .diagnostics import run_selection_comparison
    run_dir = run_selection_comparison(a.channel, load_site_config(), out_root=a.out, slug=a.slug)
    print((run_dir / "report.md").read_text())
    print(f"run directory: {run_dir}")
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
    p.add_argument("--efficiency-run", help="with --measurement all|a,b: an `ndp efficiency run` directory whose maps give the ansatz prediction")
    p.add_argument("--surrogate-kind", choices=("binned", "vbll"), default="binned",
                   help="binned: the grid's binned response; vbll: its efficiency + background with the migration from a ported VBLL model")
    _add_vbll_args(p)
    p.add_argument("--slug"); p.set_defaults(fn=_cmd_run)
    ps = sub.add_parser("surrogate", help="build / inspect detector surrogates").add_subparsers(dest="scmd", required=True)
    p = ps.add_parser("build"); p.add_argument("--channel", required=True)
    p.add_argument("--measurement", help="one name, a comma list, or `all` (every measurement YAML of the channel)")
    p.add_argument("--chunked", action="store_true", help="accumulate one playlist / file at a time (automatic with playlist products or several measurements)")
    p.add_argument("--source", choices=("mc", "artifacts"), default="mc")
    p.add_argument("--kind", choices=("binned", "parametric", "all", "vbll"), default="all",
                   help="vbll: wrap the grid's binned response with a ported VBLL model (surrogates/<channel>/<measurement>/vbll_<model>/) + closure on the MC truth")
    p.add_argument("--n-samples", type=int, default=20)
    _add_vbll_args(p); p.add_argument("--no-closure", action="store_true", help="(vbll) save the wrappers without folding the MC truth")
    p.add_argument("--out", help="output root (default surrogates/<channel>[/<measurement>])"); p.set_defaults(fn=_cmd_surrogate_build)
    p = ps.add_parser("inspect"); p.add_argument("path"); p.set_defaults(fn=_cmd_surrogate_inspect)
    p = ps.add_parser("train-vbll", help="train a VBLL surrogate on the channel's selected signal pairs (needs the ml environment)")
    p.add_argument("--channel", required=True); p.add_argument("--name", required=True, help="model name under surrogates/<channel>/_vbll/")
    p.add_argument("--inputs", default="E,px,py,pz,p,costheta"); p.add_argument("--outputs", default="E,px,py,pz,p,costheta")
    p.add_argument("--frame", default="beam", choices=("beam", "detector")); p.add_argument("--head", default="het", choices=("het", "standard"))
    p.add_argument("--d-embed", type=int, default=8); p.add_argument("--hidden", type=int, default=64); p.add_argument("--n-layers", type=int, default=3)
    p.add_argument("--noise-prior-scale", type=float, default=0.01); p.add_argument("--lr", type=float, default=1e-3)
    p.add_argument("--epochs", type=int, default=80); p.add_argument("--patience", type=int, default=10); p.add_argument("--batch-size", type=int, default=512)
    p.add_argument("--val-fraction", type=float, default=0.2); p.add_argument("--seed", type=int, default=42); p.add_argument("--threads", type=int)
    p.add_argument("--pairs", help="pairs .npz (default <data_dir>/cache/vbll_pairs_<channel>_<frame>.npz, collected if absent)")
    p.add_argument("--recollect", action="store_true"); p.set_defaults(fn=_cmd_surrogate_train_vbll)
    pd = sub.add_parser("data", help="data availability / caches").add_subparsers(dest="dcmd", required=True)
    p = pd.add_parser("status"); p.add_argument("--channel"); p.set_defaults(fn=_cmd_data_status)
    p = pd.add_parser("cache", help="build the truth/reco .npz caches from the channel's AnaTuples")
    p.add_argument("--channel", required=True); p.add_argument("--which", default="data,mc")
    p.add_argument("--reco-only", action="store_true", help="skip the (slow) Truth tree"); p.add_argument("--entry-stop", type=int)
    p.add_argument("--url", help="process this one file (local path or root:// URL) instead of the channel's files")
    p.add_argument("--kind", choices=("data", "mc")); p.add_argument("--out", help="output dir for --url (default <data_dir>/cache)")
    p.add_argument("--skim", action="store_true", help="also write derived-column skims (MC)")
    p.add_argument("--sidecar", action="store_true", help="also write manifest_<tag>.json with the channel's cutflows")
    p.set_defaults(fn=_cmd_data_cache)
    p = pd.add_parser("merge", help="merge harvested per-file products into playlist products")
    p.add_argument("--beam", default="FHC"); p.add_argument("--playlist", required=True); p.add_argument("--kind", default="data,mc")
    p.add_argument("--products-dir"); p.set_defaults(fn=_cmd_data_merge)
    pg = sub.add_parser("grid", help="grid campaign: plan worklists, track status, resubmit, harvest products").add_subparsers(dest="gcmd", required=True)
    p = pg.add_parser("harvest-pot", help="read every published file's Meta tree -> resources/minerva/opendata_pot_per_file.tsv"); p.add_argument("--workers", type=int, default=24); p.set_defaults(fn=_cmd_grid)
    p = pg.add_parser("plan"); p.add_argument("name"); p.add_argument("--channel", required=True); p.add_argument("--beams", default="FHC")
    p.add_argument("--kinds", default="data,mc"); p.add_argument("--playlists", help="comma list, e.g. 1A,1B (default: all)")
    p.add_argument("--files-per-process", help="e.g. mc=4,data=60"); p.add_argument("--pnfs-base"); p.set_defaults(fn=_cmd_grid)
    p = pg.add_parser("status"); p.add_argument("name"); p.add_argument("--no-pot-check", action="store_true"); p.set_defaults(fn=_cmd_grid)
    p = pg.add_parser("resubmit"); p.add_argument("name"); p.set_defaults(fn=_cmd_grid)
    p = pg.add_parser("harvest"); p.add_argument("name"); p.add_argument("--products-dir"); p.add_argument("--playlists")
    p.add_argument("--workers", type=int, default=4); p.add_argument("--no-archive", action="store_true", help="skip the full truth tables (keep them on PNFS)")
    p.set_defaults(fn=_cmd_grid)
    p = pg.add_parser("stage-worklists", help="upload worklists (or resubmit lists) to PNFS scratch as job inputs"); p.add_argument("name")
    p.add_argument("--files", nargs="*", help="specific files (default: every worklist of the campaign)"); p.set_defaults(fn=_cmd_grid)
    p = pg.add_parser("submit-cmd", help="print the jobsub-lite submit command for one worklist"); p.add_argument("name"); p.add_argument("worklist")
    p.add_argument("--tar-label", default="ndp-stream-v1"); p.add_argument("--max-processes", type=int)
    p.add_argument("--file", help="a staged resubmit list to submit instead of the worklist"); p.add_argument("--stem")
    p.add_argument("--disk", help="override the disk request, e.g. 8GB"); p.set_defaults(fn=_cmd_grid)
    p = sub.add_parser("flux", help="channel flux table summary / export"); p.add_argument("--channel", required=True); p.add_argument("--out"); p.set_defaults(fn=_cmd_flux)
    p = sub.add_parser("signal", help="apply a channel's truth-level signal definition to the cached MC; writes a diagnostics run")
    p.add_argument("--channel", required=True); p.add_argument("--cache", help="truth .npz (default: the channel's MC cache)")
    p.add_argument("--out"); p.add_argument("--slug"); p.set_defaults(fn=_cmd_signal)
    pg2 = sub.add_parser("gibuu", help="GiBUU generation: local smoke job, grid campaign plan/submit/status/harvest/merge").add_subparsers(dest="gcmd", required=True)
    p = pg2.add_parser("smoke", help="one local GiBUU job through the runner"); p.add_argument("model"); p.add_argument("--channel", required=True)
    p.add_argument("--ensembles", type=int, default=100); p.add_argument("--seed", type=int); p.add_argument("--job-name", default="smoke")
    p.add_argument("--stratum"); p.set_defaults(fn=_cmd_gibuu)
    p = pg2.add_parser("plan", help="campaign record + card template + flux file"); p.add_argument("name"); p.add_argument("model"); p.add_argument("--channel", required=True)
    p.add_argument("--stratum", help="stratum of the model spec (models with a `strata:` block)")
    p.add_argument("--n-jobs", type=int); p.add_argument("--ensembles", type=int); p.set_defaults(fn=_cmd_gibuu)
    p = pg2.add_parser("submit-cmd", help="print the jobsub-lite submit command"); p.add_argument("name"); p.add_argument("--tar-label", required=True)
    p.add_argument("--memory", default="2500MB"); p.add_argument("--disk", default="2GB"); p.add_argument("--lifetime", default="3h")
    p.add_argument("--n", type=int); p.add_argument("--stem"); p.add_argument("--processes", help="comma list of 4-digit process ids to rerun"); p.set_defaults(fn=_cmd_gibuu)
    p = pg2.add_parser("record", help="record a submission in campaign.json"); p.add_argument("name"); p.add_argument("--jobid", required=True)
    p.add_argument("--cluster"); p.add_argument("--tar-label", required=True); p.add_argument("--processes"); p.set_defaults(fn=_cmd_gibuu)
    p = pg2.add_parser("status"); p.add_argument("name"); p.set_defaults(fn=_cmd_gibuu)
    p = pg2.add_parser("harvest"); p.add_argument("name"); p.add_argument("--workers", type=int, default=4); p.set_defaults(fn=_cmd_gibuu)
    p = pg2.add_parser("merge"); p.add_argument("name"); p.add_argument("--channel", required=True); p.set_defaults(fn=_cmd_gibuu)
    pe = sub.add_parser("efficiency", help="selection-efficiency maps, background tables and the factorised-ansatz closure from the binned surrogates")
    pes = pe.add_subparsers(dest="ecmd")
    p = pes.add_parser("run", help="write runs/<date>_efficiency_<channel>/"); p.add_argument("--channel", required=True)
    p.add_argument("--maps", default="muon_p_costheta,proton_p_costheta,muon_costheta,proton_costheta")
    p.add_argument("--grids", default="all", help="released grids to export (comma list or `all`)")
    p.add_argument("--no-closure", action="store_true", help="skip the ansatz closure pass over the MC")
    p.add_argument("--out"); p.add_argument("--slug"); p.set_defaults(fn=_cmd_efficiency)
    p = pes.add_parser("report", help="re-render report.md of a finished efficiency run from its JSON tables")
    p.add_argument("--run", required=True); p.set_defaults(fn=_cmd_efficiency)
    p = pes.add_parser("apply", help="weight a truth sample (TruthTable npz) with the maps of an efficiency run")
    p.add_argument("--channel", required=True); p.add_argument("--run", required=True); p.add_argument("--sample", required=True)
    p.add_argument("--out", required=True); p.set_defaults(fn=_cmd_efficiency)
    p = sub.add_parser("selection", help="apply a channel's reco selection to the cached data + MC; cutflow, purity/efficiency, data-vs-MC figures")
    p.add_argument("--channel", required=True); p.add_argument("--out"); p.add_argument("--slug"); p.set_defaults(fn=_cmd_selection)
    a = ap.parse_args(argv)
    return a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
