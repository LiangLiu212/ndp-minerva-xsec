"""A GiBUU generation campaign on the FNAL grid: N independent jobs of one model spec, one seed each.

    grid/campaigns/<name>/campaign.json     spec, fingerprint, n_jobs, seed base, pnfs_out, jobs, processes
    <pnfs_out>/<%04d>/                      what grid/gibuu_worker.sh leaves per process (FinalEvents.dat.gz, the
                                            absorption cross-section file, the flux check, card, log, manifest_<%04d>.json)
    runs/_generator_cache/gibuu_<fp>/jobs/<%04d>/   harvested copies; `merge` -> truth.npz next to them

Reuses the PNFS helpers of `campaign.py` (token discovery, xrootd namespace, listing, copying).
"""
from __future__ import annotations

import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from ..io import dump_json, timestamp
from .campaign import campaign_dir, pnfs_copy, pnfs_ls, save_campaign


def load(name: str) -> dict:
    p = campaign_dir(name) / "campaign.json"
    if not p.exists():
        raise FileNotFoundError(f"no campaign {name!r} ({p}); run `ndp gibuu plan {name} <model.yaml>` first")
    c = json.loads(p.read_text())
    if c.get("kind") != "gibuu":
        raise ValueError(f"campaign {name!r} is not a GiBUU campaign")
    return c


def plan(name: str, model_path: str, channel_name: str, cfg, *, stratum: str | None = None, n_jobs: int | None = None,
         num_ensembles: int | None = None, pnfs_base: str | None = None, log=print) -> dict:
    """Prepare the cache directory (flux file, card template with the worker's placeholders) and the campaign record
    for one stratum of the model spec (`--stratum`; a spec without strata has the single stratum `all`)."""
    from ..channels import load_channel
    from ..theory.gibuu import gibuu_paths, prepare, stratum_specs
    from ..theory.models import ModelSpec
    spec_m = ModelSpec.load(model_path)
    if spec_m.kind != "gibuu":
        raise ValueError(f"{model_path}: kind {spec_m.kind!r} is not gibuu")
    strata = stratum_specs(spec_m.params)
    if stratum is None and len(strata) == 1:
        stratum = next(iter(strata))
    if stratum not in strata:
        raise ValueError(f"model {spec_m.name} has strata {list(strata)}; pass --stratum")
    g = strata[stratum]
    if n_jobs:
        g.n_jobs = int(n_jobs)
    if num_ensembles:
        g.num_ensembles = int(num_ensembles)
    ch = load_channel(channel_name)
    prep = prepare(g, ch, cfg)
    # the template for the workers: everything filled except the three run-time placeholders
    tmpl = _card_template(g, prep["dir"], cfg)
    import getpass
    user = getpass.getuser()
    base = pnfs_base or f"/pnfs/dune/scratch/users/{user}/ndp-gibuu"
    c = {"name": name, "kind": "gibuu", "created": timestamp(), "model": str(Path(model_path).resolve()), "model_name": spec_m.name,
         "stratum": stratum, "channel": ch.name, "spec": g.to_dict(), "fingerprint": prep["fingerprint"], "cache_dir": str(prep["dir"]),
         "card_template": str(tmpl), "flux_file": prep["flux_file"], "n_jobs": g.n_jobs, "seed_base": g.seed,
         "expected_events_per_job": g.events_per_job(), "pnfs_out": f"{base}/{name}", "jobs": [],
         "processes": {f"{k:04d}": {"seed": g.seed + k, "status": "planned", "attempts": 0} for k in range(g.n_jobs)},
         "gibuu_paths": {k: str(v) for k, v in gibuu_paths(cfg.repo_root).items()}}
    save_campaign(c)
    log(f"campaign {name} (stratum {stratum}): {g.n_jobs} jobs x {g.num_ensembles} ensembles x {g.num_runs} run(s) "
        f"(~{g.events_per_job()} test nucleons each), equal-weights mode {g.equal_weights_mode} ceiling {g.equal_weights_max}, fingerprint {prep['fingerprint']}")
    log(f"card template {tmpl}; flux {prep['flux_file']}; outputs -> {c['pnfs_out']}")
    return c


def _card_template(g, cache_dir: Path, cfg) -> Path:
    """The card with every spec value filled and the worker's three placeholders kept."""
    from ..theory.gibuu import substitutions
    tmpl = Path(g.template)
    if not tmpl.is_absolute():
        tmpl = cfg.repo_root / tmpl
    text = tmpl.read_text()
    subs = substitutions(g, flux_file="__FLUX_FILE__", path_to_input="__PATH_TO_INPUT__", seed=0)
    subs["__SEED__"] = "__SEED__"
    for k, v in subs.items():
        text = text.replace(k, v)
    out = cache_dir / "card.job.tmpl"
    out.write_text(text)
    return out


def submit_cmd(name: str, tar_label: str, *, memory: str = "2500MB", disk: str = "2GB", lifetime: str = "3h", n: int | None = None,
               stem: str | None = None, processes: list | None = None) -> str:
    """The jobsub-lite skill command for the campaign (or for a resubmission of the listed processes)."""
    c = load(name)
    n_proc = n or (len(processes) if processes else c["n_jobs"])
    seed_base = c["seed_base"]
    if processes:   # resubmission: the k-th new process reruns processes[k] with a fresh seed
        seed_base = c["seed_base"] + 1000 * (1 + max(c["processes"][p]["attempts"] for p in processes))
    return (f"python3 .claude/skills/jobsub-lite/scripts/jobsub.py submit --worker grid/gibuu_worker.sh -N {n_proc} --tar-label {tar_label} "
            f"--memory {memory} --disk {disk} --expected-lifetime {lifetime} --jobsub-arg=--onsite --pnfs-out {c['pnfs_out']} "
            f"--runtype ndpgibuu --stem {stem or name} -- -R @TAR_DIR@ -O @PNFS_OUT@ -C card -S {seed_base}")


def record_submission(name: str, jobid: str, cluster: str | None, tar_label: str, processes: list | None = None) -> None:
    c = load(name)
    c["jobs"].append({"jobid": jobid, "cluster": cluster, "tar_label": tar_label, "submitted": timestamp(), "processes": processes})
    for k, p in c["processes"].items():
        if processes is None or k in processes:
            p["status"] = "submitted"; p["attempts"] += 1
    save_campaign(c)


def status(name: str, log=print) -> dict:
    """Refresh every process from the sidecars on PNFS: done / failed / missing."""
    c = load(name)
    side = campaign_dir(name) / "sidecars"; side.mkdir(exist_ok=True)
    try:
        dirs = pnfs_ls(c["pnfs_out"])
    except RuntimeError as e:
        log(f"cannot list {c['pnfs_out']}: {e}"); dirs = []
    present = {Path(d).name for d in dirs}
    for k, p in c["processes"].items():
        if p["status"] == "done":
            continue
        if k not in present:
            continue
        local = side / f"manifest_{k}.json"
        if not local.exists():
            try:
                pnfs_copy(f"{c['pnfs_out']}/{k}/manifest_{k}.json", local)
            except Exception as e:  # noqa: BLE001
                log(f"{k}: sidecar not readable yet ({str(e)[:80]})"); continue
        try:
            s = json.loads(local.read_text())
        except json.JSONDecodeError:
            local.unlink(missing_ok=True); continue
        p.update({"status": "done" if s.get("status") == "ok" and s.get("n_events", 0) > 0 else "failed",
                  "n_events": s.get("n_events"), "wall_s": s.get("wall_s"), "xsec": s.get("xsec_file_1e-38cm2"), "cluster": s.get("cluster"),
                  "pnfs_dir": f"{c['pnfs_out']}/{k}"})
    save_campaign(c)
    counts = {}
    for p in c["processes"].values():
        counts[p["status"]] = counts.get(p["status"], 0) + 1
    log(f"campaign {name}: " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())))
    return counts


def harvest(name: str, workers: int = 4, log=print) -> dict:
    """Copy every done process into runs/_generator_cache/gibuu_<fp>/jobs/<%04d>/."""
    c = load(name)
    root = Path(c["cache_dir"]) / "jobs"
    names = ["FinalEvents.dat.gz", "neutrino_absorption_cross_section_ALL.dat", "neutrino_initialized_energyFlux.dat", "job.card", "gibuu.log.gz"]
    todo = [(k, p) for k, p in c["processes"].items() if p["status"] == "done" and not (root / k / "harvested.json").exists()]

    def one(item):
        k, p = item
        dst = root / k; dst.mkdir(parents=True, exist_ok=True)
        for n in names + [f"manifest_{k}.json"]:
            try:
                pnfs_copy(f"{p['pnfs_dir']}/{n}", dst / n)
            except Exception as e:  # noqa: BLE001
                if n in ("FinalEvents.dat.gz", f"manifest_{k}.json", "neutrino_absorption_cross_section_ALL.dat"):
                    raise RuntimeError(f"{k}: {n}: {e}") from None
        if (dst / "FinalEvents.dat.gz").stat().st_size == 0:
            raise RuntimeError(f"{k}: empty FinalEvents.dat.gz")
        (dst / "harvested.json").write_text(json.dumps({"harvested": timestamp(), "from": p["pnfs_dir"]}))
        return k

    n_ok, errors = 0, []
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for fu in as_completed([ex.submit(one, it) for it in todo]):
            try:
                k = fu.result(); n_ok += 1; c["processes"][k]["harvested"] = str(root / k); log(f"harvested {k}")
            except Exception as e:  # noqa: BLE001
                errors.append(str(e)[:200]); log(f"harvest failed: {str(e)[:120]}")
    save_campaign(c)
    return {"harvested": n_ok, "already": len([1 for k, p in c["processes"].items() if p["status"] == "done"]) - len(todo), "failed": len(errors), "errors": errors[:10]}


def merge(name: str, channel_name: str, cfg, log=print):
    """Merge every harvested job into the cache's truth.npz (weights / K, cross-section check per job)."""
    from ..channels import load_channel
    from ..theory.gibuu import GibuuSpec, merge_jobs
    c = load(name)
    g = GibuuSpec.from_params({"target": {"Z": c["spec"]["target_Z"], "A": c["spec"]["target_A"]}, **_spec_params(c["spec"])})
    root = Path(c["cache_dir"]) / "jobs"
    jobs = sorted(d for d in root.iterdir() if d.is_dir() and (d / "harvested.json").exists())
    if not jobs:
        raise FileNotFoundError(f"no harvested jobs under {root}")
    t = merge_jobs(jobs, g, load_channel(channel_name), cfg, log=log)
    c["merged"] = {"n_jobs": len(jobs), "n_events": int(t.n), "sigma_flux_avg_per_nucleon_cm2": t.meta["sigma_flux_avg_per_nucleon_cm2"], "when": timestamp()}
    save_campaign(c)
    return t


def _spec_params(d: dict) -> dict:
    """GibuuSpec.to_dict() -> the model-spec params layout GibuuSpec.from_params expects."""
    return {"process_ID": d["process_ID"], "flavor_ID": d["flavor_ID"], "include": d["include"],
            "card": {"num_ensembles": d["num_ensembles"], "num_time_steps": d["num_time_steps"], "delta_T": d["delta_T"], "num_runs": d["num_runs"],
                     "medium_switch": d["medium_switch"], "medium_switch_delta": d["medium_switch_delta"], "medium_switch_coll": d["medium_switch_coll"],
                     "eqs_type": d["eqs_type"], "delta_pot": d["delta_pot"], "density_switch_static": d["density_switch_static"], "apply_cuts": d["apply_cuts"],
                     "equal_weights_mode": d.get("equal_weights_mode", 0), "equal_weights_max": d.get("equal_weights_max", -1.0)},
            "flux": {"bin_width_gev": d["flux_bin_width_gev"], "e_max_gev": d["flux_e_max_gev"]}, "n_jobs": d["n_jobs"], "seed": d["seed"],
            "template": d["template"]}
