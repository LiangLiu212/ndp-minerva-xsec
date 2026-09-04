"""Run GENIE (gevgen + gntpc) for a channel's flux and target, absolutely normalised.

The runner talks to a GENIE installation through a *genie-agent environment snapshot*
(a JSON of environment variables captured from the install's setup script, see
../genie-dev/genie-agent/scripts/refresh_genie_env.py), so the platform's own Python
environment never leaks into the GENIE child process.

Normalisation: events are unweighted, so d(sigma)/dx per nucleon =
sigma_avg_per_nucleon * N_bin / (N_gen * dx), with sigma_avg the flux-averaged total
cross section of the target mix from the spline file (`splines.py`). The number of
generated events must be counted from the gst file, not assumed from `-n`.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path

import numpy as np

from ..events import TruthTable, Normalization
from ..io import sha256_text, dump_json
from ..adapters.genie_gst import read_gst
from . import flux as fluxmod
from . import splines as splmod


@dataclass
class GenieSpec:
    tune: str = "G18_02a_00_000"
    generator_list: str = "CC"
    n_events: int = 20000
    nu_pdg: int = 14
    target_mix: dict = field(default_factory=lambda: {1000060120: 0.9231, 1000010010: 0.0769})  # mass fractions
    e_min: float = 0.0
    e_max: float = 100.0
    seed: int = 1
    run_number: int = 1
    n_jobs: int = 1                 # parallel gevgen processes (seeds seed..seed+n_jobs-1)
    env_json: str | None = None     # genie-agent env snapshot; default from site config
    splines: str | None = None      # spline XML; default from site config
    gxmlpath: str | None = None     # extra tune directory for custom tunes
    extra_args: list = field(default_factory=list)

    def to_dict(self) -> dict:
        d = asdict(self)
        d["target_mix"] = {str(k): v for k, v in self.target_mix.items()}
        return d

    def fingerprint(self, flux_source: str) -> str:
        return sha256_text(json.dumps(self.to_dict(), sort_keys=True) + "|" + flux_source)[:16]


def load_genie_env(env_json: str | Path) -> dict:
    env = json.load(open(env_json))
    env = env.get("env", env)
    env = {k: str(v) for k, v in env.items()}
    env.setdefault("HOME", os.environ.get("HOME", "/tmp"))
    return env


def target_arg(mix: dict) -> str:
    items = sorted(((int(k), float(v)) for k, v in mix.items()), key=lambda kv: -kv[1])
    tot = sum(v for _, v in items)
    return ",".join(f"{k}[{v / tot:.6f}]" for k, v in items)


def _run(cmd: list[str], env: dict, cwd: Path, log_stem: Path, timeout: int) -> dict:
    t0 = time.time()
    r = subprocess.run(cmd, env=env, cwd=str(cwd), capture_output=True, text=True, timeout=timeout)
    (log_stem.with_suffix(".stdout")).write_text(r.stdout)
    (log_stem.with_suffix(".stderr")).write_text(r.stderr)
    return {"cmd": cmd, "returncode": r.returncode, "seconds": round(time.time() - t0, 1)}


def generate(spec: GenieSpec, flux_edges: np.ndarray, flux_density: np.ndarray, flux_source: str,
             workdir: str | Path, site_cfg=None, timeout: int = 6 * 3600, cache: bool = True) -> TruthTable:
    """Produce an absolutely normalised TruthTable for `spec` under the given flux.

    Output layout: <workdir>/<fingerprint>/{flux.root, job_k.ghep.root, job_k.gst.root,
    genie_run.json, truth.npz}. A finished directory is reused when `cache` is set.
    """
    env_json = spec.env_json or (str(site_cfg.genie_env_json) if site_cfg and site_cfg.genie_env_json else None)
    splines = spec.splines or (str(site_cfg.genie_splines) if site_cfg and site_cfg.genie_splines else None)
    if not env_json or not Path(env_json).exists():
        raise FileNotFoundError("no GENIE environment snapshot (spec.env_json / site genie_env_json)")
    if not splines or not Path(splines).exists():
        raise FileNotFoundError("no GENIE spline XML (spec.splines / site genie_splines)")
    workdir = Path(workdir)
    fp = spec.fingerprint(flux_source)
    out = workdir / f"genie_{spec.tune}_{fp}"
    truth_npz = out / "truth.npz"
    if cache and truth_npz.exists():
        t = TruthTable.load(truth_npz)
        t.meta["cache_hit"] = True
        return t
    out.mkdir(parents=True, exist_ok=True)
    env = load_genie_env(env_json)
    if spec.gxmlpath:
        env["GXMLPATH"] = spec.gxmlpath + (":" + env["GXMLPATH"] if env.get("GXMLPATH") else "")

    # Logging verbosity is a runtime detail, not physics: keep it out of the fingerprint.
    # GENIE ships Messenger_laconic.xml; without it gevgen writes ~800 MB of stdout per 30k events.
    runtime_args = ["--message-thresholds", "Messenger_laconic.xml"] if not spec.extra_args or \
        "--message-thresholds" not in spec.extra_args else []
    # flux histogram for gevgen (density shape only matters; units irrelevant)
    flux_root = fluxmod.write_th1_root(out / "flux.root", flux_edges, flux_density, name="flux")
    per_job = int(np.ceil(spec.n_events / max(spec.n_jobs, 1)))
    procs, logs = [], []
    for k in range(spec.n_jobs):
        stem = out / f"job_{k}"
        cmd = ["gevgen", "-n", str(per_job), "-p", str(spec.nu_pdg), "-t", target_arg(spec.target_mix),
               "-e", f"{spec.e_min},{spec.e_max}", "-f", f"{flux_root},flux",
               "--cross-sections", splines, "--tune", spec.tune, "--event-generator-list", spec.generator_list,
               "--seed", str(spec.seed + k), "-r", str(spec.run_number + k), "-o", f"{stem}.ghep.root",
               *spec.extra_args, *runtime_args]
        fo, fe = open(f"{stem}.gevgen.stdout", "w"), open(f"{stem}.gevgen.stderr", "w")
        procs.append((k, subprocess.Popen(cmd, env=env, cwd=str(out), stdout=fo, stderr=fe), fo, fe, time.time(), cmd))
    for k, p, fo, fe, t0, cmd in procs:
        rc = p.wait(timeout=timeout)
        fo.close(); fe.close()
        logs.append({"job": k, "step": "gevgen", "cmd": cmd, "returncode": rc, "seconds": round(time.time() - t0, 1)})
        if rc != 0:
            raise RuntimeError(f"gevgen job {k} failed (rc={rc}); see {out}/job_{k}.gevgen.stderr")
    tables = []
    for k in range(spec.n_jobs):
        stem = out / f"job_{k}"
        r = _run(["gntpc", "-i", f"{stem}.ghep.root", "-f", "gst", "-o", f"{stem}.gst.root", "--tune", spec.tune],
                 env, out, out / f"job_{k}.gntpc", timeout=3600)
        logs.append({"job": k, "step": "gntpc", **r})
        if r["returncode"] != 0:
            raise RuntimeError(f"gntpc job {k} failed; see {out}/job_{k}.gntpc.stderr")
        tables.append(read_gst(f"{stem}.gst.root"))
    t = TruthTable.concatenate(tables) if len(tables) > 1 else tables[0]

    # absolute normalisation from the splines + flux
    spl = splmod.total_xsec_splines(splines, spec.nu_pdg, list(spec.target_mix), "Weak[CC]" if
                                    spec.generator_list.upper().startswith("CC") else "Weak",
                                    cache_dir=workdir / "spline_cache")
    E, sig_nuc = splmod.per_nucleon_total(spl, spec.target_mix)
    sigma_avg = splmod.flux_averaged_per_nucleon(E, sig_nuc, flux_edges, flux_density, spec.e_min, spec.e_max)
    n_gen = t.n
    t.meta.update({
        "generator": f"GENIE {spec.tune} ({spec.generator_list})", "genie_spec": spec.to_dict(),
        "genie_env_json": env_json, "splines": splines, "flux_source": flux_source,
        "sigma_flux_avg_per_nucleon_cm2": sigma_avg, "n_generated": n_gen,
        "n_splines_used": {str(k): v["n_splines"] for k, v in spl.items()},
        "norm": Normalization(kind="xsec_per_nucleon", xsec_per_unit_weight=sigma_avg / float(t["weight"].sum()),
                              notes="sigma_avg from spline file x flux table; events unweighted").to_dict(),
        "genie_logs": logs, "source": str(out),
    })
    t.save(truth_npz)
    dump_json({"spec": spec.to_dict(), "flux_source": flux_source, "sigma_flux_avg_per_nucleon_cm2": sigma_avg,
               "n_generated": n_gen, "logs": logs, "fingerprint": fp}, out / "genie_run.json")
    return t
