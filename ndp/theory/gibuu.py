"""GiBUU runner: a job card from a template + the channel flux, local runs, and the merge of grid jobs
into one absolutely normalised TruthTable under runs/_generator_cache/gibuu_<fingerprint>/.

The physics of a run lives in a model spec (`models/*.yaml`, `kind: gibuu`): target, process, the
included channels, the in-medium switches and the ensemble/time-step settings; `write_card` fills
`resources/gibuu/minerva_me_numu_CC.job.tmpl` with them. The flux is the channel's flux table rebinned
to equidistant bins (GiBUU's `nuExp = 99` contract), so the flux-averaged cross section GiBUU reports
refers to the same energy range as the channel's `phi_per_pot_cm2`.

Normalisation (verified on the smoke sample, docs/decisions.md 2026-09-04): every FinalEvents row
carries the event's perweight in 1e-38 cm^2 per nucleon, and the perweights of one run sum to the
flux-averaged total cross section; the adapter divides by the number of runs. K independent jobs
(different seeds) are merged by prescaling every job's weights by 1/K, so the merged weights still sum
to the mean cross section per nucleon.
"""
from __future__ import annotations

import gzip
import json
import os
import shutil
import subprocess
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

from ..adapters.gibuu_finalevents import read_finalevents
from ..events import Normalization, TruthTable
from ..io import dump_json, sha256_file, sha256_text, timestamp
from . import flux as fluxmod

DEFAULT_TEMPLATE = "resources/gibuu/minerva_me_numu_CC.job.tmpl"
XSEC_FILE = "neutrino_absorption_cross_section_ALL.dat"
FLUX_CHECK_FILE = "neutrino_initialized_energyFlux.dat"
GIBUU_VERSION = "GiBUU Release 2025 patch 5"


@dataclass
class GibuuSpec:
    target_Z: int = 6
    target_A: int = 12
    process_ID: int = 2                      # 2 CC, 3 NC, -2 anti-CC, -3 anti-NC
    flavor_ID: int = 2                       # 1 e, 2 mu, 3 tau
    include: dict = field(default_factory=lambda: {"QE": True, "DELTA": True, "RES": True, "1pi": True,
                                                   "DIS": True, "2p2hQE": True, "2pi": True})
    num_ensembles: int = 1000
    num_time_steps: int = 150
    delta_T: float = 0.2
    num_runs: int = 1                        # num_runs_SameEnergy
    medium_switch: bool = True               # &width_Baryon mediumSwitch
    medium_switch_delta: bool = True
    medium_switch_coll: bool = False
    eqs_type: int = 5
    delta_pot: int = 1
    density_switch_static: int = 2
    apply_cuts: int = 2                      # 2: only unbound particles reach FinalEvents.dat
    equal_weights_mode: int = 0              # 0 weighted events; 2 MC rejection against equal_weights_max -> unit weights
    equal_weights_max: float = -1.0          # mode 2 ceiling (per-nucleon sigma units, 1e-38 cm^2): GiBUU aborts if an event exceeds it
    flux_bin_width_gev: float = 0.5
    flux_e_max_gev: float = 100.0
    n_jobs: int = 100
    seed: int = 20260915                     # base seed; job k uses seed + k (GiBUU `Seed` is a Fortran integer)
    template: str = DEFAULT_TEMPLATE

    KNOWN_INCLUDE = ("QE", "DELTA", "RES", "1pi", "DIS", "2p2hQE", "2pi")

    @staticmethod
    def from_params(p: dict) -> "GibuuSpec":
        """Model-spec params -> GibuuSpec (unknown keys such as `notes`, `local`, `description` are ignored)."""
        tgt = p.get("target") or {}
        card = p.get("card") or {}
        fl = p.get("flux") or {}
        inc = dict(GibuuSpec().include)
        inc.update({str(k): bool(v) for k, v in (p.get("include") or {}).items()})
        unknown = sorted(set(inc) - set(GibuuSpec.KNOWN_INCLUDE))
        if unknown:
            raise ValueError(f"unknown GiBUU include flags {unknown} (known: {GibuuSpec.KNOWN_INCLUDE})")
        kw = dict(target_Z=int(tgt.get("Z", 6)), target_A=int(tgt.get("A", 12)),
                  process_ID=int(p.get("process_ID", 2)), flavor_ID=int(p.get("flavor_ID", 2)), include=inc,
                  num_ensembles=int(card.get("num_ensembles", 1000)), num_time_steps=int(card.get("num_time_steps", 150)),
                  delta_T=float(card.get("delta_T", 0.2)), num_runs=int(card.get("num_runs", 1)),
                  medium_switch=bool(card.get("medium_switch", True)), medium_switch_delta=bool(card.get("medium_switch_delta", True)),
                  medium_switch_coll=bool(card.get("medium_switch_coll", False)), eqs_type=int(card.get("eqs_type", 5)),
                  delta_pot=int(card.get("delta_pot", 1)), density_switch_static=int(card.get("density_switch_static", 2)),
                  apply_cuts=int(card.get("apply_cuts", 2)), equal_weights_mode=int(card.get("equal_weights_mode", 0)),
                  equal_weights_max=float(card.get("equal_weights_max", -1.0)), flux_bin_width_gev=float(fl.get("bin_width_gev", 0.5)),
                  flux_e_max_gev=float(fl.get("e_max_gev", 100.0)), n_jobs=int(p.get("n_jobs", 100)),
                  seed=int(p.get("seed", 20260915)), template=str(p.get("template", DEFAULT_TEMPLATE)))
        if not (-2 ** 31 < kw["seed"] + kw["n_jobs"] < 2 ** 31):
            raise ValueError("GiBUU seeds must fit a 32-bit Fortran integer")
        if kw["equal_weights_mode"] == 2 and kw["equal_weights_max"] <= 0:
            raise ValueError("equal_weights_mode 2 needs a positive equal_weights_max (from a mode-1 pilot run)")
        return GibuuSpec(**kw)

    def to_dict(self) -> dict:
        return asdict(self)

    def fingerprint(self, flux_source: str, template_sha: str) -> str:
        return sha256_text(json.dumps(self.to_dict(), sort_keys=True) + "|" + flux_source + "|" + template_sha)[:16]

    def events_per_job(self) -> int:
        """Rule of thumb from GiBUU's card comments: events ~ A x numEnsembles x num_runs_SameEnergy."""
        return self.target_A * self.num_ensembles * self.num_runs


def _fb(b: bool) -> str:
    return "T" if b else "F"


def substitutions(spec: GibuuSpec, *, flux_file: str, path_to_input: str, seed: int, num_ensembles: int | None = None) -> dict:
    inc = spec.include
    return {"__PROCESS_ID__": str(spec.process_ID), "__FLAVOR_ID__": str(spec.flavor_ID), "__FLUX_FILE__": str(flux_file),
            "__INCLUDE_QE__": _fb(inc["QE"]), "__INCLUDE_DELTA__": _fb(inc["DELTA"]), "__INCLUDE_RES__": _fb(inc["RES"]),
            "__INCLUDE_1PI__": _fb(inc["1pi"]), "__INCLUDE_DIS__": _fb(inc["DIS"]), "__INCLUDE_2P2HQE__": _fb(inc["2p2hQE"]),
            "__INCLUDE_2PI__": _fb(inc["2pi"]), "__TARGET_Z__": str(spec.target_Z), "__TARGET_A__": str(spec.target_A),
            "__DENSITY_SWITCH_STATIC__": str(spec.density_switch_static),
            "__NUM_ENSEMBLES__": str(int(num_ensembles or spec.num_ensembles)), "__NUM_TIME_STEPS__": str(spec.num_time_steps),
            "__DELTA_T__": repr(float(spec.delta_T)), "__NUM_RUNS__": str(spec.num_runs), "__PATH_TO_INPUT__": str(path_to_input),
            "__MEDIUM_SWITCH__": ".true." if spec.medium_switch else ".false.",
            "__MEDIUM_SWITCH_DELTA__": ".true." if spec.medium_switch_delta else ".false.",
            "__MEDIUM_SWITCH_COLL__": ".true." if spec.medium_switch_coll else ".false.",
            "__EQS_TYPE__": str(spec.eqs_type), "__DELTA_POT__": str(spec.delta_pot), "__APPLY_CUTS__": str(spec.apply_cuts),
            "__EQUAL_WEIGHTS_MODE__": str(spec.equal_weights_mode), "__EQUAL_WEIGHTS_MAX__": repr(float(spec.equal_weights_max)),
            "__SEED__": str(int(seed))}


def write_card(spec: GibuuSpec, *, flux_file: str, path_to_input: str, seed: int, out: str | Path,
               repo_root: str | Path | None = None, num_ensembles: int | None = None) -> Path:
    """Fill the template with the spec; refuses to leave a placeholder behind."""
    tmpl = Path(spec.template)
    if not tmpl.is_absolute():
        tmpl = Path(repo_root or Path(__file__).resolve().parents[2]) / tmpl
    text = tmpl.read_text()
    for k, v in substitutions(spec, flux_file=flux_file, path_to_input=path_to_input, seed=seed, num_ensembles=num_ensembles).items():
        text = text.replace(k, v)
    left = sorted(set(w for w in text.split() if w.startswith("__") and w.endswith("__")))
    if left:
        raise ValueError(f"unfilled placeholders in the GiBUU card: {left}")
    out = Path(out); out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text)
    return out


def template_sha(spec: GibuuSpec, repo_root: Path) -> str:
    tmpl = Path(spec.template)
    return sha256_file(tmpl if tmpl.is_absolute() else repo_root / tmpl)


def gibuu_paths(repo_root: Path) -> dict:
    """Binary and input tables: $GIBUU / $GIBUU_INPUT if set (activate.sh), else the in-repo build."""
    ext = repo_root / "external" / "gibuu"
    gx = Path(os.environ["GIBUU"]) / "objects" / "GiBUU.x" if os.environ.get("GIBUU") else ext / "release2025" / "objects" / "GiBUU.x"
    inp = Path(os.environ["GIBUU_INPUT"]) if os.environ.get("GIBUU_INPUT") else ext / "buuinput2025"
    return {"gibuu_x": gx, "path_to_input": inp.resolve() if inp.exists() else inp}


def prepare(spec: GibuuSpec, channel, cfg) -> dict:
    """The cache directory for this spec + channel flux: writes the GiBUU flux file and records the fingerprint."""
    fl = fluxmod.load_channel_flux(channel, cfg.repo_root)
    flux_source = str(channel.normalization["flux_table"])
    tsha = template_sha(spec, cfg.repo_root)
    fp = spec.fingerprint(flux_source, tsha)
    out = cfg.runs / "_generator_cache" / f"gibuu_{fp}"
    out.mkdir(parents=True, exist_ok=True)
    centres, values = fluxmod.rebin_uniform(fl["edges"], fl["density_cm2_pot_gev"], spec.flux_bin_width_gev, spec.flux_e_max_gev)
    flux_file = fluxmod.write_gibuu_flux(out / "flux_gibuu.dat", centres, values,
                                         source=f"{flux_source} rebinned to {spec.flux_bin_width_gev} GeV bins, 0-{spec.flux_e_max_gev:g} GeV")
    info = {"fingerprint": fp, "spec": spec.to_dict(), "flux_source": flux_source, "template_sha256": tsha,
            "flux_file": str(flux_file), "flux_file_sha256": sha256_file(flux_file),
            "phi_0_100_cm2_per_pot": fl["phi_0_100_cm2_per_pot"], "flux_mean_energy_gev": float(np.sum(centres * values) / np.sum(values)),
            "gibuu_version": GIBUU_VERSION, "prepared": timestamp()}
    dump_json(info, out / "gibuu_prepare.json")
    return {"dir": out, **info}


def run_local(spec: GibuuSpec, channel, cfg, *, num_ensembles: int | None = None, seed: int | None = None,
              job_name: str = "local", timeout: int = 24 * 3600) -> Path:
    """Run GiBUU.x here (one job) in <cache>/jobs/<job_name>/; returns the job directory."""
    prep = prepare(spec, channel, cfg)
    paths = gibuu_paths(cfg.repo_root)
    if not paths["gibuu_x"].exists():
        raise FileNotFoundError(f"GiBUU binary not found at {paths['gibuu_x']} (pixi run build-gibuu)")
    job = prep["dir"] / "jobs" / job_name
    if job.exists():
        shutil.rmtree(job)
    job.mkdir(parents=True)
    card = write_card(spec, flux_file=prep["flux_file"], path_to_input=str(paths["path_to_input"]),
                      seed=spec.seed if seed is None else seed, out=job / "job.card", repo_root=cfg.repo_root, num_ensembles=num_ensembles)
    t0 = time.time()
    with open(card) as fi, open(job / "gibuu.log", "w") as fo:
        rc = subprocess.run([str(paths["gibuu_x"])], stdin=fi, stdout=fo, stderr=subprocess.STDOUT, cwd=str(job), timeout=timeout).returncode
    wall = time.time() - t0
    fe = job / "FinalEvents.dat"
    n_rows = sum(1 for ln in open(fe) if not ln.startswith("#")) if fe.exists() else 0
    dump_json({"job": job_name, "returncode": rc, "wall_s": round(wall, 1), "seed": spec.seed if seed is None else seed,
               "num_ensembles": int(num_ensembles or spec.num_ensembles), "n_finalevents_rows": n_rows, "card_sha256": sha256_file(card),
               "gibuu_x": str(paths["gibuu_x"]), "path_to_input": str(paths["path_to_input"]), "status": "ok" if rc == 0 and n_rows else "failed"},
              job / "manifest_job.json")
    if rc != 0 or not n_rows:
        raise RuntimeError(f"GiBUU failed (rc={rc}, {n_rows} FinalEvents rows); see {job / 'gibuu.log'}")
    return job


def _finalevents_path(job_dir: Path) -> Path:
    fe = job_dir / "FinalEvents.dat"
    if fe.exists():
        return fe
    gz = job_dir / "FinalEvents.dat.gz"
    if gz.exists():
        with gzip.open(gz, "rb") as fi, open(fe, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        return fe
    raise FileNotFoundError(f"{job_dir}: no FinalEvents.dat(.gz)")


def xsec_from_file(job_dir: Path) -> float | None:
    """Column 2 ("sum") of the last row of neutrino_absorption_cross_section_ALL.dat, in 1e-38 cm^2/nucleon."""
    p = job_dir / XSEC_FILE
    if not p.exists():
        return None
    rows = [ln.split() for ln in p.read_text().splitlines() if ln.strip() and not ln.startswith("#")]
    return float(rows[-1][1]) if rows else None


def read_job(job_dir: str | Path, spec: GibuuSpec, rel_tol: float = 1e-3) -> tuple[TruthTable, dict]:
    """One job's FinalEvents -> TruthTable (weights per run), cross-checked against the absorption file."""
    job_dir = Path(job_dir)
    t = read_finalevents(_finalevents_path(job_dir), target_Z=spec.target_Z, target_A=spec.target_A, n_runs=spec.num_runs)
    sw = float(t["weight"].sum())
    sig = xsec_from_file(job_dir)
    info = {"job": job_dir.name, "n_events": int(t.n), "sum_weights_1e-38cm2": sw, "xsec_file_1e-38cm2": sig}
    if sig is not None:
        info["rel_dev"] = (sw - sig) / sig if sig else None
        if sig and abs(sw - sig) > rel_tol * abs(sig):
            raise ValueError(f"{job_dir}: sum of weights {sw:.5g} != {XSEC_FILE} total {sig:.5g} (rel {abs(sw - sig) / sig:.2e})")
    return t, info


def merge_jobs(job_dirs: list, spec: GibuuSpec, channel, cfg, out_dir: str | Path | None = None, log=print) -> TruthTable:
    """K jobs -> one table whose weights sum to the mean cross section per nucleon (1e-38 cm^2 per unit weight)."""
    prep = prepare(spec, channel, cfg)
    out_dir = Path(out_dir) if out_dir else prep["dir"]
    tables, infos = [], []
    for jd in sorted(Path(j) for j in job_dirs):
        t, info = read_job(jd, spec)
        tables.append(t); infos.append(info)
        log(f"{jd.name}: {t.n} events, sigma {info['sum_weights_1e-38cm2']:.4f}e-38 cm^2/nucleon")
    K = len(tables)
    if K == 0:
        raise ValueError("no GiBUU jobs to merge")
    for t in tables:
        t.columns["weight"] = t["weight"] / K
    merged = TruthTable.concatenate(tables) if K > 1 else tables[0]
    sigmas = np.array([i["sum_weights_1e-38cm2"] for i in infos])
    sigma_avg = float(sigmas.mean()) * 1e-38
    merged.meta.update({
        "generator": GIBUU_VERSION, "gibuu_spec": spec.to_dict(), "frame": "beam", "has_geometry": False,
        "flux_source": prep["flux_source"], "flux_file": prep["flux_file"], "flux_file_sha256": prep["flux_file_sha256"],
        "template_sha256": prep["template_sha256"], "fingerprint": prep["fingerprint"],
        "n_jobs_merged": K, "jobs": infos, "sigma_flux_avg_per_nucleon_cm2": sigma_avg,
        "sigma_per_job_1e-38cm2": {"mean": float(sigmas.mean()), "std": float(sigmas.std(ddof=1)) if K > 1 else 0.0,
                                   "min": float(sigmas.min()), "max": float(sigmas.max())},
        "n_generated": int(merged.n), "source": str(out_dir),
        "norm": Normalization(kind="xsec_per_nucleon", xsec_per_unit_weight=1e-38,
                              notes=f"perweight/num_runs/K over K={K} jobs: weights sum to the mean flux-averaged sigma_CC per nucleon").to_dict(),
    })
    merged.save(out_dir / "truth.npz")
    dump_json({"fingerprint": prep["fingerprint"], "spec": spec.to_dict(), "flux_source": prep["flux_source"], "n_jobs_merged": K,
               "n_generated": int(merged.n), "sigma_flux_avg_per_nucleon_cm2": sigma_avg, "sigma_per_job_1e-38cm2": merged.meta["sigma_per_job_1e-38cm2"],
               "jobs": infos, "merged": timestamp(), "gibuu_version": GIBUU_VERSION}, out_dir / "gibuu_run.json")
    log(f"merged {K} job(s): {merged.n} events, sigma_CC = {sigma_avg:.4e} cm^2/nucleon -> {out_dir / 'truth.npz'}")
    return merged
