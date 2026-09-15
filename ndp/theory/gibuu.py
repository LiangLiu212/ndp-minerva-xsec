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
    mode: str = "flux_mc"                    # "flux_mc": nuXsectionMode 16 (weighted events); "energy_scan": nuXsectionMode 0 per energy of a flux-weighted grid
    energy_e_min_gev: float = 2.0            # energy_scan grid: bin centres e_min + step/2, ..., < e_max
    energy_e_max_gev: float = 60.0
    energy_step_gev: float = 0.5
    total_ensembles: int = 85000             # energy_scan: ensembles summed over the grid, allocated ~ flux x E (at least min_ensembles each)
    min_ensembles: int = 100
    max_ensembles_per_job: int = 1500        # energy_scan: a point with more ensembles is split over several jobs (merged as n_jobs_k)
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
                  flux_e_max_gev=float(fl.get("e_max_gev", 100.0)), mode=str(p.get("mode", "flux_mc")),
                  energy_e_min_gev=float((p.get("energy_grid") or {}).get("e_min_gev", 2.0)), energy_e_max_gev=float((p.get("energy_grid") or {}).get("e_max_gev", 60.0)),
                  energy_step_gev=float((p.get("energy_grid") or {}).get("step_gev", 0.5)), total_ensembles=int(p.get("total_ensembles", 85000)),
                  min_ensembles=int(p.get("min_ensembles", 100)), max_ensembles_per_job=int(p.get("max_ensembles_per_job", 1500)),
                  n_jobs=int(p.get("n_jobs", 100)),
                  seed=int(p.get("seed", 20260915)), template=str(p.get("template", DEFAULT_TEMPLATE)))
        if not (-2 ** 31 < kw["seed"] + kw["n_jobs"] < 2 ** 31):
            raise ValueError("GiBUU seeds must fit a 32-bit Fortran integer")
        if kw["equal_weights_mode"] == 2 and kw["equal_weights_max"] <= 0:
            raise ValueError("equal_weights_mode 2 needs a positive equal_weights_max (from a mode-1 pilot run)")
        if kw["mode"] not in ("flux_mc", "energy_scan"):
            raise ValueError(f"unknown GiBUU mode {kw['mode']!r} (flux_mc | energy_scan)")
        return GibuuSpec(**kw)

    def energy_centres(self) -> np.ndarray:
        """Bin centres of the energy_scan grid."""
        n = int(round((self.energy_e_max_gev - self.energy_e_min_gev) / self.energy_step_gev))
        return self.energy_e_min_gev + self.energy_step_gev * (np.arange(n) + 0.5)

    def to_dict(self) -> dict:
        return asdict(self)

    def fingerprint(self, flux_source: str, template_sha: str) -> str:
        return sha256_text(json.dumps(self.to_dict(), sort_keys=True) + "|" + flux_source + "|" + template_sha)[:16]

    def events_per_job(self) -> int:
        """Rule of thumb from GiBUU's card comments: events ~ A x numEnsembles x num_runs_SameEnergy."""
        return self.target_A * self.num_ensembles * self.num_runs


def stratum_specs(params: dict) -> dict:
    """{stratum name: GibuuSpec}. A model spec may split the generation into strata (e.g. `qe` and `rest`,
    each with its own `include` flags, equal-weights ceiling, ensembles, jobs and seed) whose weights are
    normalised per stratum and simply concatenated: sigma_total = sum of the strata's cross sections.
    Without a `strata:` block the whole spec is one stratum named `all`."""
    strata = params.get("strata")
    if not strata:
        return {"all": GibuuSpec.from_params(params)}
    out = {}
    for name, over in strata.items():
        merged = {k: v for k, v in params.items() if k != "strata"}
        over = over or {}
        for k, v in over.items():
            if isinstance(v, dict) and isinstance(merged.get(k), dict):
                merged[k] = {**merged[k], **v}
            else:
                merged[k] = v
        out[str(name)] = GibuuSpec.from_params(merged)
    incl = [tuple(sorted(k for k, v in g.include.items() if v)) for g in out.values()]
    if len(set(incl)) != len(incl):
        raise ValueError("GiBUU strata must include disjoint sets of channels")
    for a in range(len(incl)):
        for b in range(a + 1, len(incl)):
            if set(incl[a]) & set(incl[b]):
                raise ValueError(f"GiBUU strata overlap in channels {sorted(set(incl[a]) & set(incl[b]))}")
    return out


def concatenate_strata(tables: dict, log=print) -> TruthTable:
    """Strata tables (each normalised to its own cross section) -> one table; weights are NOT rescaled."""
    names = list(tables)
    if len(names) == 1:
        return tables[names[0]]
    for t in tables.values():
        t.columns["stratum"] = None
    parts = []
    for i, n in enumerate(names):
        t = tables[n]
        t.columns.pop("stratum", None)
        t.columns["gibuu_stratum"] = np.full(t.n, i, dtype=np.int64)
        parts.append(t)
    merged = TruthTable.concatenate(parts)
    sig = {n: float(tables[n]["weight"].sum()) * 1e-38 for n in names}
    merged.meta.update({"gibuu_strata": {n: {"index": i, "sigma_flux_avg_per_nucleon_cm2": sig[n], "n_events": int(tables[n].n),
                                            "fingerprint": tables[n].meta.get("fingerprint"), "gibuu_spec": tables[n].meta.get("gibuu_spec")}
                                        for i, n in enumerate(names)},
                        "sigma_flux_avg_per_nucleon_cm2": float(sum(sig.values())), "n_generated": int(merged.n),
                        "source": "; ".join(str(tables[n].meta.get("source")) for n in names)})
    log(f"strata {names}: sigma_CC = {sum(sig.values()):.4e} cm^2/nucleon over {merged.n} events")
    return merged


def _fb(b: bool) -> str:
    return "T" if b else "F"


def substitutions(spec: GibuuSpec, *, flux_file: str, path_to_input: str, seed: int, num_ensembles: int | None = None,
                  energy: float | str | None = None) -> dict:
    inc = spec.include
    scan = spec.mode == "energy_scan"
    enu = energy if energy is not None else (1.0 if not scan else "__ENU__")
    return {"__PROCESS_ID__": str(spec.process_ID), "__FLAVOR_ID__": str(spec.flavor_ID), "__FLUX_FILE__": str(flux_file),
            "__NU_XSECTION_MODE__": "0" if scan else "16", "__NU_EXP__": "0" if scan else "99",
            "__ENU__": enu if isinstance(enu, str) else repr(float(enu)), "__DELTA_ENU__": "0.0",
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
               repo_root: str | Path | None = None, num_ensembles: int | None = None, energy: float | None = None) -> Path:
    """Fill the template with the spec; refuses to leave a placeholder behind."""
    tmpl = Path(spec.template)
    if not tmpl.is_absolute():
        tmpl = Path(repo_root or Path(__file__).resolve().parents[2]) / tmpl
    text = tmpl.read_text()
    if spec.mode == "energy_scan" and energy is None:
        raise ValueError("energy_scan cards need the energy of the point")
    for k, v in substitutions(spec, flux_file=flux_file, path_to_input=path_to_input, seed=seed, num_ensembles=num_ensembles, energy=energy).items():
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


def energy_allocation(spec: GibuuSpec, fl: dict) -> list[dict]:
    """energy_scan points: {k, energy, flux_fraction, n_ensembles}. flux_fraction = the flux integral of the point's
    bin over the FULL flux integral (0-100 GeV), so events of a point weighted by perweight x flux_fraction sum to that
    bin's share of the flux-averaged cross section. Ensembles are allocated ~ flux_fraction x E (sigma_CC roughly linear
    in E above 2 GeV) so that the merged event weights are as uniform as the allocation allows; every point gets at
    least min_ensembles (GiBUU refuses fewer than 100)."""
    E = spec.energy_centres()
    half = spec.energy_step_gev / 2
    phi_tot = fluxmod.integrated_flux(fl["edges"], fl["density_cm2_pot_gev"], 0.0, 100.0)
    frac = np.array([fluxmod.integrated_flux(fl["edges"], fl["density_cm2_pot_gev"], e - half, e + half) / phi_tot for e in E])
    p = frac * E
    p = p / p.sum()
    n = np.maximum(np.round(p * spec.total_ensembles).astype(int), spec.min_ensembles)
    out = []
    for k in range(len(E)):
        n_jobs = int(np.ceil(n[k] / spec.max_ensembles_per_job))
        per = int(np.ceil(n[k] / n_jobs))
        out.append({"k": int(k), "energy": float(E[k]), "flux_fraction": float(frac[k]), "n_ensembles": int(per * n_jobs),
                    "n_jobs": n_jobs, "ensembles_per_job": per})
    return out


def energy_jobs(points: list[dict], seed: int) -> list[dict]:
    """One grid job per (point, split): {j, k, energy, flux_fraction, n_ensembles, seed}; seed = base + j."""
    jobs = []
    for q in points:
        for _ in range(q["n_jobs"]):
            j = len(jobs)
            jobs.append({"j": j, "k": q["k"], "energy": q["energy"], "flux_fraction": q["flux_fraction"], "n_ensembles": q["ensembles_per_job"], "seed": seed + j})
    return jobs


def prepare(spec: GibuuSpec, channel, cfg) -> dict:
    """The cache directory for this spec + channel flux: writes the GiBUU flux file (and, for energy_scan, the
    energy table) and records the fingerprint."""
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
    if spec.mode == "energy_scan":
        pts = energy_allocation(spec, fl)
        jobs = energy_jobs(pts, spec.seed)
        (out / "energies.txt").write_text("# j k energy_gev flux_fraction n_ensembles seed   (one line per grid job; j = process number)\n" +
                                          "".join(f"{q['j']} {q['k']} {q['energy']:.4f} {q['flux_fraction']:.6e} {q['n_ensembles']} {q['seed']}\n" for q in jobs))
        info.update({"energy_points": pts, "energy_jobs": jobs, "n_energy_points": len(pts), "n_energy_jobs": len(jobs),
                     "flux_fraction_covered": float(sum(q["flux_fraction"] for q in pts)), "energies_file": str(out / "energies.txt")})
    dump_json(info, out / "gibuu_prepare.json")
    return {"dir": out, **info}


def run_local(spec: GibuuSpec, channel, cfg, *, num_ensembles: int | None = None, seed: int | None = None,
              job_name: str = "local", timeout: int = 24 * 3600, energy_point: int | None = None) -> Path:
    """Run GiBUU.x here (one job) in <cache>/jobs/<job_name>/; returns the job directory. For energy_scan specs
    `energy_point` is the index k of the energy table (default: the point nearest the flux mean)."""
    prep = prepare(spec, channel, cfg)
    paths = gibuu_paths(cfg.repo_root)
    if not paths["gibuu_x"].exists():
        raise FileNotFoundError(f"GiBUU binary not found at {paths['gibuu_x']} (pixi run build-gibuu)")
    job = prep["dir"] / "jobs" / job_name
    if job.exists():
        shutil.rmtree(job)
    job.mkdir(parents=True)
    energy = None; point = None
    if spec.mode == "energy_scan":
        pts = prep["energy_points"]
        k = energy_point if energy_point is not None else int(np.argmin([abs(q["energy"] - prep["flux_mean_energy_gev"]) for q in pts]))
        point = {kk: v for kk, v in pts[k].items() if kk in ("k", "energy", "flux_fraction", "n_ensembles")}; energy = point["energy"]
        num_ensembles = num_ensembles or pts[k]["ensembles_per_job"]
        seed = spec.seed + k if seed is None else seed
    card = write_card(spec, flux_file=prep["flux_file"], path_to_input=str(paths["path_to_input"]),
                      seed=spec.seed if seed is None else seed, out=job / "job.card", repo_root=cfg.repo_root, num_ensembles=num_ensembles, energy=energy)
    t0 = time.time()
    with open(card) as fi, open(job / "gibuu.log", "w") as fo:
        rc = subprocess.run([str(paths["gibuu_x"])], stdin=fi, stdout=fo, stderr=subprocess.STDOUT, cwd=str(job), timeout=timeout).returncode
    wall = time.time() - t0
    fe = job / "FinalEvents.dat"
    n_rows = sum(1 for ln in open(fe) if not ln.startswith("#")) if fe.exists() else 0
    dump_json({"job": job_name, "returncode": rc, "wall_s": round(wall, 1), "seed": spec.seed if seed is None else seed,
               "num_ensembles": int(num_ensembles or spec.num_ensembles), "n_finalevents_rows": n_rows, "card_sha256": sha256_file(card),
               "mode": spec.mode, "energy_point": point,
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


def _job_energy_point(job_dir: Path) -> dict | None:
    for name in ("manifest_job.json",) + tuple(p.name for p in job_dir.glob("manifest_*.json")):
        q = job_dir / name
        if q.exists():
            m = json.loads(q.read_text())
            if m.get("energy_point"):
                return m["energy_point"]
            if m.get("energy_k") is not None:
                return {"k": int(m["energy_k"]), "energy": float(m["energy_gev"]), "flux_fraction": float(m["flux_fraction"]),
                        "n_ensembles": int(m.get("num_ensembles", 0))}
    return None


def merge_energy_scan(job_dirs: list, spec: GibuuSpec, channel, cfg, out_dir: str | Path | None = None, log=print) -> TruthTable:
    """energy_scan jobs (one fixed energy each, unit-like weights) -> one flux-averaged table.

    An event of energy point k gets weight perweight x flux_fraction_k / n_jobs_k, so the merged weights sum to
    sum_k f_k sigma_CC(E_k) = the flux-averaged cross section per nucleon over the covered energy range
    (1e-38 cm^2 per unit weight); the per-point cross sections are kept in the meta as sigma_CC(E).
    """
    prep = prepare(spec, channel, cfg)
    out_dir = Path(out_dir) if out_dir else prep["dir"]
    out_dir.mkdir(parents=True, exist_ok=True)
    by_k: dict[int, list] = {}
    for jd in sorted(Path(j) for j in job_dirs):
        pt = _job_energy_point(jd)
        if pt is None:
            raise ValueError(f"{jd}: no energy point in its manifest (not an energy_scan job?)")
        by_k.setdefault(int(pt["k"]), []).append((jd, pt))
    tables, points, missing = [], [], []
    for q in prep["energy_points"]:
        jobs = by_k.get(q["k"], [])
        if not jobs:
            missing.append(q["k"]); continue
        sig = []
        for jd, pt in jobs:
            t, info = read_job(jd, spec)
            t.columns["weight"] = t["weight"] * (q["flux_fraction"] / len(jobs))
            t.columns["gibuu_energy_k"] = np.full(t.n, q["k"], dtype=np.int64)
            tables.append(t); sig.append(info["sum_weights_1e-38cm2"])
        points.append({**{kk: v for kk, v in q.items() if kk in ("k", "energy", "flux_fraction", "n_ensembles")}, "n_jobs_planned": q["n_jobs"],
                       "n_jobs": len(jobs), "sigma_cc_1e-38cm2": float(np.mean(sig)), "n_events": int(sum(tb.n for tb in tables[-len(jobs):]))})
        log(f"E = {q['energy']:.2f} GeV: {len(jobs)} job(s), {points[-1]['n_events']} events, sigma_CC {points[-1]['sigma_cc_1e-38cm2']:.4f}e-38, flux fraction {q['flux_fraction']:.4e}")
    if not tables:
        raise ValueError("no energy_scan jobs to merge")
    merged = TruthTable.concatenate(tables) if len(tables) > 1 else tables[0]
    sigma_avg = float(merged["weight"].sum()) * 1e-38
    covered = float(sum(q["flux_fraction"] for q in points))
    merged.meta.update({
        "generator": GIBUU_VERSION, "gibuu_spec": spec.to_dict(), "frame": "beam", "has_geometry": False,
        "flux_source": prep["flux_source"], "template_sha256": prep["template_sha256"], "fingerprint": prep["fingerprint"],
        "mode": "energy_scan", "energy_points": points, "energy_points_missing": missing, "flux_fraction_covered": covered,
        "n_jobs_merged": int(sum(q["n_jobs"] for q in points)), "sigma_flux_avg_per_nucleon_cm2": sigma_avg,
        "n_generated": int(merged.n), "source": str(out_dir),
        "norm": Normalization(kind="xsec_per_nucleon", xsec_per_unit_weight=1e-38,
                              notes=f"energy_scan: perweight x flux_fraction / n_jobs per energy point; weights sum to the flux-averaged sigma_CC over "
                                    f"{covered:.4f} of the 0-100 GeV flux ({len(points)} points)").to_dict(),
    })
    merged.save(out_dir / "truth.npz")
    dump_json({"fingerprint": prep["fingerprint"], "spec": spec.to_dict(), "flux_source": prep["flux_source"], "mode": "energy_scan",
               "energy_points": points, "energy_points_missing": missing, "flux_fraction_covered": covered, "n_generated": int(merged.n),
               "sigma_flux_avg_per_nucleon_cm2": sigma_avg, "merged": timestamp(), "gibuu_version": GIBUU_VERSION}, out_dir / "gibuu_run.json")
    log(f"merged {len(points)} energy points ({sum(q['n_jobs'] for q in points)} jobs, {len(missing)} missing): {merged.n} events, "
        f"sigma_CC = {sigma_avg:.4e} cm^2/nucleon over {covered:.4f} of the flux -> {out_dir / 'truth.npz'}")
    return merged


def merge_jobs(job_dirs: list, spec: GibuuSpec, channel, cfg, out_dir: str | Path | None = None, log=print) -> TruthTable:
    """K jobs -> one table whose weights sum to the mean cross section per nucleon (1e-38 cm^2 per unit weight)."""
    if spec.mode == "energy_scan":
        return merge_energy_scan(job_dirs, spec, channel, cfg, out_dir=out_dir, log=log)
    prep = prepare(spec, channel, cfg)
    out_dir = Path(out_dir) if out_dir else prep["dir"]
    out_dir.mkdir(parents=True, exist_ok=True)
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
