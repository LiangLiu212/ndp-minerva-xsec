"""GiBUU runner pieces that need no GiBUU binary: flux rebinning, card writing, the energy-scan allocation,
and the merge of FinalEvents jobs (weights / K; flux-fraction weights for energy points; the cross-section check)."""
import json
import tempfile
from pathlib import Path

import numpy as np

from ndp.theory import flux as fluxmod
from ndp.theory.gibuu import GibuuSpec, energy_allocation, merge_energy_scan, merge_jobs, read_job, stratum_specs, substitutions, write_card
from _helpers import site

ROOT = Path(__file__).resolve().parents[1]

import ndp.theory.gibuu as gibuu_mod
gibuu_mod.CACHE_ROOT_OVERRIDE = Path(tempfile.mkdtemp())      # keep the tests' sample caches out of runs/_generator_cache

FINALEVENTS_HEADER = "# 1:Run 2:Event 3:ID 4:Charge 5:perweight 6:position(1) 7:position(2) 8:position(3) 9:momentum(0) 10:momentum(1) 11:momentum(2) 12:momentum(3) 13:history 14:production_ID 15:enu\n"


def _finalevents(path: Path, events: list, run: int = 1):
    """events: list of (perweight, production_id, enu, [(gibuu_id, charge, E, px, py, pz), ...]) — a muon row is added to each."""
    lines = [FINALEVENTS_HEADER]
    for i, (w, pid, enu, parts) in enumerate(events, start=1):
        rows = [(902, -1, 0.5 * enu, 0.05, 0.0, 0.5 * enu)] + list(parts)
        for gid, q, E, px, py, pz in rows:
            lines.append(f"{run} {i} {gid} {q} {w:.6e} 0.0 0.0 0.0 {E:.6f} {px:.6f} {py:.6f} {pz:.6f} 0 {pid} {enu:.4f}\n")
    path.write_text("".join(lines))


def _xsec_file(path: Path, total: float):
    path.write_text("# 1:var 2:sum 3:QE 4:Delta 5:highRES 6:1pi 7:DIS 8:2p2h-QE  9:2p2h-Delta 10:2pi\n"
                    f"   6.0000    {total:.5f}   {total:.5f}   0.0 0.0 0.0 0.0 0.0 0.0 0.0\n")


def test_rebin_uniform_conserves_the_integral_and_gibuu_flux_file_is_equidistant():
    edges = np.array([0.0, 0.5, 1.0, 1.2, 1.4, 2.0, 3.5, 5.0])
    dens = np.array([1.0, 2.0, 3.0, 3.5, 2.5, 1.0, 0.2])
    c, v = fluxmod.rebin_uniform(edges, dens, width=0.5, e_max=5.0)
    assert len(c) == 10 and np.allclose(c[:3], [0.25, 0.75, 1.25])
    assert np.isclose(np.sum(v) * 0.5, fluxmod.integrated_flux(edges, dens, 0.0, 5.0))
    assert np.isclose(v[0], 1.0) and np.isclose(v[2], (0.2 * 3.0 + 0.2 * 3.5 + 0.1 * 2.5) / 0.5)
    p = Path(tempfile.mkdtemp()) / "flux.dat"
    fluxmod.write_gibuu_flux(p, c, v, source="toy")
    rows = [ln.split() for ln in p.read_text().splitlines() if not ln.startswith("#")]
    assert len(rows) == 10 and np.allclose(np.diff([float(r[0]) for r in rows]), 0.5)


def test_card_has_no_placeholder_left_and_reflects_the_spec():
    spec = GibuuSpec(num_ensembles=250, seed=7, medium_switch=False, equal_weights_mode=2, equal_weights_max=600.0)
    out = Path(tempfile.mkdtemp()) / "job.card"
    write_card(spec, flux_file="/x/flux.dat", path_to_input="/x/buuinput", seed=7, out=out, repo_root=ROOT)
    txt = out.read_text()
    assert "__" not in txt.replace("__file__", "")
    assert "numEnsembles=250" in txt and "SEED=7" in txt and "mediumSwitch=.false." in txt and "equalWeights_Max  = 600.0" in txt
    assert "nuXsectionMode  = 16" in txt and "nuExp           = 99" in txt
    scan = GibuuSpec(mode="energy_scan")
    out2 = Path(tempfile.mkdtemp()) / "scan.card"
    write_card(scan, flux_file="/x/flux.dat", path_to_input="/x/buuinput", seed=1, out=out2, repo_root=ROOT, energy=6.25, num_ensembles=100)
    t2 = out2.read_text()
    assert "nuXsectionMode  = 0 " in t2 and "nuExp           = 0 " in t2 and "enu=6.25" in t2
    try:
        write_card(scan, flux_file="f", path_to_input="p", seed=1, out=out2, repo_root=ROOT)
        raise AssertionError("energy_scan without an energy must be refused")
    except ValueError:
        pass
    d = substitutions(GibuuSpec(), flux_file="f", path_to_input="p", seed=3)
    assert d["__INCLUDE_QE__"] == "T" and d["__SEED__"] == "3"


def test_strata_are_disjoint_and_override_the_base_spec():
    base = {"target": {"Z": 6, "A": 12}, "card": {"num_ensembles": 1000}, "seed": 5,
            "strata": {"qe": {"include": {k: False for k in ("DELTA", "RES", "1pi", "DIS", "2p2hQE", "2pi")}, "card": {"equal_weights_mode": 2, "equal_weights_max": 1500}},
                       "rest": {"include": {"QE": False}, "card": {"equal_weights_mode": 2, "equal_weights_max": 600}}}}
    st = stratum_specs(base)
    assert set(st) == {"qe", "rest"} and st["qe"].include["QE"] and not st["rest"].include["QE"] and st["rest"].num_ensembles == 1000
    assert st["qe"].equal_weights_max == 1500 and st["rest"].equal_weights_max == 600
    bad = dict(base); bad["strata"] = {"a": {}, "b": {}}
    try:
        stratum_specs(bad); raise AssertionError("overlapping strata must be refused")
    except ValueError:
        pass


def test_energy_allocation_covers_the_flux_and_respects_the_minimum():
    cfg = site()
    from ndp.channels import load_channel
    fl = fluxmod.load_channel_flux(load_channel("minerva_me_ccqelike_1mu1p"), cfg.repo_root)
    spec = GibuuSpec(mode="energy_scan", total_ensembles=50000, energy_e_min_gev=2.0, energy_e_max_gev=60.0, energy_step_gev=0.5)
    pts = energy_allocation(spec, fl)
    assert len(pts) == 116 and all(q["n_ensembles"] >= 100 for q in pts)
    from ndp.theory.gibuu import energy_jobs
    jobs = energy_jobs(pts, 5)
    assert len(jobs) == sum(q["n_jobs"] for q in pts) and all(j["n_ensembles"] <= spec.max_ensembles_per_job for j in jobs)
    assert [j["j"] for j in jobs] == list(range(len(jobs))) and jobs[0]["seed"] == 5
    frac = np.array([q["flux_fraction"] for q in pts])
    assert 0.8 < frac.sum() < 1.0                       # 2-60 GeV holds most of the 0-100 GeV flux
    assert abs(frac.sum() - fluxmod.integrated_flux(fl["edges"], fl["density_cm2_pot_gev"], 2.0, 60.0) / fluxmod.integrated_flux(fl["edges"], fl["density_cm2_pot_gev"], 0.0, 100.0)) < 1e-9
    peak = max(pts, key=lambda q: q["n_ensembles"])
    assert 3.0 < peak["energy"] < 8.0                   # the ME flux peak


def test_merge_jobs_prescales_by_the_job_count_and_checks_the_cross_section():
    cfg = site()
    from ndp.channels import load_channel
    ch = load_channel("minerva_me_ccqelike_1mu1p")
    spec = GibuuSpec(num_ensembles=100)
    root = Path(tempfile.mkdtemp())
    jobs = []
    for k, (w, tot) in enumerate(((0.5, 1.0), (0.25, 1.0))):
        d = root / f"job{k}"; d.mkdir()
        _finalevents(d / "FinalEvents.dat", [(w, 1, 6.0, [(1, 1, 1.2, 0.1, 0.0, 0.7)]), (w, 34, 6.0, [(101, 1, 0.4, 0.0, 0.1, 0.3)])] * 2)
        _xsec_file(d / "neutrino_absorption_cross_section_ALL.dat", 4 * w)
        jobs.append(d)
    t = merge_jobs(jobs, spec, ch, cfg, out_dir=root / "merged", log=lambda *a: None)
    assert t.n == 8 and np.isclose(t["weight"].sum(), (4 * 0.5 + 4 * 0.25) / 2)        # mean of the two jobs' cross sections
    assert t.meta["norm"]["xsec_per_unit_weight"] == 1e-38 and t.meta["frame"] == "beam" and t.meta["n_jobs_merged"] == 2
    assert (root / "merged" / "truth.npz").exists() and (root / "merged" / "gibuu_run.json").exists()
    bad = root / "bad"; bad.mkdir()
    _finalevents(bad / "FinalEvents.dat", [(0.5, 1, 6.0, [(1, 1, 1.2, 0.1, 0.0, 0.7)])])
    _xsec_file(bad / "neutrino_absorption_cross_section_ALL.dat", 0.7)                   # disagrees with the summed weights (0.5)
    try:
        read_job(bad, spec); raise AssertionError("cross-section mismatch must be refused")
    except ValueError:
        pass


def test_merge_energy_scan_weights_by_flux_fraction():
    cfg = site()
    from ndp.channels import load_channel
    ch = load_channel("minerva_me_ccqelike_1mu1p")
    spec = GibuuSpec(mode="energy_scan", total_ensembles=20000, energy_e_min_gev=2.0, energy_e_max_gev=4.0, energy_step_gev=1.0)   # points 2.5, 3.5
    root = Path(tempfile.mkdtemp())
    fl = fluxmod.load_channel_flux(ch, cfg.repo_root)
    pts = energy_allocation(spec, fl)
    jobs = []
    for q in pts:
        d = root / f"E{q['k']}"; d.mkdir()
        _finalevents(d / "FinalEvents.dat", [(1.0, 1, q["energy"], [(1, 1, 1.2, 0.1, 0.0, 0.7)])] * 3)
        _xsec_file(d / "neutrino_absorption_cross_section_ALL.dat", 3.0)
        (d / "manifest_job.json").write_text(json.dumps({"energy_point": q}))
        jobs.append(d)
    t = merge_energy_scan(jobs, spec, ch, cfg, out_dir=root / "merged", log=lambda *a: None)
    assert t.n == 6
    expected = sum(3.0 * q["flux_fraction"] for q in pts)
    assert np.isclose(t["weight"].sum(), expected) and np.isclose(t.meta["sigma_flux_avg_per_nucleon_cm2"], expected * 1e-38)
    assert t.meta["mode"] == "energy_scan" and t.meta["energy_points_missing"] == [] and "gibuu_energy_k" in t
