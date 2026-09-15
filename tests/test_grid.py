"""Grid processing pieces: URL handling, derived-column skims, POT-summing concatenation, vectorised
select, worklist slicing / resubmission logic, playlist merge, and the products loader fallback."""
import json
import os
import tempfile
from pathlib import Path

import numpy as np
from _helpers import legacy_channel, site, toy_truth, have_mc_cache, have_data_file, skip
from ndp.channels import load_channel
from ndp.channels import signal as sig
from ndp.channels import observables as obs
from ndp.events import TruthTable, Normalization
from ndp.io import is_url
from ndp.adapters.minerva_anatuple import cache_tag, derive_fs_columns, skim_truth, skim_mask, DEFAULT_SKIM_BOX
from ndp.grid.process_file import swap_door

URL = "root://fndcadoor.fnal.gov:1095//pnfs/fnal.gov/usr/minerva/persistent/OpenData/MediumEnergy_FHC/MC/StandardMC/Playlist1A/MasterAnaDev_mc_AnaTuple_run00110040_Playlist.root"


def _fs_truth(n=4000, seed=3):
    """Toy numu CC sample with a realistic final-state mix (protons in and out of the window, pions, photons)."""
    rng = np.random.default_rng(seed)
    t = toy_truth(n, seed=seed, with_fs=True)
    cols = dict(t.columns)
    # replace the toy's collinear FS with random directions and a richer PDG mix
    off = cols["fs_offsets"]; tot = int(off[-1])
    pdg = rng.choice([2212, 2112, 211, -211, 111, 22, 3122, -2212], size=tot, p=[0.45, 0.25, 0.08, 0.06, 0.06, 0.06, 0.02, 0.02])
    mass = np.select([pdg == 2212, pdg == 2112, np.abs(pdg) == 211, pdg == 111, pdg == 22, pdg == 3122, pdg == -2212],
                     [0.93827, 0.93957, 0.13957, 0.13498, 0.0, 1.11568, 0.93827], 0.5)
    p = np.where(pdg == 22, rng.exponential(0.02, tot), rng.uniform(0.2, 1.4, tot))
    th = np.deg2rad(rng.uniform(0, 100, tot)); ph = rng.uniform(0, 2 * np.pi, tot)
    cols.update({"fs_pdg": pdg, "fs_E": np.sqrt(p * p + mass * mass), "fs_px": p * np.sin(th) * np.cos(ph),
                 "fs_py": p * np.sin(th) * np.sin(ph), "fs_pz": p * np.cos(th)})
    cols["vtx_z"] = rng.uniform(5500, 8800, n); cols["vtx_x"] = rng.uniform(-950, 950, n); cols["vtx_y"] = rng.uniform(-950, 950, n)
    cols["current"] = rng.choice([1, 2], size=n, p=[0.85, 0.15])
    meta = dict(t.meta); meta["norm"] = Normalization(kind="pot", pot=1.0e18).to_dict(); meta["has_geometry"] = True
    return TruthTable(cols, meta)


def test_urls_never_go_through_path():
    from ndp.grid.campaign import xrootd_path
    assert xrootd_path("/pnfs/dune/scratch/users/x/y") == "/pnfs/fs/usr/dune/scratch/users/x/y"
    assert xrootd_path("/pnfs/fs/usr/dune/z") == "/pnfs/fs/usr/dune/z"
    assert is_url(URL) and not is_url("/exp/dune/data/x.root")
    assert cache_tag(URL) == "mc110040"
    assert swap_door(URL, "fndca1.fnal.gov:1095").startswith("root://fndca1.fnal.gov:1095//pnfs/") and swap_door("/a/b.root", "x") == "/a/b.root"


def test_derived_skim_reproduces_the_fs_path():
    ch = load_channel("minerva_me_ccqelike_1mu1p")
    t = _fs_truth()
    d = derive_fs_columns(t, ch)
    assert d.meta["derived"]["frame"] == ch.frame and "signal_minerva_ccqelike_1mu1p" in d
    assert np.array_equal(d["signal_minerva_ccqelike_1mu1p"], ch.is_signal(t))
    ka = skim_truth(d, keep_all_rows=True)
    assert ka.n == t.n and not ka.has_fs and ka.meta["skim"]["keep_all_rows"]
    s = skim_truth(d, DEFAULT_SKIM_BOX)
    m = skim_mask(d, DEFAULT_SKIM_BOX)
    assert not s.has_fs and s.n == int(m.sum()) and s.meta["skim"]["n_after"] == s.n and (s["current"] == 1).all()
    full = t.select(m)
    # every truth-level quantity the channel uses agrees between the skim and the FS computation
    assert np.array_equal(ch.is_signal(s), ch.is_signal(full))
    for name in ("proton_p", "proton_theta", "proton_pT", "n_protons_in_window", "dpT", "dpTx", "dpTy", "dalphaT", "dphiT", "dpL", "pn", "E_avail", "E_had_fs"):
        a, b = ch.evaluate(name, s), ch.evaluate(name, full)
        assert np.allclose(a, b, equal_nan=True), name
    steps_s, steps_f = dict(ch.signal_cutflow(s)), dict(ch.signal_cutflow(full))
    assert all(np.array_equal(steps_s[k], steps_f[k]) for k in steps_f)
    # a skim derived for another window / frame is refused, not silently reused
    s.meta["derived"]["proton_window"]["p_max_gev"] = 1.3
    try:
        ch.evaluate("proton_p", s); raise AssertionError("window mismatch must raise")
    except ValueError:
        pass


def test_concatenate_sums_pot_and_select_is_vectorised():
    a = _fs_truth(300, seed=1); b = _fs_truth(500, seed=2)
    b.meta["norm"] = Normalization(kind="pot", pot=3.0e18).to_dict()
    c = TruthTable.concatenate([a, b])
    assert abs(c.norm.pot - 4.0e18) < 1 and len(c.meta["sources"]) == 2 and c.n == 800
    shape = _fs_truth(50, seed=4); shape.meta["norm"] = Normalization(kind="shape").to_dict()
    try:
        TruthTable.concatenate([a, shape]); raise AssertionError("mixed kinds must raise")
    except ValueError:
        pass
    # vectorised select == explicit per-event loop
    mask = np.random.default_rng(0).random(a.n) < 0.4
    sel = a.select(mask)
    off = a["fs_offsets"]; idx = np.concatenate([np.arange(off[i], off[i + 1]) for i in np.flatnonzero(mask)])
    assert np.array_equal(sel["fs_pdg"], a["fs_pdg"][idx]) and np.array_equal(sel["fs_E"], a["fs_E"][idx])
    assert sel.n == int(mask.sum()) and sel["fs_offsets"][-1] == len(idx)


def test_worklist_slicing_and_resubmit(tmp_path=None):
    from ndp.grid import campaign as cp
    tmp = Path(tempfile.mkdtemp())
    old = cp.CAMPAIGNS; cp.CAMPAIGNS = tmp
    try:
        c = {"name": "t", "worklists": {"FHC_mc_1A": {"file": str(tmp / "w.txt"), "kind": "mc", "playlist": "1A", "n_files": 5,
                                                       "files_per_process": 2, "n_processes": 3, "pnfs_out": "/pnfs/x"}},
             "files": {f"mc{i}": {"tag": f"mc{i}", "url": f"root://h//f{i}.root", "beam": "FHC", "kind": "mc", "playlist": "1A",
                                  "worklist": "FHC_mc_1A", "line": i, "process": i // 2, "status": s, "attempts": 0}
                       for i, s in enumerate(["done", "failed", "done", "incomplete", "planned"])}}
        cp.save_campaign(c)
        out = cp.resubmit("t", log=lambda *a: None)
        lines = Path(out["FHC_mc_1A"]).read_text().split()
        assert lines == ["root://h//f1.root", "root://h//f3.root", "root://h//f4.root"]
        c2 = cp.load_campaign("t")
        assert c2["files"]["mc1"]["status"] == "planned" and c2["files"]["mc0"]["status"] == "done"
        # the worker's slice rule: process p takes lines [p*n, (p+1)*n) (1-based sed range p*n+1 .. p*n+n)
        n = 2; per = {p: lines[p * n:(p + 1) * n] for p in range(2)}
        assert per == {0: lines[0:2], 1: lines[2:3]}
    finally:
        cp.CAMPAIGNS = old


def test_merge_playlist_from_per_file_products():
    from ndp.products import merge_playlist, load_reco_npz
    ch = load_channel("minerva_me_ccqelike_1mu1p")
    root = Path(tempfile.mkdtemp())
    pots = {"mc1": 1.5e18, "mc2": 2.5e18}
    for tag, pot in pots.items():
        d = root / "FHC" / "1A" / "files" / tag; d.mkdir(parents=True)
        t = _fs_truth(200, seed=int(tag[-1])); t.meta["norm"] = Normalization(kind="pot", pot=pot).to_dict()
        skim_truth(derive_fs_columns(t, ch)).save(d / f"truth_{tag}_skim.npz")
        skim_truth(derive_fs_columns(t, ch), keep_all_rows=True).save(d / f"reco_{tag}_truthcols_skim.npz")
        meta = {"cache_version": 3, "is_mc": True, "pot": {"pot_used": pot}}
        np.savez_compressed(d / f"reco_{tag}.npz", __meta__=json.dumps(meta), passed=np.ones(200, bool), reco_p=np.arange(200.0))
        (d / f"manifest_{tag}.json").write_text(json.dumps({"tag": tag, "kind": "mc", "status": "ok", "pot_used": pot, "pot_total": pot * 1.1, "n_reco": 10,
                                                            "selection_cutflow": {"all_candidates": 200, "IsoBlobs": 3},
                                                            "signal_cutflow": {"numu_cc_muon": {"n": 100, "n_fiducial": 40}}}))
    dd = root / "FHC" / "1A" / "files" / "data7"; dd.mkdir(parents=True)          # a data file sharing the directory is ignored
    np.savez_compressed(dd / "reco_data7.npz", __meta__=json.dumps({"cache_version": 3, "is_mc": False}), passed=np.ones(4, bool), reco_p=np.arange(4.0))
    (dd / "manifest_data7.json").write_text(json.dumps({"tag": "data7", "kind": "data", "status": "ok", "pot_used": 7.0e17, "n_reco": 4}))
    mf = merge_playlist(root, "FHC", "1A", "mc", log=lambda *a: None)
    md = merge_playlist(root, "FHC", "1A", "data", log=lambda *a: None)
    assert md["tags"] == ["data7"] and abs(json.loads((root / "FHC/1A/pot_1A_data.json").read_text())["pot_used"] - 7.0e17) < 1
    pl = json.loads((root / "FHC/1A/pot_1A_mc.json").read_text())
    assert abs(pl["pot_used"] - 4.0e18) < 1 and pl["n_files"] == 2 and mf["selection_cutflow_sum"]["IsoBlobs"] == 6
    r = load_reco_npz(root / "FHC/1A/reco_1A_mc.npz")
    assert len(r["reco_p"]) == 400 and abs(r["__meta__"]["pot"] - 4.0e18) < 1
    t = TruthTable.load(root / "FHC/1A/truth_1A_skim.npz")
    assert abs(t.norm.pot - 4.0e18) < 1 and not t.has_fs and "lp_p" in t
    assert mf["signal_cutflow_sum"]["numu_cc_muon"] == {"n": 200, "n_fiducial": 80}
    # chunked access: one playlist at a time, same tables and POT as the merged products
    from types import SimpleNamespace
    from ndp.products import iter_mc_chunks, iter_reco_chunks
    fake = SimpleNamespace(name="fake", data={"beam": "FHC", "playlists": {"mc": ["1A"], "data": ["1A"]}, "products_dir": str(root)})
    chunks = list(iter_reco_chunks(None, fake, "mc"))
    assert [c[0] for c in chunks] == ["FHC/1A"] and abs(chunks[0][2] - 4.0e18) < 1 and len(chunks[0][1]["reco_p"]) == 400
    label, rm, rt, tt, pot, srcs = next(iter_mc_chunks(None, fake))
    assert label == "FHC/1A" and rt.n == 400 and tt.n == t.n and abs(pot - 4.0e18) < 1 and len(srcs) == 3
    dchunks = list(iter_reco_chunks(None, fake, "data"))
    assert len(dchunks) == 1 and len(dchunks[0][1]["reco_p"]) == 4 and abs(dchunks[0][2] - 7.0e17) < 1


def test_products_loader_legacy_matches_direct_caches():
    if not (have_mc_cache() and have_data_file()):
        skip("needs the cached data + MC tables")
    from ndp.products import load_reco, load_truth, load_reco_truth
    from ndp.adapters.minerva_anatuple import load_reco_cache
    cfg = site(); ch = legacy_channel("minerva_me_ccqelike_1mu1p")
    rd, pot_d, _ = load_reco(cfg, ch, "data")
    rm, pot_m, _ = load_reco(cfg, ch, "mc")
    assert abs(pot_d - 2.0497721920490272e17) / pot_d < 1e-12 and abs(pot_m - 9.988796749584837e18) / pot_m < 1e-12
    direct = load_reco_cache(cfg.data_dir / "cache/reco_mc110040.npz")
    assert np.array_equal(rm["reco_p"], direct["reco_p"]) and int(rd["passed"].sum()) == 844
    assert load_truth(cfg, ch).n == 544600 and load_reco_truth(cfg, ch).has_fs


def test_streamed_build_matches_local_cache_on_first_entries():
    """Stream the first entries of the local MC file's published twin and compare with a local build."""
    if not have_mc_cache():
        skip("needs the local MC file")
    try:
        import XRootD  # noqa: F401
    except ImportError:
        skip("no XRootD bindings")
    from ndp.adapters.minerva_anatuple import build_cache, load_reco_cache
    cfg = site(); ch = load_channel("minerva_me_ccqelike_1mu1p")
    local = cfg.data_dir / "MasterAnaDev_mc_AnaTuple_run00110040_Playlist.root"
    a, b = Path(tempfile.mkdtemp()), Path(tempfile.mkdtemp())
    try:
        ra = build_cache(str(local), a, is_mc=True, entry_stop=2000, channel=ch, skim=True, sidecar=True, log=lambda *x: None)
        rb = build_cache(URL, b, is_mc=True, entry_stop=2000, channel=ch, skim=True, sidecar=True, log=lambda *x: None)
    except Exception as e:  # network not available on this node
        skip(f"streaming unavailable: {type(e).__name__}: {e}")
    assert ra["fingerprint"]["sha256_partial"] == rb["fingerprint"]["sha256_partial"] and ra["pot"] == rb["pot"]
    for name in ("reco_mc110040.npz",):
        x, y = load_reco_cache(a / name), load_reco_cache(b / name)
        for k in x:
            if k != "__meta__":
                assert np.array_equal(x[k], y[k], equal_nan=True), k
    for name in ("truth_mc110040.npz", "truth_mc110040_skim.npz", "reco_mc110040_truthcols_skim.npz"):
        x, y = TruthTable.load(a / name), TruthTable.load(b / name)
        assert x.n == y.n and all(np.array_equal(x[k], y[k], equal_nan=True) for k in x.columns)
    sa, sb = json.loads((a / "manifest_mc110040.json").read_text()), json.loads((b / "manifest_mc110040.json").read_text())
    assert sa["selection_cutflow"] == sb["selection_cutflow"] and sa["signal_cutflow"] == sb["signal_cutflow"]
    assert sa["status"] == "ok" and "n_truth_skim" in sa and sa["total_truth_entries"] == 544600
    assert TruthTable.load(b / "reco_mc110040_truthcols_skim.npz").n == len(load_reco_cache(b / "reco_mc110040.npz")["passed"])
