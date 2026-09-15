"""Playlist-level data products and the one place consumers load MINERvA tables from.

Two layouts are supported, chosen by the channel manifest's `data:` block:

1. **Playlist products** (grid campaigns, `ndp grid ...`): the channel lists
   `playlists: {mc: [1A, ...], data: [1A, ...]}` (+ optional `beam`, default FHC, and
   `products_dir`, default `<data_dir>/products`). Under `<products_dir>/<beam>/<playlist>/`
   a merge (`merge_playlist`) has written

       reco_<pl>_data.npz / reco_<pl>_mc.npz     the concatenated reco tables (cache v3 columns), per kind
       truth_<pl>_skim.npz                       TruthTable skim of the MC Truth trees (derived columns, no fs_*)
       reco_<pl>_truthcols_skim.npz              TruthTable of the MC reco rows' truth (derived columns)
       pot_<pl>_data.json / pot_<pl>_mc.json     {pot_used, pot_total, n_files, files: [{tag, pot_used, ...}]}
       merge_<pl>_data.json / merge_<pl>_mc.json fingerprints of the merged per-file products

   and the per-file products live in `<products_dir>/<beam>/<playlist>/files/<tag>/`.
   Several playlists are concatenated on load (POT summed).

2. **Per-file caches** (the original single-file layout, `ndp data cache`): `reco_data_files` /
   `reco_mc_files` under `<data_dir>` with `truth_<tag>.npz`, `reco_<tag>.npz`,
   `reco_<tag>_truthcols.npz` in `<data_dir>/cache/`. POT comes from the reco cache meta
   (caches built after 2026-09-14), the truth cache meta, a `manifest_<tag>.json` sidecar, or —
   last — the original AnaTuple's Meta tree when it is present locally.

Every loader returns the POT with the table so nobody reopens an AnaTuple to normalise.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from .events import TruthTable


# ---- manifest helpers -------------------------------------------------------------------------
def has_products(channel) -> bool:
    return bool((channel.data or {}).get("playlists"))


def products_root(cfg, channel) -> Path:
    d = (channel.data or {}).get("products_dir")
    return Path(d) if d else cfg.require("data_dir") / "products"


def playlists(channel, kind: str) -> list[tuple[str, str]]:
    """[(beam, playlist), ...] for kind 'mc' | 'data'."""
    pl = (channel.data or {}).get("playlists") or {}
    beam = str((channel.data or {}).get("beam", "FHC"))
    out = []
    for it in pl.get(kind) or []:
        if isinstance(it, dict):
            out.append((str(it.get("beam", beam)), str(it["playlist"])))
        else:
            out.append((beam, str(it)))
    return out


def playlist_dir(cfg, channel, beam: str, playlist: str) -> Path:
    return products_root(cfg, channel) / beam / playlist


def read_json(p: Path) -> dict:
    return json.loads(Path(p).read_text())


# ---- reco tables ------------------------------------------------------------------------------
def load_reco_npz(path: Path) -> dict:
    from .adapters.minerva_anatuple import load_reco_cache
    return load_reco_cache(path)


def concat_reco(tables: list[dict], sources: list[str], pots: list[float]) -> dict:
    """Concatenate reco tables column by column; refuse differing column sets / cache versions."""
    cols0 = sorted(k for k in tables[0] if k != "__meta__")
    v0 = int(tables[0].get("__meta__", {}).get("cache_version", 1))
    for src, t in zip(sources, tables):
        cols = sorted(k for k in t if k != "__meta__")
        v = int(t.get("__meta__", {}).get("cache_version", 1))
        if cols != cols0 or v != v0:
            raise ValueError(f"reco table {src} has columns/version ({len(cols)}, v{v}) != first ({len(cols0)}, v{v0}); "
                             f"differing columns: {sorted(set(cols) ^ set(cols0))}")
    out = {k: np.concatenate([t[k] for t in tables]) for k in cols0}
    m0 = dict(tables[0].get("__meta__", {}))
    out["__meta__"] = {"cache_version": v0, "columns": cols0, "is_mc": m0.get("is_mc"), "units": m0.get("units"),
                       "selection": m0.get("selection"), "n_entries": int(len(out[cols0[0]])),
                       "pot": float(math.fsum(pots)) if pots else None,
                       "sources": [{"source": s, "n_entries": int(len(t[cols0[0]])), "pot": p} for s, t, p in zip(sources, tables, pots)]}
    return out


# ---- legacy per-file caches -------------------------------------------------------------------
def legacy_files(channel, kind: str) -> list[str]:
    return list((channel.data or {}).get("reco_mc_files" if kind == "mc" else "reco_data_files", []))


def legacy_pot(cfg, tag: str, original: Path | None, reco_meta: dict | None = None, truth_meta: dict | None = None) -> float:
    """POT_Used of one per-file cache: reco meta, truth meta, sidecar, or the original file's Meta tree."""
    if reco_meta and isinstance(reco_meta.get("pot"), dict) and reco_meta["pot"].get("pot_used") is not None:
        return float(reco_meta["pot"]["pot_used"])
    if truth_meta and (truth_meta.get("norm") or {}).get("kind") == "pot" and truth_meta["norm"].get("pot"):
        return float(truth_meta["norm"]["pot"])
    cache = cfg.require("data_dir") / "cache"
    side = cache / f"manifest_{tag}.json"
    if side.exists():
        return float(read_json(side)["pot_used"])
    if original is not None and Path(original).exists():
        from .adapters.minerva_anatuple import read_pot
        return float(read_pot(original)["pot_used"])
    raise FileNotFoundError(f"no POT for cache tag {tag}: no reco/truth meta, no manifest_{tag}.json, no original file")


# ---- public loaders ---------------------------------------------------------------------------
def load_reco(cfg, channel, kind: str) -> tuple[dict, float, list[str]]:
    """(reco table, POT_Used, sources) for kind 'mc' | 'data', products or legacy caches."""
    from .adapters.minerva_anatuple import cache_tag
    if has_products(channel):
        tables, srcs, pots = [], [], []
        for beam, pl in playlists(channel, kind):
            d = playlist_dir(cfg, channel, beam, pl)
            p = d / f"reco_{pl}_{kind}.npz"
            if not p.exists():
                raise FileNotFoundError(f"missing playlist product {p} (run `ndp data merge --beam {beam} --playlist {pl} --kind {kind}`)")
            t = load_reco_npz(p)
            pj = d / f"pot_{pl}_{kind}.json"
            pot = read_json(pj)["pot_used"] if pj.exists() else t["__meta__"].get("pot")
            tables.append(t); srcs.append(str(p)); pots.append(float(pot))
        if not tables:
            raise FileNotFoundError(f"channel {channel.name} lists no {kind} playlists")
        r = tables[0] if len(tables) == 1 else concat_reco(tables, srcs, pots)
        if len(tables) == 1:
            r["__meta__"] = dict(r.get("__meta__", {}), pot=pots[0])
        return r, float(math.fsum(pots)), srcs
    data_dir = cfg.require("data_dir"); cache = data_dir / "cache"
    tables, srcs, pots = [], [], []
    for fn in legacy_files(channel, kind):
        tag = cache_tag(fn)
        p = cache / f"reco_{tag}.npz"
        if not p.exists():
            raise FileNotFoundError(f"missing reco cache {p}: `python -m ndp data cache --channel {channel.name}`")
        t = load_reco_npz(p)
        tm = cache / f"truth_{tag}.npz"
        truth_meta = TruthTable.load(tm).meta if (kind == "mc" and tm.exists()) else None
        tables.append(t); srcs.append(str(p)); pots.append(legacy_pot(cfg, tag, data_dir / fn, t.get("__meta__"), truth_meta))
    if not tables:
        raise FileNotFoundError(f"channel {channel.name} lists no {kind} files")
    r = tables[0] if len(tables) == 1 else concat_reco(tables, srcs, pots)
    return r, float(math.fsum(pots)), srcs


def load_truth(cfg, channel) -> TruthTable:
    """The MC truth table(s): playlist skims (derived columns, POT summed) or the legacy truth cache(s)."""
    from .adapters.minerva_anatuple import cache_tag
    if has_products(channel):
        tabs = []
        for beam, pl in playlists(channel, "mc"):
            p = playlist_dir(cfg, channel, beam, pl) / f"truth_{pl}_skim.npz"
            if not p.exists():
                raise FileNotFoundError(f"missing playlist product {p} (run `ndp data merge --beam {beam} --playlist {pl} --kind mc`)")
            tabs.append(TruthTable.load(p))
        if not tabs:
            raise FileNotFoundError(f"channel {channel.name} lists no mc playlists")
        return tabs[0] if len(tabs) == 1 else TruthTable.concatenate(tabs)
    data_dir = cfg.require("data_dir"); cache = data_dir / "cache"
    tabs = []
    for fn in legacy_files(channel, "mc"):
        p = cache / f"truth_{cache_tag(fn)}.npz"
        if not p.exists():
            raise FileNotFoundError(f"missing truth cache {p}: `python -m ndp data cache --channel {channel.name}`")
        tabs.append(TruthTable.load(p))
    if not tabs:
        raise FileNotFoundError(f"channel {channel.name} lists no reco_mc_files")
    return tabs[0] if len(tabs) == 1 else TruthTable.concatenate(tabs)


def load_reco_truth(cfg, channel) -> TruthTable:
    """The reco rows' truth (MC): playlist skims or the legacy `reco_<tag>_truthcols.npz`."""
    from .adapters.minerva_anatuple import cache_tag
    if has_products(channel):
        tabs = []
        for beam, pl in playlists(channel, "mc"):
            p = playlist_dir(cfg, channel, beam, pl) / f"reco_{pl}_truthcols_skim.npz"
            if not p.exists():
                raise FileNotFoundError(f"missing playlist product {p}")
            tabs.append(TruthTable.load(p))
        return tabs[0] if len(tabs) == 1 else TruthTable.concatenate(tabs)
    cache = cfg.require("data_dir") / "cache"
    tabs = [TruthTable.load(cache / f"reco_{cache_tag(fn)}_truthcols.npz") for fn in legacy_files(channel, "mc")]
    if not tabs:
        raise FileNotFoundError(f"channel {channel.name} lists no reco_mc_files")
    return tabs[0] if len(tabs) == 1 else TruthTable.concatenate(tabs)


def total_pot(cfg, channel, kind: str) -> float:
    """POT_Used summed over the channel's inputs of `kind`, from the merged `pot_<pl>_<kind>.json` files
    (or the legacy caches' meta) without loading any table."""
    from .adapters.minerva_anatuple import cache_tag
    if has_products(channel):
        tot = 0.0
        for beam, pl in playlists(channel, kind):
            pj = playlist_dir(cfg, channel, beam, pl) / f"pot_{pl}_{kind}.json"
            if not pj.exists():
                raise FileNotFoundError(f"missing {pj} (run `ndp data merge --beam {beam} --playlist {pl} --kind {kind}`)")
            tot += float(read_json(pj)["pot_used"])
        return tot
    data_dir = cfg.require("data_dir"); cache = data_dir / "cache"
    tot = 0.0
    for fn in legacy_files(channel, kind):
        tag = cache_tag(fn)
        t = load_reco_npz(cache / f"reco_{tag}.npz")
        tm = cache / f"truth_{tag}.npz"
        truth_meta = TruthTable.load(tm).meta if (kind == "mc" and tm.exists()) else None
        tot += float(legacy_pot(cfg, tag, data_dir / fn, t.get("__meta__"), truth_meta))
    return tot


# ---- chunked access (one playlist / legacy file at a time) ------------------------------------
def iter_reco_chunks(cfg, channel, kind: str):
    """Yield (label, reco table, POT_Used, source) one playlist product (or legacy cache) at a time.

    Consumers that only accumulate counts and histograms use this instead of `load_reco`, so the peak
    memory is one playlist, not the whole campaign (the 12 ME FHC playlists are ~90 GB of arrays).
    """
    from .adapters.minerva_anatuple import cache_tag
    if has_products(channel):
        n = 0
        for beam, pl in playlists(channel, kind):
            d = playlist_dir(cfg, channel, beam, pl)
            p = d / f"reco_{pl}_{kind}.npz"
            if not p.exists():
                raise FileNotFoundError(f"missing playlist product {p} (run `ndp data merge --beam {beam} --playlist {pl} --kind {kind}`)")
            t = load_reco_npz(p)
            pj = d / f"pot_{pl}_{kind}.json"
            pot = float(read_json(pj)["pot_used"] if pj.exists() else t["__meta__"].get("pot"))
            t["__meta__"] = dict(t.get("__meta__", {}), pot=pot)
            n += 1
            yield f"{beam}/{pl}", t, pot, str(p)
        if n == 0:
            raise FileNotFoundError(f"channel {channel.name} lists no {kind} playlists")
        return
    data_dir = cfg.require("data_dir"); cache = data_dir / "cache"
    files = legacy_files(channel, kind)
    if not files:
        raise FileNotFoundError(f"channel {channel.name} lists no {kind} files")
    for fn in files:
        tag = cache_tag(fn)
        p = cache / f"reco_{tag}.npz"
        if not p.exists():
            raise FileNotFoundError(f"missing reco cache {p}: `python -m ndp data cache --channel {channel.name}`")
        t = load_reco_npz(p)
        tm = cache / f"truth_{tag}.npz"
        truth_meta = TruthTable.load(tm).meta if (kind == "mc" and tm.exists()) else None
        pot = float(legacy_pot(cfg, tag, data_dir / fn, t.get("__meta__"), truth_meta))
        t["__meta__"] = dict(t.get("__meta__", {}), pot=pot)
        yield tag, t, pot, str(p)


def iter_mc_chunks(cfg, channel):
    """Yield (label, reco MC table, reco-side truth, MC truth, POT_Used, sources) per playlist product
    (or legacy MC file). The reco table and its reco-side truth come from the same product and are
    aligned row by row; the truth table is the playlist skim (or the full legacy truth cache)."""
    from .adapters.minerva_anatuple import cache_tag
    if has_products(channel):
        for label, rm, pot, src in iter_reco_chunks(cfg, channel, "mc"):
            beam, pl = label.split("/", 1)
            d = playlist_dir(cfg, channel, beam, pl)
            prt, pt = d / f"reco_{pl}_truthcols_skim.npz", d / f"truth_{pl}_skim.npz"
            for q in (prt, pt):
                if not q.exists():
                    raise FileNotFoundError(f"missing playlist product {q} (run `ndp data merge --beam {beam} --playlist {pl} --kind mc`)")
            rt = TruthTable.load(prt)
            if rt.n != len(rm[next(k for k in rm if k != "__meta__")]):
                raise ValueError(f"{prt}: {rt.n} rows but the reco table has {len(rm[next(k for k in rm if k != '__meta__')])}")
            yield label, rm, rt, TruthTable.load(pt), pot, [src, str(prt), str(pt)]
        return
    cache = cfg.require("data_dir") / "cache"
    for tag, rm, pot, src in iter_reco_chunks(cfg, channel, "mc"):
        prt, pt = cache / f"reco_{tag}_truthcols.npz", cache / f"truth_{tag}.npz"
        for q in (prt, pt):
            if not q.exists():
                raise FileNotFoundError(f"missing cache {q}: `python -m ndp data cache --channel {channel.name}`")
        yield tag, rm, TruthTable.load(prt), TruthTable.load(pt), pot, [src, str(prt), str(pt)]


def sources_fingerprints(cfg, channel) -> list[dict]:
    """Fingerprints of the local product / cache files a run reads (for manifests)."""
    from .io import cheap_fingerprint
    from .adapters.minerva_anatuple import cache_tag
    out = []
    if has_products(channel):
        for kind in ("data", "mc"):
            for beam, pl in playlists(channel, kind):
                d = playlist_dir(cfg, channel, beam, pl)
                for name in ([f"reco_{pl}_{kind}.npz", f"pot_{pl}_{kind}.json"] + ([f"truth_{pl}_skim.npz", f"reco_{pl}_truthcols_skim.npz"] if kind == "mc" else [])):
                    if (d / name).exists():
                        out.append({"role": f"{kind}_product", **cheap_fingerprint(d / name)})
        return out
    cache = cfg.require("data_dir") / "cache"
    for kind in ("data", "mc"):
        for fn in legacy_files(channel, kind):
            tag = cache_tag(fn)
            for name in ([f"reco_{tag}.npz"] + ([f"truth_{tag}.npz", f"reco_{tag}_truthcols.npz"] if kind == "mc" else [])):
                if (cache / name).exists():
                    out.append({"role": f"{kind}_cache", **cheap_fingerprint(cache / name)})
    return out


# ---- merge (per-file products -> playlist products) ---------------------------------------------
def merge_playlist(products_dir: Path, beam: str, playlist: str, kind: str, log=print) -> dict:
    """Merge `<products_dir>/<beam>/<playlist>/files/<tag>/` into the playlist products.

    kind 'mc': truth skims, reco tables, reco-truth skims; kind 'data': reco tables. POT is the
    fsum of every file's `manifest_<tag>.json` pot_used; a file without a sidecar is refused.
    """
    from .io import cheap_fingerprint, timestamp
    if kind not in ("data", "mc"):
        raise ValueError(f"kind must be 'data' or 'mc', not {kind!r}")
    d = Path(products_dir) / beam / playlist
    fdir = d / "files"
    all_tags = sorted(p.name for p in fdir.iterdir() if p.is_dir()) if fdir.exists() else []
    tags, sides = [], []
    for tag in all_tags:                                  # data and MC files share files/<tag>/: pick by kind
        sp = fdir / tag / f"manifest_{tag}.json"
        if not sp.exists():
            raise FileNotFoundError(f"{sp} missing: harvest the sidecar before merging")
        s = read_json(sp)
        if s.get("kind", "mc" if tag.startswith("mc") else "data") != kind:
            continue
        if s.get("status") != "ok":
            raise ValueError(f"file {tag} has status {s.get('status')!r}; resubmit it before merging")
        tags.append(tag); sides.append(s)
    if not tags:
        raise FileNotFoundError(f"no per-file {kind} products under {fdir}")
    pots = [float(s["pot_used"]) for s in sides]
    written = []
    # reco
    recos = [load_reco_npz(fdir / tag / f"reco_{tag}.npz") for tag in tags]
    r = concat_reco(recos, [str(fdir / tag / f"reco_{tag}.npz") for tag in tags], pots) if len(recos) > 1 else recos[0]
    if len(recos) == 1:
        r["__meta__"] = dict(r.get("__meta__", {}), pot=pots[0], sources=[{"source": str(fdir / tags[0]), "pot": pots[0]}])
    meta = r.pop("__meta__")
    np.savez_compressed(d / f"reco_{playlist}_{kind}.npz", __meta__=json.dumps(meta, default=str), **r)
    written.append(str(d / f"reco_{playlist}_{kind}.npz"))
    if kind == "mc":
        for stem, out_name in ((f"truth_{{tag}}_skim", f"truth_{playlist}_skim"), (f"reco_{{tag}}_truthcols_skim", f"reco_{playlist}_truthcols_skim")):
            tabs = [TruthTable.load(fdir / tag / (stem.format(tag=tag) + ".npz")) for tag in tags]
            if stem.startswith("reco_"):                       # the reco-side truth must align with the reco rows, file by file
                for tag, t, rt in zip(tags, tabs, recos):
                    n_reco = int(len(next(v for k, v in rt.items() if k != "__meta__")))
                    if t.n != n_reco:
                        raise ValueError(f"{tag}: reco_{tag}_truthcols_skim has {t.n} rows but reco_{tag} has {n_reco}; rebuild that file")
            t = tabs[0] if len(tabs) == 1 else TruthTable.concatenate(tabs)
            if stem.startswith("truth_"):
                t.meta["norm"] = {"kind": "pot", "pot": float(math.fsum(pots)), "xsec_per_unit_weight": None,
                                  "notes": f"sum of {len(tags)} files' POT_Used (manifest sidecars)"}
            t.save(d / f"{out_name}.npz"); written.append(str(d / f"{out_name}.npz"))
    pot = {"beam": beam, "playlist": playlist, "kind": kind, "n_files": len(tags), "pot_used": float(math.fsum(pots)),
           "pot_total": float(math.fsum(float(s.get("pot_total", 0.0)) for s in sides)),
           "files": [{"tag": s["tag"], "source": s.get("source"), "pot_used": s["pot_used"], "pot_total": s.get("pot_total"),
                      "n_reco": s.get("n_reco"), "n_truth": s.get("n_truth")} for s in sides]}
    (d / f"pot_{playlist}_{kind}.json").write_text(json.dumps(pot, indent=2))
    mf = {"merged": timestamp(), "beam": beam, "playlist": playlist, "kind": kind, "tags": tags,
          "inputs": [cheap_fingerprint(fdir / tag / f"manifest_{tag}.json") for tag in tags],
          "outputs": [cheap_fingerprint(w) for w in written], "pot": pot["pot_used"],
          "selection_cutflow_sum": _sum_cutflows(sides, "selection_cutflow"),
          "signal_cutflow_sum": _sum_signal_cutflows(sides)}
    (d / f"merge_{playlist}_{kind}.json").write_text(json.dumps(mf, indent=2, default=str))
    log(f"merged {beam}/{playlist} {kind}: {len(tags)} files, POT_Used {pot['pot_used']:.4e} -> {d}")
    return mf


def _sum_cutflows(sides: list[dict], key: str) -> dict:
    out: dict = {}
    for s in sides:
        for k, v in (s.get(key) or {}).items():
            out[k] = out.get(k, 0) + int(v)
    return out


def _sum_signal_cutflows(sides: list[dict]) -> dict:
    out: dict = {}
    for s in sides:
        for k, v in (s.get("signal_cutflow") or {}).items():
            cur = out.setdefault(k, {"n": 0, "n_fiducial": 0})
            cur["n"] += int(v.get("n", 0)); cur["n_fiducial"] += int(v.get("n_fiducial", 0))
    return out
