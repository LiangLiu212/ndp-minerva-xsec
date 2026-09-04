"""Total cross sections from GENIE spline XML files (streaming parse, cached).

GENIE stores one spline per (algorithm, initial state, process). The total CC cross
section for a nuclide is the sum of every spline whose name carries `nu:<pdg>;tgt:<pdg>`
and `proc:Weak[CC]`. Knot values are in GENIE natural units (1/GeV^2); converted here to
cm^2 with 1 GeV^-2 = 0.3893793e-27 cm^2.
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np

GEV2_TO_CM2 = 0.3893793e-27

_NAME = re.compile(r"nu:(-?\d+);tgt:(\d+);.*proc:([^;]*);")


def _cache_key(xml_path: Path, nu_pdg: int, targets: list[int], proc: str) -> str:
    st = xml_path.stat()
    raw = f"{xml_path}|{st.st_size}|{int(st.st_mtime)}|{nu_pdg}|{sorted(targets)}|{proc}"
    return hashlib.sha256(raw.encode()).hexdigest()[:16]


def total_xsec_splines(xml_path: str | Path, nu_pdg: int, targets: list[int], proc_contains: str = "Weak[CC]",
                       cache_dir: str | Path | None = None, e_grid: np.ndarray | None = None) -> dict:
    """{target_pdg: {"E": GeV grid, "sigma_cm2": total, "n_splines": int}} summed over processes.

    Splines are re-sampled by log-log interpolation onto a common grid before summing,
    so algorithms with different knot placements add correctly.
    """
    xml_path = Path(xml_path)
    targets = [int(t) for t in targets]
    if cache_dir is not None:
        cache = Path(cache_dir) / f"splines_{_cache_key(xml_path, nu_pdg, targets, proc_contains)}.npz"
        if cache.exists():
            z = np.load(cache, allow_pickle=False)
            return {int(t): {"E": z["E"], "sigma_cm2": z[f"sig_{t}"], "n_splines": int(z[f"n_{t}"])} for t in targets}
    e_grid = np.geomspace(0.01, 1000.0, 400) if e_grid is None else np.asarray(e_grid, float)
    sums = {t: np.zeros_like(e_grid) for t in targets}
    counts = {t: 0 for t in targets}
    want = set(targets)
    for _, el in ET.iterparse(str(xml_path), events=("end",)):
        if el.tag != "spline":
            continue
        name = el.get("name", "")
        m = _NAME.search(name)
        if m and int(m.group(1)) == nu_pdg and int(m.group(2)) in want and proc_contains in m.group(3):
            E = np.array([float(k.find("E").text) for k in el.iter("knot")])
            x = np.array([float(k.find("xsec").text) for k in el.iter("knot")])
            pos = x > 0
            if pos.sum() >= 2:
                logs = np.interp(np.log(e_grid), np.log(E[pos]), np.log(x[pos]), left=-np.inf, right=-np.inf)
                # do not extrapolate below the first positive knot (threshold) or above the last knot
                inside = (e_grid >= E[pos][0]) & (e_grid <= E[pos][-1])
                sums[int(m.group(2))] += np.where(inside, np.exp(logs), 0.0)
                counts[int(m.group(2))] += 1
        el.clear()
    out = {t: {"E": e_grid, "sigma_cm2": sums[t] * GEV2_TO_CM2, "n_splines": counts[t]} for t in targets}
    if cache_dir is not None:
        Path(cache_dir).mkdir(parents=True, exist_ok=True)
        np.savez_compressed(cache, E=e_grid, **{f"sig_{t}": out[t]["sigma_cm2"] for t in targets},
                            **{f"n_{t}": counts[t] for t in targets})
    return out


def per_nucleon_total(splines: dict, mass_fractions: dict) -> tuple[np.ndarray, np.ndarray]:
    """sigma_per_nucleon(E) for a target mix given by mass fractions w_i (sum 1).

    Nuclei of type i per unit mass ~ w_i / A_i, nucleons per unit mass ~ w_i, so
    sigma/nucleon = sum_i (w_i/A_i) sigma_i(E) / sum_i w_i.
    """
    E = None
    tot = None
    wsum = 0.0
    for pdg, w in mass_fractions.items():
        pdg = int(pdg)
        A = (pdg // 10) % 1000
        s = splines[pdg]
        if E is None:
            E, tot = s["E"], np.zeros_like(s["E"])
        tot = tot + (w / A) * s["sigma_cm2"]
        wsum += w
    return E, tot / wsum


def flux_averaged_per_nucleon(E: np.ndarray, sigma: np.ndarray, flux_edges: np.ndarray, flux_density: np.ndarray,
                              e_min: float = 0.0, e_max: float | None = None) -> float:
    from .flux import flux_average
    def sig(e):
        return np.interp(e, E, sigma, left=0.0, right=sigma[-1])
    return flux_average(flux_edges, flux_density, sig, e_min=e_min, e_max=e_max)


def check_spline_coverage(xml_path: str | Path, nu_pdg: int, targets: list[int], chunk: int = 1 << 24) -> list[int]:
    """Return the target PDGs for which the spline file holds no `nu:<pdg>;tgt:<t>;` entry.

    A streaming byte search (a few seconds on a 500 MB file), run before launching gevgen:
    GENIE aborts with `Assertion fUseSplines failed` when any nuclide of a target mix has no
    splines for the tune, and says nothing about which one.
    """
    needles = {int(t): f"nu:{nu_pdg};tgt:{int(t)};".encode() for t in targets}
    found = set()
    with open(xml_path, "rb") as fh:
        carry = b""
        while True:
            block = fh.read(chunk)
            if not block:
                break
            data = carry + block
            for t, n in needles.items():
                if t not in found and n in data:
                    found.add(t)
            if len(found) == len(needles):
                break
            carry = data[-64:]
    return [t for t in needles if t not in found]


def spline_tune_name(xml_path: str | Path) -> str | None:
    """The `<genie_tune name="...">` header of a spline file (first 4 kB), or None."""
    head = Path(xml_path).open("rb").read(4096).decode(errors="ignore")
    m = re.search(r'genie_tune name="([^"]+)"', head)
    return m.group(1) if m else None
