"""Neutrino flux tables: parsing, unit conversion, integrals, and GENIE-ready histograms."""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

_ROW = re.compile(r"^\[\s*([0-9.]+)\s*-\s*([0-9.]+)\)\s+([0-9.eE+-]+)\s*$")


def load_flux_table(path: str | Path, table: str = "Table I") -> dict:
    """Parse a MINERvA supplemental flux table ("[lo - hi)<tab>value" rows) into edges/values.

    Returns {"edges": GeV, "values": as printed, "units": header line, "n_rows": int}.
    The arXiv:2110.13372 Table I is nu_mu / m^2 / 1e5 POT / GeV.
    """
    lines = Path(path).read_text().splitlines()
    in_table, units, edges, vals = False, "", [], []
    for ln in lines:
        if ln.strip().startswith(table):
            in_table = True
            continue
        if in_table and ln.strip().startswith("Table") and not ln.strip().startswith(table):
            break
        if not in_table:
            continue
        if ln.startswith("#"):
            units = ln.strip("# \n")
            continue
        m = _ROW.match(ln.strip())
        if m:
            lo, hi, v = float(m.group(1)), float(m.group(2)), float(m.group(3))
            if edges and abs(edges[-1] - lo) > 1e-9:
                raise ValueError(f"non-contiguous flux bins at {lo}")
            if not edges:
                edges.append(lo)
            edges.append(hi)
            vals.append(v)
    if not vals:
        raise ValueError(f"no flux rows found in {path} under {table!r}")
    return {"edges": np.array(edges), "values": np.array(vals), "units": units, "n_rows": len(vals),
            "source": str(path)}


def per_cm2_per_pot_per_gev(values: np.ndarray, units: str) -> np.ndarray:
    """Convert a printed flux density to nu / cm^2 / POT / GeV."""
    u = units.replace(" ", "").lower()
    if "m^{2}" in u and "10^{5}p.o.t" in u:
        return values * 1e-4 / 1e5
    raise ValueError(f"unrecognised flux units {units!r}; add a conversion")


def integrated_flux(edges: np.ndarray, density: np.ndarray, e_min: float = 0.0, e_max: float = 100.0) -> float:
    """Integral of a per-GeV density over [e_min, e_max] (partial bins clipped)."""
    lo = np.clip(edges[:-1], e_min, e_max)
    hi = np.clip(edges[1:], e_min, e_max)
    return float(np.sum(density * np.maximum(hi - lo, 0.0)))


def flux_average(edges: np.ndarray, density: np.ndarray, sigma_of_E, e_min=0.0, e_max=None) -> float:
    """<sigma> = int phi(E) sigma(E) dE / int phi(E) dE with sigma evaluated at bin centres
    (finer sub-sampling inside each bin so wide bins do not bias the average)."""
    e_max = edges[-1] if e_max is None else e_max
    num = den = 0.0
    for lo, hi, d in zip(edges[:-1], edges[1:], density):
        lo, hi = max(lo, e_min), min(hi, e_max)
        if hi <= lo or d <= 0:
            continue
        sub = np.linspace(lo, hi, 9)
        mid = 0.5 * (sub[:-1] + sub[1:])
        num += float(np.sum(d * sigma_of_E(mid) * np.diff(sub)))
        den += d * (hi - lo)
    return num / den if den > 0 else float("nan")


def write_th1_root(path: str | Path, edges: np.ndarray, density: np.ndarray, name: str = "flux",
                   title: str = "nu_mu flux") -> Path:
    """Write a TH1D (uproot) that `gevgen -f file.root,name` can consume."""
    import uproot
    path = Path(path)
    with uproot.recreate(path) as f:
        f[name] = (np.asarray(density, float), np.asarray(edges, float))
    return path


def load_channel_flux(channel, repo_root: Path) -> dict:
    """Resolve a channel's flux_table into converted arrays + the integrated flux."""
    rel = channel.normalization.get("flux_table")
    if not rel:
        raise KeyError("channel has no normalization.flux_table")
    tab = load_flux_table(repo_root / rel)
    dens = per_cm2_per_pot_per_gev(tab["values"], tab["units"])
    return {"edges": tab["edges"], "density_cm2_pot_gev": dens, "units_in": tab["units"], "source": tab["source"],
            "phi_0_100_cm2_per_pot": integrated_flux(tab["edges"], dens, 0.0, 100.0)}
