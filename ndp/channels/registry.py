"""Channel manifests: the physics contract a model is tested against.

A channel YAML declares everything the pipeline needs that is *not* the theorist's
model: experiment, signal definition, true phase space, observables + binning,
the reconstruction-level selection, where the data lives, normalisation constants,
and the default detector surrogate. All physics choices in a channel file are the
analyst's — the code only reads them. Fields marked `status:` in the YAML record
whether a value is ratified or a platform default awaiting the physicist's sign-off.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..io import load_yaml_or_json
from ..events import TruthTable
from .binning import Binning
from . import observables as obs

_CHANNELS_DIR = Path(__file__).resolve().parents[2] / "channels"


@dataclass
class ChannelSpec:
    name: str
    experiment: str
    description: str
    signal: dict                 # {nu_pdg: [14], current: "CC"}
    phase_space: dict            # {frame, vertex: {z_min_mm, z_max_mm, apothem_mm}, theta_max_deg, pz_min_gev, p_min_gev, p_max_gev}
    binning: Binning
    selection: dict              # {name, version, source}
    data: dict                   # {paper_manifest, reco_data_file, reco_mc_file, pot_data, pot_mc, ...}
    normalization: dict          # {phi_per_pot_cm2, n_nucleons, flux_table, target_mix, ...}
    surrogate: dict = field(default_factory=dict)
    raw: dict = field(default_factory=dict)
    path: Path | None = None

    # ---- truth-level predicates ---------------------------------------------------
    def is_signal(self, t: TruthTable) -> np.ndarray:
        nu = np.isin(t["nu_pdg"], self.signal.get("nu_pdg", [14]))
        cur = {"CC": 1, "NC": 2}.get(str(self.signal.get("current", "CC")).upper(), 1)
        return nu & (t["current"] == cur)

    def in_phase_space(self, t: TruthTable) -> np.ndarray:
        """True-kinematics + true-vertex phase space (the efficiency denominator)."""
        ps = self.phase_space
        frame = ps.get("frame", "detector")
        ok = np.ones(t.n, dtype=bool)
        v = ps.get("vertex")
        if v and t.meta.get("has_geometry", "vtx_z" in t):
            if "vtx_z" not in t:
                raise ValueError("sample claims has_geometry but has no vtx_* columns")
            z = t["vtx_z"]
            ok &= (z >= v["z_min_mm"]) & (z <= v["z_max_mm"])
            if "apothem_mm" in v:
                ok &= hex_apothem_mask(t["vtx_x"], t["vtx_y"], v["apothem_mm"])
        # A sample without geometry (e.g. a GENIE point-target run) is *defined* on the fiducial
        # target: the vertex requirement selects the target mass (n_nucleons), not kinematics, so
        # it is not applied. Samples with geometry (the experiment's own truth) do get it.
        if "theta_max_deg" in ps:
            ok &= obs.lep_theta(t, frame) <= np.deg2rad(ps["theta_max_deg"])
        if "pz_min_gev" in ps:
            ok &= obs.lep_pz(t, frame) >= ps["pz_min_gev"]
        if "p_min_gev" in ps:
            ok &= obs.lep_p(t, frame) >= ps["p_min_gev"]
        if "p_max_gev" in ps:
            ok &= obs.lep_p(t, frame) <= ps["p_max_gev"]
        return ok

    def observables(self, t: TruthTable):
        frame = self.phase_space.get("frame", "detector")
        return (obs.evaluate(self.binning.x_name, t, frame=frame),
                obs.evaluate(self.binning.y_name, t, frame=frame))

    def truth_cells(self, t: TruthTable, weights: np.ndarray | None = None, require_signal=True):
        """Histogram signal-and-in-phase-space events in true cells -> (sumw, sumw2, n_out, mask)."""
        mask = self.in_phase_space(t)
        if require_signal:
            mask &= self.is_signal(t)
        x, y = self.observables(t)
        w = t["weight"] if weights is None else np.asarray(weights)
        sumw, sumw2, n_out = self.binning.histogram(x[mask], y[mask], w[mask])
        return sumw, sumw2, n_out, mask

    def to_dict(self) -> dict:
        return dict(self.raw)


def hex_apothem_mask(x: np.ndarray, y: np.ndarray, apothem: float) -> np.ndarray:
    """MINERvA flat-top hexagonal fiducial (same inequality as tools/cc_inclusive_selector)."""
    slope = -1.0 / np.sqrt(3.0)
    intercept = 2.0 * apothem / np.sqrt(3.0)
    ax, ay = np.abs(x), np.abs(y)
    return (ax < apothem) & (ay < slope * ax + intercept)


def load_channel(name_or_path: str | Path) -> ChannelSpec:
    p = Path(name_or_path)
    if not p.exists():
        cand = _CHANNELS_DIR / f"{name_or_path}.yaml"
        if not cand.exists():
            raise FileNotFoundError(f"no channel {name_or_path!r} (looked for {cand})")
        p = cand
    raw = load_yaml_or_json(p)
    # PyYAML only parses floats with a dot AND a signed exponent ("3.23e30" stays a string);
    # coerce every numeric-looking normalisation / phase-space scalar so no caller has to.
    for sect in ("normalization", "phase_space"):
        for k, v in list(raw.get(sect, {}).items()):
            if isinstance(v, str):
                try:
                    raw[sect][k] = float(v)
                except ValueError:
                    pass
        if isinstance(raw.get(sect, {}).get("vertex"), dict):
            raw[sect]["vertex"] = {k: float(v) for k, v in raw[sect]["vertex"].items()}
    b = raw["binning"]
    binning = Binning(b["x"]["observable"], b["y"]["observable"], tuple(b["x"]["edges"]),
                      tuple(b["y"]["edges"]), b.get("global_cell_formula", "ix*n_y + iy"))
    return ChannelSpec(name=raw["name"], experiment=raw["experiment"], description=raw.get("description", ""),
                       signal=raw["signal"], phase_space=raw["phase_space"], binning=binning,
                       selection=raw.get("selection", {}), data=raw.get("data", {}),
                       normalization=raw.get("normalization", {}), surrogate=raw.get("surrogate", {}),
                       raw=raw, path=p)


def list_channels() -> list[str]:
    return sorted(p.stem for p in _CHANNELS_DIR.glob("*.yaml"))
