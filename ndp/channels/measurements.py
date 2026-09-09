"""Measurements: the observable pair + binning a model is compared on, at truth AND reco level.

A channel fixes what is *not* the theorist's model and not the analyst's choice of
observable: signal definition, true phase space, reconstruction-level selection, data
files, normalisation constants. A measurement adds the analyst's choice: which two
observables (one function of the truth, one function of the reconstruction), which bins,
and — if the experiment published an unfolded result on exactly this grid — which
release. Every channel has a `published` measurement (its own `binning`, the paper's
grid); an analyst adds their own as `measurements/<channel>/<name>.yaml` without
touching the channel's decided values:

    name: muon_p_theta
    channel: minerva_me_cc_inclusive_ptpz
    x: {observable: lep_p,          reco: reco_p,          units: GeV/c, edges: [...]}
    y: {observable: lep_theta_deg,  reco: reco_theta_deg,  units: deg,   edges: [...]}
    # omit y for a one-dimensional measurement

`observable` is a truth observable name or expression (`ndp.channels.observables`), `reco`
a reco observable name or expression over the cached reco columns
(`ndp.channels.reco_observables`). A user measurement has no published unfolded result,
so it is compared in folded space only: the model is pushed through a surrogate learned on
this measurement's grid from the experiment's paired MC and compared with the selected data
counts. That is the forward-folding workflow; nothing is unfolded.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..io import load_yaml_or_json
from ..events import TruthTable
from .binning import Binning
from . import observables as obs
from . import reco_observables as robs

_MEAS_DIR = Path(__file__).resolve().parents[2] / "measurements"

#: truth observable -> the reco observable a channel's published grid is measured in, when the
#: channel YAML does not say (`binning.x.reco`).
DEFAULT_RECO = {"lep_pT": "reco_pT", "lep_pz": "reco_pz", "lep_p": "reco_p", "lep_theta": "reco_theta",
                "lep_theta_deg": "reco_theta_deg", "lep_E": "reco_E_mu", "unit": "unit"}


@dataclass(frozen=True)
class Axis:
    truth: str                 # observable name or expression on the TruthTable
    reco: str                  # observable name or expression on the reco table
    edges: tuple
    units: str = ""
    log: bool = False
    label: str | None = None

    def axis_label(self, level: str) -> str:
        base = self.label or (self.truth if level == "truth" else self.reco)
        if level == "reco" and self.label:
            base = f"reco {self.label}"
        elif level == "truth" and self.label:
            base = f"true {self.label}"
        return f"{base} [{self.units}]" if self.units else base


@dataclass
class Measurement:
    name: str
    channel: str
    x: Axis
    y: Axis
    binning: Binning
    description: str = ""
    author: str = ""
    status: str = "user"
    release: str | None = None          # paper manifest id when the experiment published this exact grid
    surrogate: dict = field(default_factory=dict)
    notes: list = field(default_factory=list)
    is_1d: bool = False
    raw: dict = field(default_factory=dict)
    path: Path | None = None

    # ---- evaluation ---------------------------------------------------------------------
    def truth_observables(self, channel, t: TruthTable):
        frame = channel.phase_space.get("frame", "detector")
        return obs.evaluate(self.x.truth, t, frame=frame), obs.evaluate(self.y.truth, t, frame=frame)

    def reco_observables(self, r: dict):
        return robs.evaluate(self.x.reco, r), robs.evaluate(self.y.reco, r)

    def truth_cells(self, channel, t: TruthTable, weights=None, require_signal=True):
        """Histogram signal-and-in-phase-space events in true cells -> (sumw, sumw2, n_out, mask)."""
        mask = channel.in_phase_space(t)
        if require_signal:
            mask &= channel.is_signal(t)
        x, y = self.truth_observables(channel, t)
        w = t["weight"] if weights is None else np.asarray(weights)
        sumw, sumw2, n_out = self.binning.histogram(x[mask], y[mask], w[mask])
        return sumw, sumw2, n_out, mask

    # ---- bookkeeping ----------------------------------------------------------------------
    @property
    def is_published(self) -> bool:
        return self.release is not None

    def slug(self) -> str:
        return "" if self.name == "published" else self.name

    def surrogate_root(self, cfg) -> Path:
        """Where this measurement's surrogates live: surrogates/<channel>[/<measurement>]."""
        root = cfg.surrogates / self.channel
        return root if self.name == "published" else root / self.name

    def to_dict(self) -> dict:
        return {"name": self.name, "channel": self.channel, "description": self.description, "author": self.author,
                "status": self.status, "release": self.release, "is_1d": self.is_1d,
                "x": {"observable": self.x.truth, "reco": self.x.reco, "edges": list(self.x.edges), "units": self.x.units, "log": self.x.log},
                "y": {"observable": self.y.truth, "reco": self.y.reco, "edges": list(self.y.edges), "units": self.y.units, "log": self.y.log},
                "global_cell_formula": self.binning.formula, "surrogate": dict(self.surrogate), "notes": list(self.notes),
                "path": str(self.path) if self.path else None}


def _axis(d: dict, default_reco: bool = True) -> Axis:
    truth = str(d["observable"])
    reco = d.get("reco")
    if reco is None:
        if truth in DEFAULT_RECO:
            reco = DEFAULT_RECO[truth]
        else:
            raise ValueError(f"axis {truth!r} needs a `reco` observable (name or expression over the reco cache)")
    edges = tuple(float(v) for v in d["edges"])
    return Axis(truth=truth, reco=str(reco), edges=edges, units=str(d.get("units", "")), log=bool(d.get("log", False)),
                label=d.get("label"))


_UNIT_AXIS = Axis(truth="unit", reco="unit", edges=(0.0, 1.0), units="", label="(1D)")


def _from_raw(raw: dict, channel_name: str, *, name: str, path: Path | None, release: str | None,
              surrogate: dict | None, status: str) -> Measurement:
    x = _axis(raw["x"])
    is_1d = raw.get("y") in (None, {}, "none", "None")
    y = _UNIT_AXIS if is_1d else _axis(raw["y"])
    formula = raw.get("global_cell_formula", "ix*n_y + iy")
    binning = Binning(x.truth, y.truth, x.edges, y.edges, formula)
    return Measurement(name=name, channel=channel_name, x=x, y=y, binning=binning,
                       description=str(raw.get("description", "")), author=str(raw.get("author", "")),
                       status=status, release=release, surrogate=dict(surrogate or {}),
                       notes=list(raw.get("notes", []) if isinstance(raw.get("notes"), list) else ([raw["notes"]] if raw.get("notes") else [])),
                       is_1d=is_1d, raw=raw, path=path)


def published_measurement(channel) -> Measurement:
    """The channel's own grid (the paper's), with the release attached for the unfolded comparison."""
    b = dict(channel.raw["binning"])
    raw = {"x": b["x"], "y": b["y"], "global_cell_formula": b.get("global_cell_formula", "ix*n_y + iy"),
           "description": f"published grid of {channel.name} (arXiv:{channel.data.get('paper_manifest')})"}
    m = _from_raw(raw, channel.name, name="published", path=channel.path, release=channel.data.get("paper_manifest"),
                  surrogate=channel.surrogate, status=str(b.get("status", "decided")))
    if m.binning != channel.binning:
        raise ValueError("published measurement must reproduce the channel binning")
    return m


def measurement_files(channel_name: str) -> dict[str, Path]:
    d = _MEAS_DIR / channel_name
    return {p.stem: p for p in sorted(d.glob("*.y*ml"))} if d.exists() else {}


def list_measurements(channel) -> list[str]:
    inline = list((channel.raw.get("measurements") or {}).keys())
    return ["published", *inline, *[k for k in measurement_files(channel.name) if k not in inline]]


def load_measurement(channel, name_or_path: str | Path | None = None) -> Measurement:
    if name_or_path in (None, "", "published"):
        return published_measurement(channel)
    p = Path(name_or_path)
    inline = channel.raw.get("measurements") or {}
    if p.exists() and p.is_file():
        raw = load_yaml_or_json(p); name = raw.get("name", p.stem); path = p
    elif str(name_or_path) in inline:
        raw = dict(inline[str(name_or_path)]); name = str(name_or_path); path = channel.path
    elif str(name_or_path) in measurement_files(channel.name):
        path = measurement_files(channel.name)[str(name_or_path)]
        raw = load_yaml_or_json(path); name = raw.get("name", path.stem)
    else:
        raise FileNotFoundError(f"no measurement {name_or_path!r} for channel {channel.name} "
                                f"(known: {list_measurements(channel)}; or give a YAML path)")
    if raw.get("channel") not in (None, channel.name):
        raise ValueError(f"measurement {name} belongs to channel {raw['channel']!r}, not {channel.name!r}")
    release = raw.get("release")
    return _from_raw(raw, channel.name, name=name, path=path, release=release, surrogate=raw.get("surrogate"),
                     status=str(raw.get("status", "user")))
