"""The common truth-level event table every adapter produces and every stage consumes.

A `TruthTable` is a dict of equal-length numpy columns plus a metadata dict. Scalar
columns are per event; final-state particles are stored *jagged* as a flat array plus
`fs_offsets` (length n+1, CSR style), so the core needs nothing beyond numpy.

Required scalar columns (all float64 unless noted):
    nu_pdg (int)      incoming neutrino PDG code
    E_nu              neutrino energy                             [GeV]
    lep_pdg (int)     primary final-state lepton PDG code
    lep_px lep_py lep_pz lep_E   primary lepton 4-momentum         [GeV]
    current (int)     1 = CC, 2 = NC
    int_type (int)    NDP interaction code (see INT_TYPES)
    target_Z target_A (int)   struck nucleus
    Q2                true four-momentum transfer squared          [GeV^2]
    W                 true hadronic invariant mass                 [GeV]
    weight            per-event weight (1 for unweighted generators)
Optional:
    vtx_x vtx_y vtx_z  true vertex                                 [mm, detector frame]
    fs_offsets (int64, n+1) + fs_pdg (int) fs_E fs_px fs_py fs_pz   final-state particles [GeV]

`meta["norm"]` says how weights turn into a cross section — see `Normalization`.
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from pathlib import Path
import json

import numpy as np

#: NDP interaction-type codes. Generator codes are mapped onto these by the adapters.
INT_TYPES = {0: "other", 1: "QE", 2: "RES", 3: "DIS", 4: "COH", 5: "MEC"}
INT_CODE = {v: k for k, v in INT_TYPES.items()}

REQUIRED = ("nu_pdg", "E_nu", "lep_pdg", "lep_px", "lep_py", "lep_pz", "lep_E", "current",
            "int_type", "target_Z", "target_A", "Q2", "W", "weight")
FS_COLUMNS = ("fs_pdg", "fs_E", "fs_px", "fs_py", "fs_pz")
INT_COLUMNS = {"nu_pdg", "lep_pdg", "current", "int_type", "target_Z", "target_A", "fs_pdg"}

M_MU = 0.1056583755  # GeV
M_P = 0.93827208816
M_N = 0.93956542052


@dataclass
class Normalization:
    """How to turn a sum of weights in a bin into a cross section (or an event rate).

    kind = "xsec_per_nucleon":
        d(sigma)/dx [cm^2/nucleon] = xsec_per_unit_weight * sum(w in bin) / dx
        i.e. `xsec_per_unit_weight = sigma_flux_avg_per_nucleon / sum(all weights)`.
        Absolute — what a generator run with a known flux-averaged total gives.
    kind = "pot":
        the sample corresponds to `pot` protons on target hitting the modelled target
        mass. Absolute cross sections need flux (phi_per_pot, nu/cm^2/POT) and
        n_nucleons from the channel; event rates need only the POT ratio.
    kind = "shape":
        no absolute normalisation — only shapes are meaningful.
    """
    kind: str = "shape"
    xsec_per_unit_weight: float | None = None   # cm^2/nucleon per unit weight
    pot: float | None = None
    notes: str = ""

    def to_dict(self) -> dict:
        return asdict(self)

    @staticmethod
    def from_dict(d: dict | None) -> "Normalization":
        d = dict(d or {})
        return Normalization(**{k: d.get(k) for k in ("kind", "xsec_per_unit_weight", "pot", "notes")
                                if d.get(k) is not None} | ({} if "kind" in d else {}))


class TruthTable:
    def __init__(self, columns: dict, meta: dict | None = None):
        cols = {k: np.asarray(v) for k, v in columns.items()}
        n = None
        for k in REQUIRED:
            if k not in cols:
                raise KeyError(f"TruthTable missing required column {k!r}")
        for k, v in cols.items():
            if k in FS_COLUMNS or k == "fs_offsets":
                continue
            if n is None:
                n = len(v)
            elif len(v) != n:
                raise ValueError(f"column {k!r} has length {len(v)}, expected {n}")
        if "fs_offsets" in cols:
            off = cols["fs_offsets"]
            if len(off) != n + 1:
                raise ValueError("fs_offsets must have length n+1")
            for k in FS_COLUMNS:
                if k in cols and len(cols[k]) != int(off[-1]):
                    raise ValueError(f"{k} length {len(cols[k])} != fs_offsets[-1] {off[-1]}")
        for k in INT_COLUMNS & set(cols):
            cols[k] = cols[k].astype(np.int64)
        self.columns = cols
        self.n = int(n if n is not None else 0)
        self.meta = dict(meta or {})
        self.meta.setdefault("norm", Normalization().to_dict())

    # ---- access -----------------------------------------------------------------
    def __getitem__(self, key: str) -> np.ndarray:
        return self.columns[key]

    def __contains__(self, key: str) -> bool:
        return key in self.columns

    def __len__(self) -> int:
        return self.n

    @property
    def norm(self) -> Normalization:
        return Normalization.from_dict(self.meta.get("norm"))

    @property
    def has_fs(self) -> bool:
        return "fs_offsets" in self.columns and "fs_pdg" in self.columns

    def with_weights(self, w: np.ndarray, note: str = "") -> "TruthTable":
        w = np.asarray(w, dtype=float)
        if w.shape != (self.n,):
            raise ValueError("weights must be one per event")
        cols = dict(self.columns)
        cols["weight"] = w
        meta = json.loads(json.dumps(self.meta, default=str))
        if note:
            meta.setdefault("reweight_history", []).append(note)
        return TruthTable(cols, meta)

    def select(self, mask: np.ndarray) -> "TruthTable":
        mask = np.asarray(mask, dtype=bool)
        cols = {}
        for k, v in self.columns.items():
            if k in FS_COLUMNS or k == "fs_offsets":
                continue
            cols[k] = v[mask]
        if self.has_fs:
            off = self.columns["fs_offsets"]
            starts, stops = off[:-1][mask], off[1:][mask]
            counts = stops - starts
            idx = np.concatenate([np.arange(s, e) for s, e in zip(starts, stops)]) if counts.size else np.zeros(0, int)
            cols["fs_offsets"] = np.concatenate([[0], np.cumsum(counts)]).astype(np.int64)
            for k in FS_COLUMNS:
                if k in self.columns:
                    cols[k] = self.columns[k][idx]
        return TruthTable(cols, self.meta)

    # ---- final-state helpers --------------------------------------------------------
    def fs_event_index(self) -> np.ndarray:
        """Event index of every final-state particle (length fs_offsets[-1])."""
        off = self.columns["fs_offsets"]
        return np.repeat(np.arange(self.n), np.diff(off))

    def fs_sum(self, values: np.ndarray, particle_mask: np.ndarray | None = None) -> np.ndarray:
        """Per-event sum of a per-particle array (optionally masked)."""
        if particle_mask is not None:
            values = np.where(particle_mask, values, 0.0)
        out = np.zeros(self.n)
        np.add.at(out, self.fs_event_index(), values)
        return out

    # ---- I/O -----------------------------------------------------------------------
    def save(self, path: str | Path) -> Path:
        path = Path(path)
        np.savez_compressed(path, __meta__=json.dumps(self.meta, default=str), **self.columns)
        return path if path.suffix == ".npz" else path.with_suffix(path.suffix + ".npz")

    @staticmethod
    def load(path: str | Path) -> "TruthTable":
        z = np.load(path, allow_pickle=False)
        meta = json.loads(str(z["__meta__"]))
        cols = {k: z[k] for k in z.files if k != "__meta__"}
        return TruthTable(cols, meta)

    @staticmethod
    def concatenate(tables: list["TruthTable"]) -> "TruthTable":
        if not tables:
            raise ValueError("nothing to concatenate")
        keys = set(tables[0].columns)
        for t in tables[1:]:
            if set(t.columns) != keys:
                raise ValueError("tables have different columns")
        cols = {}
        for k in keys:
            if k == "fs_offsets":
                parts, base = [tables[0].columns[k]], int(tables[0].columns[k][-1])
                for t in tables[1:]:
                    parts.append(t.columns[k][1:] + base)
                    base += int(t.columns[k][-1])
                cols[k] = np.concatenate(parts)
            else:
                cols[k] = np.concatenate([t.columns[k] for t in tables])
        meta = dict(tables[0].meta)
        meta["concatenated_from"] = [t.meta.get("source", "?") for t in tables]
        return TruthTable(cols, meta)

    def summary(self) -> dict:
        w = self.columns["weight"]
        it = self.columns["int_type"]
        return {
            "n_events": self.n, "sum_weights": float(w.sum()),
            "E_nu_median_gev": float(np.median(self.columns["E_nu"])) if self.n else None,
            "int_type_fractions": {INT_TYPES.get(int(k), str(k)): float(w[it == k].sum() / w.sum())
                                   for k in np.unique(it)} if self.n else {},
            "norm": self.meta.get("norm"),
        }
