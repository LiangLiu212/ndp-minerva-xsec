"""Standardized schema for analysis cross-check diagnostic payloads.

Experiment-agnostic: any neutrino cross-section analysis (MINERvA, DUNE,
MicroBooNE, T2K, ...) emits a JSON file conforming to this schema. The
diagnostic tool consumes the validated payload without knowing anything
about ROOT files, branch names, or detector geometry.

Usage:
    from scripts.diagnostic_schema import load_payload
    payload = load_payload("path/to/payload.json")
"""

import json
import math
from dataclasses import dataclass, field
from pathlib import Path


SCHEMA_VERSION = "1.3"
# Versions whose payloads this module can load. 1.1 adds the optional
# `cut_variables` section; 1.2 adds the optional per-stage `eff_numerator`
# (truth-binned efficiency numerator); 1.3 adds the optional top-level
# `cutflow_crosscheck` section (Python-vs-external integrated stage-total
# comparison). Older payloads remain valid (the newer sections absent).
SUPPORTED_VERSIONS = ("1.0", "1.1", "1.2", "1.3")


@dataclass
class Metadata:
    experiment: str
    channel: str
    variable: str
    variable_label: str
    units: str
    signal_definition: str
    pot_data: float
    pot_mc: float


@dataclass
class Binning:
    edges: list[float]

    @property
    def n_bins(self) -> int:
        return len(self.edges) - 1

    @property
    def centres(self):
        return [(self.edges[i] + self.edges[i + 1]) / 2
                for i in range(self.n_bins)]

    @property
    def widths(self):
        return [self.edges[i + 1] - self.edges[i]
                for i in range(self.n_bins)]


@dataclass
class Truth:
    counts: list[float]
    n_total: int


@dataclass
class Stage:
    name: str
    cumulative_cuts: list[str]
    mc_signal: list[float]
    mc_background: list[float]
    data: list[float]
    n_mc_signal: int
    n_mc_background: int
    n_data: int
    # Optional (schema v1.2+). Truth-binned efficiency numerator: MC events
    # that are truth-signal-in-phase-space AND pass the reco selection through
    # this stage, binned by the TRUE kinematic variable. This is the correct
    # efficiency numerator (paired with truth.counts as denominator) — distinct
    # from `mc_signal`, which is reco-binned reco-selected signal for purity /
    # data-MC. Absent in 1.0/1.1 payloads (efficiency falls back to mc_signal).
    eff_numerator: list[float] | None = None


@dataclass
class CutVariable:
    """Distribution of a single cut variable, data vs MC.

    Experiment-agnostic: `name`/`label`/`units`/`edges` describe an
    arbitrary observable the analysis cuts on (vertex Z, apothem radius,
    a BDT score, a boolean pass/fail with 2-bin edges, ...). Histograms
    are filled at the selection stage named by `stage`. Observable-
    independent, so a given analysis emits the same `cut_variables`
    block in every per-observable payload.
    """
    name: str
    label: str
    units: str
    stage: str
    edges: list[float]
    mc_signal: list[float]
    mc_background: list[float]
    data: list[float]

    @property
    def n_bins(self) -> int:
        return len(self.edges) - 1

    @property
    def centres(self):
        return [(self.edges[i] + self.edges[i + 1]) / 2
                for i in range(self.n_bins)]

    @property
    def widths(self):
        return [self.edges[i + 1] - self.edges[i]
                for i in range(self.n_bins)]


@dataclass
class CutflowCrosscheck:
    """Cross-check of integrated per-stage totals vs an external source.

    Optional (schema v1.3+). Records a comparison of the Python companion
    pass's own scalar stage totals against the C++ `runEventLoop` `mycuts`
    cut-summary table (dumped to `cutflow_summary_*.json`). A `status` of
    "fail" means the two selection implementations diverge — a correctness
    bug — and `diagnostic_plots.py` escalates it to a hard sanity violation.
    """
    source: str
    status: str
    comparisons: list[dict] = field(default_factory=list)
    rel_tol: float | None = None


@dataclass
class DiagnosticPayload:
    schema_version: str
    metadata: Metadata
    binning: Binning
    truth: Truth
    stages: list[Stage] = field(default_factory=list)
    cut_variables: list[CutVariable] = field(default_factory=list)
    cutflow_crosscheck: CutflowCrosscheck | None = None


class PayloadError(ValueError):
    pass


def _check(condition: bool, msg: str):
    if not condition:
        raise PayloadError(msg)


def _check_hist(vals, label: str):
    """Every bin of a histogram array must be a finite, non-negative number.

    Rejects non-numeric dtypes (str, bool), NaN/inf, and negatives — the
    three per-bin failure modes a malformed producer can emit.
    """
    for i, v in enumerate(vals):
        _check(isinstance(v, (int, float)) and not isinstance(v, bool),
               f"{label}[{i}] must be a number, got {type(v).__name__}")
        _check(math.isfinite(v),
               f"{label}[{i}] must be finite, got {v}")
        _check(v >= 0,
               f"{label}[{i}] must be >= 0, got {v}")


def validate_payload(payload: DiagnosticPayload):
    _check(payload.schema_version in SUPPORTED_VERSIONS,
           f"schema_version must be one of {SUPPORTED_VERSIONS}, "
           f"got '{payload.schema_version}'")

    n = payload.binning.n_bins
    _check(n >= 1, f"binning must have >= 2 edges, got {len(payload.binning.edges)}")

    edges = payload.binning.edges
    for i in range(n):
        _check(edges[i] < edges[i + 1],
               f"bin edges must be strictly increasing: "
               f"edges[{i}]={edges[i]} >= edges[{i+1}]={edges[i+1]}")

    _check(len(payload.truth.counts) == n,
           f"truth.counts length {len(payload.truth.counts)} != {n} bins")
    _check(payload.truth.n_total >= 0,
           f"truth.n_total must be >= 0, got {payload.truth.n_total}")

    _check(len(payload.stages) >= 1, "at least one stage required")

    stage_names = [s.name for s in payload.stages]
    _check(len(set(stage_names)) == len(stage_names),
           f"stage names must be unique, got {stage_names}")
    _check("no_cuts" in stage_names,
           "a stage named 'no_cuts' is required (the pre-cut baseline); "
           f"got stages {stage_names}")

    for si, stage in enumerate(payload.stages):
        prefix = f"stages[{si}] ({stage.name!r})"
        _check(len(stage.mc_signal) == n,
               f"{prefix}: mc_signal length {len(stage.mc_signal)} != {n}")
        _check(len(stage.mc_background) == n,
               f"{prefix}: mc_background length {len(stage.mc_background)} != {n}")
        _check(len(stage.data) == n,
               f"{prefix}: data length {len(stage.data)} != {n}")
        _check_hist(stage.mc_signal, f"{prefix}: mc_signal")
        _check_hist(stage.mc_background, f"{prefix}: mc_background")
        _check_hist(stage.data, f"{prefix}: data")
        _check(stage.n_mc_signal >= 0,
               f"{prefix}: n_mc_signal must be >= 0")
        _check(stage.n_mc_background >= 0,
               f"{prefix}: n_mc_background must be >= 0")
        _check(stage.n_data >= 0,
               f"{prefix}: n_data must be >= 0")
        if stage.eff_numerator is not None:
            _check(len(stage.eff_numerator) == n,
                   f"{prefix}: eff_numerator length "
                   f"{len(stage.eff_numerator)} != {n}")
            _check(all(v >= 0 for v in stage.eff_numerator),
                   f"{prefix}: eff_numerator has negative entries")

    for ci, cv in enumerate(payload.cut_variables):
        prefix = f"cut_variables[{ci}] ({cv.name!r})"
        cn = cv.n_bins
        _check(cn >= 1,
               f"{prefix}: edges must have >= 2 values, got {len(cv.edges)}")
        for i in range(cn):
            _check(cv.edges[i] < cv.edges[i + 1],
                   f"{prefix}: edges must be strictly increasing: "
                   f"edges[{i}]={cv.edges[i]} >= edges[{i+1}]={cv.edges[i+1]}")
        for fld in ("mc_signal", "mc_background", "data"):
            vals = getattr(cv, fld)
            _check(len(vals) == cn,
                   f"{prefix}: {fld} length {len(vals)} != {cn} bins")
            _check(all(v >= 0 for v in vals),
                   f"{prefix}: {fld} has negative entries")

    cc = payload.cutflow_crosscheck
    if cc is not None:
        _check(cc.status in ("pass", "fail"),
               f"cutflow_crosscheck.status must be 'pass'/'fail', got {cc.status!r}")
        for ci, comp in enumerate(cc.comparisons):
            _check(isinstance(comp, dict),
                   f"cutflow_crosscheck.comparisons[{ci}] must be an object")
            for k in ("stage", "quantity", "match"):
                _check(k in comp,
                       f"cutflow_crosscheck.comparisons[{ci}] missing '{k}'")
        # status must be consistent with the comparisons it reports.
        any_fail = any(not c.get("match", True) for c in cc.comparisons)
        _check((cc.status == "fail") == any_fail,
               f"cutflow_crosscheck.status '{cc.status}' inconsistent with "
               f"comparisons (any mismatch: {any_fail})")

    m = payload.metadata
    _check(m.pot_data > 0, f"pot_data must be > 0, got {m.pot_data}")
    _check(m.pot_mc > 0, f"pot_mc must be > 0, got {m.pot_mc}")


def _parse_metadata(d: dict) -> Metadata:
    required = ("experiment", "channel", "variable", "variable_label",
                "units", "signal_definition", "pot_data", "pot_mc")
    for k in required:
        _check(k in d, f"metadata missing required field '{k}'")
    return Metadata(**{k: d[k] for k in required})


def _parse_stage(d: dict) -> Stage:
    required = ("name", "cumulative_cuts", "mc_signal", "mc_background",
                "data", "n_mc_signal", "n_mc_background", "n_data")
    for k in required:
        _check(k in d, f"stage missing required field '{k}'")
    kwargs = {k: d[k] for k in required}
    # Optional (schema v1.2+).
    if "eff_numerator" in d:
        kwargs["eff_numerator"] = d["eff_numerator"]
    return Stage(**kwargs)


def _parse_cut_variable(d: dict) -> CutVariable:
    required = ("name", "label", "units", "stage", "edges",
                "mc_signal", "mc_background", "data")
    for k in required:
        _check(k in d, f"cut_variable missing required field '{k}'")
    return CutVariable(**{k: d[k] for k in required})


def _parse_cutflow_crosscheck(d: dict) -> CutflowCrosscheck:
    for k in ("source", "status"):
        _check(k in d, f"cutflow_crosscheck missing required field '{k}'")
    return CutflowCrosscheck(
        source=d["source"],
        status=d["status"],
        comparisons=d.get("comparisons", []),
        rel_tol=d.get("rel_tol"),
    )


def load_payload(path) -> DiagnosticPayload:
    path = Path(path)
    _check(path.exists(), f"payload file not found: {path}")

    with open(path) as f:
        raw = json.load(f)

    _check("schema_version" in raw, "missing 'schema_version'")
    _check("metadata" in raw, "missing 'metadata'")
    _check("binning" in raw, "missing 'binning'")
    _check("truth" in raw, "missing 'truth'")
    _check("stages" in raw, "missing 'stages'")

    payload = DiagnosticPayload(
        schema_version=raw["schema_version"],
        metadata=_parse_metadata(raw["metadata"]),
        binning=Binning(edges=raw["binning"]["edges"]),
        truth=Truth(counts=raw["truth"]["counts"],
                    n_total=raw["truth"]["n_total"]),
        stages=[_parse_stage(s) for s in raw["stages"]],
        cut_variables=[_parse_cut_variable(c)
                       for c in raw.get("cut_variables", [])],
        cutflow_crosscheck=(
            _parse_cutflow_crosscheck(raw["cutflow_crosscheck"])
            if "cutflow_crosscheck" in raw else None),
    )
    validate_payload(payload)
    return payload
