"""Reusable per-stage histogram accumulator for DiagnosticPayload emission.

Experiment-agnostic: any neutrino cross-section analysis creates a
StageAccumulator, calls fill_* during its event loop, and gets a
validated DiagnosticPayload dict ready for json.dumps.

Usage:
    from scripts.diagnostic_accumulator import StageAccumulator

    acc = StageAccumulator(
        bin_edges=[0.0, 0.5, 1.0, 1.5, 2.0],
        stage_names=["no_cuts", "after_FV", "after_all_cuts"],
        cut_labels=["FV", "AllCuts"],
    )

    # MC reco loop — event passes first 2 cuts (stages 0, 1, 2):
    acc.fill_mc_through(max_stage=2, bin_index=3, is_signal=True)

    # Data loop:
    acc.fill_data_through(max_stage=2, bin_index=3)

    # Truth loop (no stages):
    acc.fill_truth(bin_index=1)

    # After loops:
    acc.set_truth_total(n_total=100000)
    payload = acc.to_payload(
        experiment="MINERvA", channel="CC-inclusive",
        variable="pT", variable_label="Muon p_T", units="GeV/c",
        signal_definition="numu CC", pot_data=2.47e18, pot_mc=4.97e19,
    )
"""

import numpy as np

SCHEMA_VERSION = "1.3"


class StageAccumulator:
    """Accumulates per-stage binned histograms for diagnostic payloads.

    Parameters
    ----------
    bin_edges : list[float]
        Bin edges (N+1 values for N bins), strictly increasing.
    stage_names : list[str]
        Ordered stage names. Stage 0 is typically "no_cuts".
    cut_labels : list[str] or None
        If provided, must have len(stage_names) - 1 entries.
        Stage i (for i >= 1) gets cumulative_cuts = cut_labels[:i].
        Stage 0 always gets cumulative_cuts = [].
    """

    def __init__(self, bin_edges, stage_names, cut_labels=None):
        self._edges = list(bin_edges)
        self._n_bins = len(self._edges) - 1
        if self._n_bins < 1:
            raise ValueError("bin_edges must have >= 2 values")
        self._stage_names = list(stage_names)
        self._n_stages = len(self._stage_names)
        if self._n_stages < 1:
            raise ValueError("at least one stage required")

        if cut_labels is not None:
            if len(cut_labels) != self._n_stages - 1:
                raise ValueError(
                    f"cut_labels length {len(cut_labels)} != "
                    f"n_stages - 1 ({self._n_stages - 1})"
                )
            self._cut_labels = list(cut_labels)
        else:
            self._cut_labels = None

        self._mc_signal = np.zeros((self._n_stages, self._n_bins))
        self._mc_background = np.zeros((self._n_stages, self._n_bins))
        self._data = np.zeros((self._n_stages, self._n_bins))
        # Truth-binned efficiency numerator (schema v1.3): truth-signal-in-
        # phase-space MC events that pass the reco selection through each stage,
        # binned by the TRUE kinematic variable. Paired with self._truth as the
        # efficiency denominator. Distinct from self._mc_signal (reco-binned).
        self._eff_numerator = np.zeros((self._n_stages, self._n_bins))
        self._truth = np.zeros(self._n_bins)

        self._n_mc_signal = np.zeros(self._n_stages, dtype=int)
        self._n_mc_background = np.zeros(self._n_stages, dtype=int)
        self._n_data = np.zeros(self._n_stages, dtype=int)
        self._truth_total = 0

        # Cut-variable distributions (observable-independent). Keyed by
        # name, preserving registration order for stable payload output.
        self._cut_vars = {}

        # Optional cross-check of the integrated per-stage totals against an
        # external source (e.g. the C++ runEventLoop `mycuts` cut-summary
        # table dumped to cutflow_summary_*.json). None until attached.
        self._cutflow_crosscheck = None

    def fill_mc(self, stage_index, bin_index, is_signal):
        """Fill a single MC event into one stage."""
        if bin_index < 0 or bin_index >= self._n_bins:
            return
        if is_signal:
            self._mc_signal[stage_index, bin_index] += 1
        else:
            self._mc_background[stage_index, bin_index] += 1

    def fill_mc_through(self, max_stage, bin_index, is_signal):
        """Fill a single MC event into stages 0..max_stage (inclusive)."""
        if bin_index < 0 or bin_index >= self._n_bins:
            return
        end = min(max_stage + 1, self._n_stages)
        if is_signal:
            self._mc_signal[:end, bin_index] += 1
        else:
            self._mc_background[:end, bin_index] += 1

    def count_mc_through(self, max_stage, is_signal):
        """Increment scalar MC counters for stages 0..max_stage."""
        end = min(max_stage + 1, self._n_stages)
        if is_signal:
            self._n_mc_signal[:end] += 1
        else:
            self._n_mc_background[:end] += 1

    def fill_eff_numerator_through(self, max_stage, true_bin_index):
        """Fill the efficiency numerator into stages 0..max_stage (inclusive).

        Call once per MC event that is truth-signal-in-phase-space, with
        `max_stage` = the number of reco cuts the event passes and
        `true_bin_index` = the bin of its TRUE kinematic value. The event is
        recorded in every stage up to and including the last reco cut it
        survives, so the per-stage efficiency (eff_numerator / truth) tracks
        acceptance loss cut-by-cut. Out-of-range true bins are dropped (they
        are likewise excluded from the truth denominator).
        """
        if true_bin_index < 0 or true_bin_index >= self._n_bins:
            return
        end = min(max_stage + 1, self._n_stages)
        self._eff_numerator[:end, true_bin_index] += 1

    def fill_data(self, stage_index, bin_index):
        """Fill a single data event into one stage."""
        if bin_index < 0 or bin_index >= self._n_bins:
            return
        self._data[stage_index, bin_index] += 1

    def fill_data_through(self, max_stage, bin_index):
        """Fill a single data event into stages 0..max_stage (inclusive)."""
        if bin_index < 0 or bin_index >= self._n_bins:
            return
        end = min(max_stage + 1, self._n_stages)
        self._data[:end, bin_index] += 1

    def count_data_through(self, max_stage):
        """Increment scalar data counters for stages 0..max_stage."""
        end = min(max_stage + 1, self._n_stages)
        self._n_data[:end] += 1

    def fill_truth(self, bin_index):
        """Fill a single truth event (no stage dimension)."""
        if bin_index < 0 or bin_index >= self._n_bins:
            return
        self._truth[bin_index] += 1

    def set_truth_total(self, n_total):
        """Set the total number of truth entries processed."""
        self._truth_total = n_total

    def register_cut_variable(self, name, edges, *, label, units, stage):
        """Declare a cut-variable distribution to accumulate.

        Experiment-agnostic: `edges` may be any strictly increasing list
        (use 2-bin edges such as [-0.5, 0.5, 1.5] for a boolean pass/fail
        cut). `stage` records the selection stage at which the histogram
        is filled (a free-text label, e.g. "entering_Apothem").
        """
        if name in self._cut_vars:
            raise ValueError(f"cut variable '{name}' already registered")
        edges = list(edges)
        if len(edges) < 2:
            raise ValueError(f"cut variable '{name}': edges need >= 2 values")
        n_bins = len(edges) - 1
        self._cut_vars[name] = {
            "label": label,
            "units": units,
            "stage": stage,
            "edges": edges,
            "mc_signal": np.zeros(n_bins),
            "mc_background": np.zeros(n_bins),
            "data": np.zeros(n_bins),
        }

    def fill_cut_variable(self, name, value, kind):
        """Fill one entry into a registered cut variable's histogram.

        Parameters
        ----------
        name : str
            Registered cut-variable name.
        value : float
            Observable value; binned against the variable's edges.
            Out-of-range values are dropped (overflow/underflow ignored).
        kind : str
            One of 'mc_signal', 'mc_background', 'data'.
        """
        cv = self._cut_vars.get(name)
        if cv is None:
            raise KeyError(f"cut variable '{name}' not registered")
        if kind not in ("mc_signal", "mc_background", "data"):
            raise ValueError(
                f"kind must be mc_signal/mc_background/data, got {kind!r}")
        edges = cv["edges"]
        if value < edges[0] or value >= edges[-1]:
            return
        bin_index = int(np.searchsorted(edges, value, side="right") - 1)
        cv[kind][bin_index] += 1

    def _python_stage_total(self, stage_index, quantity):
        """Python-side integrated total for one stage, for cross-check."""
        if quantity == "n_data":
            return float(self._n_data[stage_index])
        if quantity == "n_mc":
            return float(self._n_mc_signal[stage_index]
                         + self._n_mc_background[stage_index])
        if quantity == "n_mc_signal":
            return float(self._n_mc_signal[stage_index])
        if quantity == "n_mc_background":
            return float(self._n_mc_background[stage_index])
        raise ValueError(
            f"unknown cross-check quantity {quantity!r}; expected one of "
            "n_data / n_mc / n_mc_signal / n_mc_background")

    def attach_cutflow_crosscheck(self, external_totals, *, source,
                                  compare=("n_data",), rel_tol=1e-6):
        """Cross-check integrated per-stage totals against an external source.

        The external source is the C++ `runEventLoop` `mycuts` cut-summary
        table (per-stage event counts / % Eff / % Purity for MC, total
        entries for data), dumped to a `cutflow_summary_*.json` file. This
        method diffs those integrated totals against this accumulator's own
        scalar counters (`count_mc_through` / `count_data_through`), which are
        a cheap side-product of building the per-bin histograms. A
        disagreement means the Python companion selection and the compiled C++
        selection diverge — a real correctness bug, surfaced downstream by
        `diagnostic_plots.py` with the same severity as a failed sanity status.

        Parameters
        ----------
        external_totals : dict or list
            Either `{stage_name: {quantity: value, ...}, ...}` or a list of
            dicts each carrying a `name` (or `stage`) key plus quantity keys.
            Stage names must match this accumulator's `stage_names`.
        source : str
            Human-readable provenance (e.g. the cutflow JSON path).
        compare : tuple[str]
            Quantities to diff. Default `("n_data",)`: data counts are
            unweighted integers and must match exactly. Add `"n_mc"` only when
            the companion pass reproduces the C++ CV weighting (otherwise a
            weighted-vs-count difference would false-positive).
        rel_tol : float
            Relative tolerance for a match (default 1e-6, i.e. effectively
            exact for integer counts).

        Notes
        -----
        Only quantities present in a given external stage entry are compared;
        a compared stage name absent from this accumulator (or vice versa) is
        itself recorded as a mismatch (a structural divergence).
        """
        # Normalize external_totals to {stage_name: {quantity: value}}.
        norm = {}
        if isinstance(external_totals, dict):
            norm = {str(k): dict(v) for k, v in external_totals.items()}
        else:
            for entry in external_totals:
                entry = dict(entry)
                name = entry.pop("name", None) or entry.pop("stage", None)
                if name is None:
                    raise ValueError(
                        "each external_totals entry needs a 'name'/'stage' key")
                norm[str(name)] = entry

        comparisons = []
        for stage_name, quantities in norm.items():
            if stage_name not in self._stage_names:
                comparisons.append({
                    "stage": stage_name, "quantity": "(structure)",
                    "python": None, "external": None, "rel_diff": None,
                    "match": False,
                    "note": "external stage not present in Python accumulator",
                })
                continue
            si = self._stage_names.index(stage_name)
            for q in compare:
                if q not in quantities:
                    continue
                ext = float(quantities[q])
                py = self._python_stage_total(si, q)
                denom = max(abs(ext), 1e-9)
                rel_diff = abs(py - ext) / denom
                comparisons.append({
                    "stage": stage_name, "quantity": q,
                    "python": py, "external": ext,
                    "rel_diff": rel_diff, "match": rel_diff <= rel_tol,
                })

        status = "fail" if any(not c["match"] for c in comparisons) else "pass"
        self._cutflow_crosscheck = {
            "source": source,
            "rel_tol": rel_tol,
            "status": status,
            "comparisons": comparisons,
        }
        return self._cutflow_crosscheck

    def _cumulative_cuts(self, stage_index):
        if self._cut_labels is None:
            return []
        return self._cut_labels[:stage_index]

    def to_payload(self, *, experiment, channel, variable, variable_label,
                   units, signal_definition, pot_data, pot_mc):
        """Return a dict conforming to DiagnosticPayload schema v1.3."""
        stages = []
        for i in range(self._n_stages):
            stages.append({
                "name": self._stage_names[i],
                "cumulative_cuts": self._cumulative_cuts(i),
                "mc_signal": self._mc_signal[i].tolist(),
                "mc_background": self._mc_background[i].tolist(),
                "data": self._data[i].tolist(),
                "n_mc_signal": int(self._n_mc_signal[i]),
                "n_mc_background": int(self._n_mc_background[i]),
                "n_data": int(self._n_data[i]),
                "eff_numerator": self._eff_numerator[i].tolist(),
            })

        cut_variables = []
        for name, cv in self._cut_vars.items():
            cut_variables.append({
                "name": name,
                "label": cv["label"],
                "units": cv["units"],
                "stage": cv["stage"],
                "edges": cv["edges"],
                "mc_signal": cv["mc_signal"].tolist(),
                "mc_background": cv["mc_background"].tolist(),
                "data": cv["data"].tolist(),
            })

        payload = {
            "schema_version": SCHEMA_VERSION,
            "metadata": {
                "experiment": experiment,
                "channel": channel,
                "variable": variable,
                "variable_label": variable_label,
                "units": units,
                "signal_definition": signal_definition,
                "pot_data": pot_data,
                "pot_mc": pot_mc,
            },
            "binning": {"edges": self._edges},
            "truth": {
                "counts": self._truth.tolist(),
                "n_total": self._truth_total,
            },
            "stages": stages,
            "cut_variables": cut_variables,
        }
        if self._cutflow_crosscheck is not None:
            payload["cutflow_crosscheck"] = self._cutflow_crosscheck
        return payload
