# ndp-minerva-xsec

Reproduction of the MINERvA ME FHC inclusive CC ν_μ cross section
(arXiv:2106.16210) from the MINERvA Open Data MasterAnaDev AnaTuples,
in pure Python (uproot / awkward / numpy; PyROOT available).

## Workflow (four stages)

1. **Input files** — dataset specs in `config/datasets/*.json` declare every external
   input (AnaTuples, flux/reweight files) with its xrootd URL and expected streamed
   fingerprints (POT, entry counts); `summarize_inputs.py` verifies and summarizes
   them (POT, data-taking time span, branch contract). No path lives in analysis code.
2. **Cuts** — `xsec/cuts.py` (reco selection, 6 cuts), `xsec/signal.py` (truth signal
   definition + phase space), `xsec/kinematics.py` (θ3d, p_T, p_∥, beam rotation).
3. **Raw ingredients** — `make_ingredients.py`: one event-loop pass per dataset fills
   the six cross-section ingredients (data hist, background, migration matrix,
   efficiency numerator/denominator, POT ledger) into an intermediate `.npz`.
4. **Unfolding + normalization** — `extract_xsec.py`: background subtraction →
   D'Agostini unfolding → efficiency → flux Φ, targets T, POT, bin widths → dσ;
   `compare_published.py` validates against the paper's ancillary tables.

POT and target count T are scalar normalizations outside the loop.

## Environment

Managed by the parent pixi workspace at `/home/feanor/ndp-genesis-agent/pixi.toml`
(auto-discovered — run `pixi run python …` from anywhere inside this directory).
Python 3.14 / ROOT 6.40 / uproot / awkward / numpy / xrootd.

Run metadata: every finalized script wraps its body in `runlog_tools.RunLog`
(JSON logs in `~/log/ndp-minerva-xsec/`).

Physics constants (PDG IDs, masses, charges) come from the official `pdg` package —
never hardcoded. MINERvA-specific constants live in `config/constants.py` with
provenance comments.

## Data

**Streaming-only:** AnaTuples are never copied to local disk. Everything is read
via xrootd streaming (`root://fndcadoor.fnal.gov:1095/...`), with uproot fetching
only the branches the analysis declares in `config/branches.json`. File identity is
verified by streamed fingerprints (POT, entry counts), not checksums. The single
local-fetch exception is the FluxAndReweightFiles tarball (consumed at the weights
stage). Stage-3 intermediate `.npz` files are the only local data artifacts.

Reference material: certified selection, frozen 1D dσ/dp_T result, and the paper's
ancillary answer key are in the sibling repo
`../ndp-minerva-data-release-exploration` (see `papers/minerva/2106.16210/anc/`).

## Plans

Canonical plans are committed under `.claude/plans/`. Execution is stepwise:
each step runs only on explicit user go-ahead and ends with its gate evaluated.
