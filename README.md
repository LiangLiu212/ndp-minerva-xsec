# Neutrino Discovery Platform (NDP)

**A theorist provides a model. The platform turns it into generator-level events, pushes them
through a learned detector surrogate, and compares them with neutrino-scattering data — in both
the experiment's unfolded space and the reconstructed-event space — with every physics choice in
a manifest and every number in a manifest-backed run directory.**

```
 theory model ──► generator sample ──► detector surrogate ──► comparison with data
 (YAML spec)      (TruthTable)         p(reco | truth)         unfolded: published d²σ + covariance
                                                               folded:   predicted reco counts vs data
```

| stage | what exists today | module |
|---|---|---|
| **theory → generator** | GENIE runs (any tune / custom tune dir) on the channel's flux and target mix, absolutely normalised from the spline file; reweighting of a reference MC by a formula or a Python function; external event files; shipped generator curves | `ndp/theory/` |
| **channel** | signal definition, true phase space, observables, binning, selection, data references, normalisation constants — one YAML per measurement | `channels/`, `ndp/channels/` |
| **surrogate** | binned response (efficiency × migration + background) and a parametric smearing model, both learned from the experiment's paired truth/reco MC and certified by exact closure | `ndp/surrogate/` |
| **data** | MINERvA Open Data AnaTuples (reco-level data counts) and the published cross-section releases with covariances via the certified MINERvA benchmark engine | `ndp/adapters/`, `ndp/compare/` |
| **orchestration** | `ndp run model.yaml --channel …` → `runs/<id>/{manifest,scorecard,report,figs}` | `ndp/pipeline.py`, `ndp/cli.py` |

First channel: **MINERvA ME FHC inclusive CC νμ, d²σ/dpT dp∥** (arXiv:2106.16210), with the
open-data me1A files and the paper's 224-cell grid.

## Quick start

```bash
cd ndp-platform
python -m ndp channels                                   # channels and their status
python -m ndp models                                     # example model specs (validated)
python -m ndp data status                                # which inputs are present
python -m ndp run models/reweight_mec_x1p5.yaml --channel minerva_me_cc_inclusive_ptpz
cat runs/<the new run dir>/report.md
```

Requirements: Python ≥ 3.10 with numpy, PyYAML, uproot, awkward, matplotlib (scipy optional;
pytest for the tests, or use `python tests/run_tests.py`). `pixi.toml` pins an environment if you
want one. Site paths are in `ndp.yaml` (or `NDP_*` environment variables, see `ndp/config.py`):

- `minerva_repo` — the `ndp-minerva-data-release-exploration` checkout (paper releases, the
  `benchmark/` χ² engine, the certified selection tool). Read-only from here.
- `data_dir` — the MINERvA open-data AnaTuples (`xrdcp` from
  `root://fndcadoor.fnal.gov:1095//pnfs/fnal.gov/usr/minerva/persistent/OpenData/...`; the
  `README` of the MINERvA repo lists the files). `data_dir/cache/` holds the vectorised truth/reco
  tables the surrogate builder reads; build them once with the snippet in `docs/architecture.md`.
- `genie_env_json` — a genie-agent environment snapshot for a GENIE installation.
- `genie_splines` — a GENIE cross-section spline XML covering the channel's nuclides.

## Writing a model

```yaml
name: mec_x1p5
kind: reweight                 # shipped_curve | reference_mc | reweight | genie | external
base: reference_mc
weight_expr: "where(int_type == MEC, 1.5, 1.0)"
description: 2p2h scaled by 1.5
```

```yaml
name: genie_G18_02a
kind: genie
tune: G18_02a_00_000
generator_list: CC
n_events: 180000
n_jobs: 6
# gxmlpath: path/to/custom/tunes   # a theorist's own GENIE tune directory
```

The weight expression sees every truth column (`E_nu`, `Q2`, `W`, `int_type`, `target_A`, …),
every observable (`lep_pT`, `lep_pz`, `q0`, `q3`, `E_avail`, …), the codes `QE RES DIS COH MEC`
and numpy. A Python function of the `TruthTable` works too (`weight_module: my_model.py:weight`).
External samples are read from GENIE `gst` files or the platform's `.npz` format.

## What a run tells you

`report.md` gives the model's sample summary, the **unfolded** rows (total χ²/ndf, shape χ²/ndf
with its profiled α, normalisation offset — beside every generator curve the paper shipped), the
**folded** result (data/prediction, −2lnL/ndf, Pearson χ²/ndf with MC statistics, which surrogate),
the figures, and every caveat. `scorecard.json` holds all numbers; `manifest.json` the inputs,
fingerprints, versions and git state.

## Documentation

- `docs/architecture.md` — data contracts, normalisation conventions, how to extend
- `docs/decisions.md` — platform defaults and the evidence behind them
- `docs/open_questions.md` — what the physicist still has to decide
- `docs/roadmap.md` — done / next
- `CLAUDE.md`, `.claude/skills/ndp-model/` — how the coding agent operates the platform
- `archive/ndp-minerva-xsec/` — the earlier pure-Python 2106.16210 cross-section reproduction with systematics (history preserved; see `archive/README.md`)
