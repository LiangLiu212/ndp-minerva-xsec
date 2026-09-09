# Neutrino Discovery Platform (NDP)

**A theorist provides a model. The platform turns it into generator-level events, pushes them
forward through a detector surrogate learned from the experiment's official MC, and compares them
with the experiment's reconstructed data — on the experiment's published grid or on any observable
the analyst defines — with every physics choice in a manifest and every number in a manifest-backed
run directory. The model is folded to the data; the data are never unfolded here.**

```
 theory model ──► generator sample ──► detector surrogate ──► comparison with data
 (YAML spec)      (TruthTable)         p(reco | truth),        folded:   predicted reco counts vs selected data,
                                       learned from the                  on the published grid or a user measurement
                                       official MC             unfolded: published d²σ + covariance (context, where it exists)
```

Three inputs play three roles: **data** is what was recorded (reco-level candidates, plus the
published cross section where one exists); the **official MC** only teaches the surrogate what the
detector does to true signal; the **theorist's model** is what gets tested.

| stage | what exists today | module |
|---|---|---|
| **theory → generator** | GENIE runs (any tune / custom tune dir) on the channel's flux and target mix, absolutely normalised from the spline file; reweighting of a reference MC by a formula or a Python function; external event files; shipped generator curves | `ndp/theory/` |
| **channel** | signal definition, true phase space, observables, binning, selection, data references, normalisation constants — one YAML per measurement | `channels/`, `ndp/channels/` |
| **measurement** | the observable pair a model is compared on — truth side and reco side, each a registry name or an expression — with its binning; the channel's published grid plus any number of user-defined ones (`measurements/<channel>/*.yaml`) | `measurements/`, `ndp/channels/measurements.py` |
| **surrogate** | binned response (efficiency × migration + background) and a parametric smearing model, learned per measurement from the experiment's paired truth/reco MC and certified by exact closure | `ndp/surrogate/` |
| **data** | MINERvA Open Data AnaTuples (reco-level data counts) and the published cross-section releases with covariances via the certified MINERvA benchmark engine | `ndp/adapters/`, `ndp/compare/` |
| **orchestration** | `ndp run model.yaml --channel …` → `runs/<id>/{manifest,scorecard,report,figs}` | `ndp/pipeline.py`, `ndp/cli.py` |

First channel: **MINERvA ME FHC inclusive CC νμ, d²σ/dpT dp∥** (arXiv:2106.16210), with the
open-data me1A files and the paper's 224-cell grid.

## Quick start

```bash
cd ndp-platform
python -m ndp channels                                   # channels, their status and measurements
python -m ndp models                                     # example model specs (validated)
python -m ndp data status                                # which inputs / caches are present
python -m ndp run models/reweight_mec_x1p5.yaml --channel minerva_me_cc_inclusive_ptpz
cat runs/<the new run dir>/report.md
```

### Your own observable (forward folding on a user measurement)

```bash
python -m ndp data cache --channel minerva_me_cc_inclusive_ptpz          # once: reco/truth tables from the AnaTuples
python -m ndp measurements --channel minerva_me_cc_inclusive_ptpz         # published + user measurements, surrogate status
python -m ndp surrogate build --channel minerva_me_cc_inclusive_ptpz --measurement muon_p_theta
python -m ndp run models/genie_g18_02a_me_tracker.yaml --channel minerva_me_cc_inclusive_ptpz --measurement muon_p_theta
```

A measurement is a small YAML (`measurements/<channel>/<name>.yaml`):

```yaml
name: enu_calorimetric
channel: minerva_me_cc_inclusive_ptpz
x: {observable: E_nu, reco: "reco_E_mu + reco_recoil_E", units: GeV, edges: [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 30, 50]}
# y: {observable: lep_theta_deg, reco: reco_theta_deg, units: deg, edges: [...]}   # omit for 1D
```

`observable` is a truth quantity (registry name or expression over the truth columns), `reco` a
reconstructed quantity (registry name or expression over the cached reco columns —
`python -m ndp data status` shows the cache; `ndp/channels/reco_observables.py` the names). The
surrogate is learned on exactly this grid from the official MC and certified by closure; the
comparison is folded only, since nothing was published on it.

Requirements: Python ≥ 3.10 with numpy, PyYAML, uproot, awkward, matplotlib (scipy optional;
pytest for the tests, or use `python tests/run_tests.py`). Site paths are in `ndp.yaml` (or
`NDP_*` environment variables, see `ndp/config.py`):

- `minerva_repo` — the `ndp-minerva-data-release-exploration` checkout (paper releases, the
  `benchmark/` χ² engine, the certified selection tool). Read-only from here.
- `data_dir` — the MINERvA open-data AnaTuples (`xrdcp` from
  `root://fndcadoor.fnal.gov:1095//pnfs/fnal.gov/usr/minerva/persistent/OpenData/...`; the
  `README` of the MINERvA repo lists the files). `data_dir/cache/` holds the vectorised truth/reco
  tables the surrogate builder reads; build them once with the snippet in `docs/architecture.md`.
- `genie_env_json` — a genie-agent environment snapshot for a GENIE installation.
- `genie_splines` — a GENIE cross-section spline XML covering the channel's nuclides.

## Environment and the in-repository GENIE

`pixi.toml` describes the complete environment from conda-forge — the Python stack plus ROOT,
GCC/gfortran, GSL, log4cpp, libxml2 and LHAPDF 6 — and `activate.sh` exports the GENIE variables
whenever you `pixi shell` or `pixi run`:

```bash
pixi install                 # ~3 GB under .pixi/
pixi run build-pythia6       # ROOTEGPythia6: Pythia6 + ROOT's TPythia6 interface (conda ROOT ships neither)
pixi run build-genie         # GENIE Generator $GENIE_VERSION (default R-3_06_02) built in external/genie/Generator
pixi run snapshot-genie-env  # -> external/genie_env.json, the environment the pipeline hands to gevgen/gntpc
pixi run test-genie          # 50-event gevgen + gntpc smoke test
pixi run make-splines -- --tune <tune>   # cross-section splines for a (custom) tune, one gmkspl per nuclide
pixi run test                # platform tests under pixi's pytest
```

The other generators build the same way, into the same environment:

```bash
pixi run build-nuwro         # NuWro (GitHub master, 25.11.1 at the time of writing) with ROOTEGPythia6
pixi run build-gibuu         # GiBUU release 2025 + buuinput (hepforge tarballs; GIBUU_TARBALLS=<dir> reuses local copies)
pixi run build-achilles      # ACHILLES v0.3.1 (CMake; HDF5/zlib from pixi, the rest fetched by CPM)
pixi run build-generators    # all four
pixi run test-generators     # NuWro / GiBUU / ACHILLES smoke runs (GENIE: test-genie)
```

`activate.sh` exports `NUWRO`, `GIBUU`, `GIBUU_INPUT`, `ACHILLES` and puts every binary on `PATH`.
Their outputs enter the pipeline through the `external` model kind: `format: nuwro_root`
(`treeout`), `format: gibuu_finalevents` (`FinalEvents.dat`, needs `target_Z`/`target_A`) and
`format: nuhepmc` (ACHILLES / any NuHepMC file). Everything under `external/` is gitignored and
rebuildable from `scripts/`. `ndp.yaml`'s `genie_env_json` points at the snapshot, so `kind: genie` models run on
the in-repo build; `genie_splines` names the cross-section splines (the CVMFS `G18_02a_00_000` set
by default when mounted). A custom tune directory goes in via a model's `gxmlpath`. Any other GENIE
install works the same way: point `genie_env_json` at a snapshot of its environment.

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

`report.md` names the measurement, gives the model's sample summary, the **folded** result
(data/prediction, −2lnL/ndf, Pearson χ²/ndf with MC statistics, which surrogate), then — for the
published grid — the **unfolded** rows (total χ²/ndf, shape χ²/ndf with its profiled α, normalisation
offset, beside every generator curve the paper shipped), the figures, and every caveat. `scorecard.json` holds all numbers; `manifest.json` the inputs,
fingerprints, versions and git state.

## Documentation

- `docs/architecture.md` — data contracts, normalisation conventions, how to extend
- `docs/decisions.md` — platform defaults and the evidence behind them
- `docs/open_questions.md` — what the physicist still has to decide
- `docs/roadmap.md` — done / next
- `CLAUDE.md`, `.claude/skills/ndp-model/` — how the coding agent operates the platform
- The earlier pure-Python 2106.16210 cross-section reproduction with systematics (`ndp-minerva-xsec`) lives in this repository's history: tag `ndp-minerva-xsec-final`, merged at `8a1785a`
