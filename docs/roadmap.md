# Roadmap — where the vision stands

Vision: *a theorist provides a model; the platform converts it to a generator sample, smears the
generator truth into detector reconstruction with a surrogate model learned from the experiment's
official MC, and compares with data — the published cross section where one exists, or any
observable the analyst defines; theory → generator, channel, measurement, surrogate and data are all
implemented from manifests. The model is folded forward to the data; the data are never unfolded.*

## Done (2026-09-04)

| stage | status | evidence |
|---|---|---|
| theory → generator | GENIE runs on the channel flux/target with absolute normalisation (180k-event G18_02a CC sample in 230 s on 6 cores); reweighting by formula or Python; external gst/npz samples; shipped curves | `runs/2026-09-04_genie_G18_02a_00_000__*`, `models/` |
| channel | `minerva_me_cc_inclusive_ptpz` fully specified with per-field status; low-recoil channel drafted | `channels/` |
| surrogate | binned response + parametric smearing learned from the me1A MC, exact closure, beam-frame truth | `surrogates/minerva_me_cc_inclusive_ptpz/`, tests |
| data | reco-level open data (844 selected / 835 in grid, 2.05e17 POT) and the 2106.16210 release + covariance through the certified benchmark engine (Tune v1 row 33.03 reproduced) | `runs/2026-09-04_MINERvA_Tune_v1_shipped__*_2` |
| orchestration | `ndp run` → manifest-backed run dir, report, figures; CLI; agent skill; 24 tests | `ndp/pipeline.py`, `tests/` |
| environment | pixi project (conda-forge) with ROOT 6.40, GCC 15.2, GSL/log4cpp/libxml2/LHAPDF 6 and an in-repository GENIE R-3_06_02 (Pythia6 via ROOTEGPythia6); the pipeline drives it through `external/genie_env.json` (12k-event check: data/pred 1.065, same σ_avg as the spack build) | `pixi.toml`, `scripts/`, `runs/2026-09-04_genie_G18_02a_inrepo_12k__*` |
| generators | NuWro 25.11.1, GiBUU release 2025 (patch 5) + buuinput, ACHILLES v0.3.1 (e02d266) built in the same pixi environment (`pixi run build-generators`), smoke-tested, and readable through the `external` model formats `nuwro_root` / `gibuu_finalevents` / `nuhepmc` | `scripts/build_*.sh`, `ndp/adapters/`, `tests/test_generator_adapters.py` |

## Done (2026-09-09) — forward folding on user-defined observables

| stage | status | evidence |
|---|---|---|
| measurement | `Measurement` = (truth observable ↔ reco observable) × 2 + edges [+ release]; registry names or expressions on both sides; 1D or 2D; the channel's `published` grid is one of them, user ones live in `measurements/<channel>/` | `ndp/channels/measurements.py`, `reco_observables.py`, `measurements/minerva_me_cc_inclusive_ptpz/{muon_p_theta,enu_calorimetric,q2_calorimetric}.yaml` |
| reco cache | versioned per-AnaTuple reco tables (muon, MINOS, recoil-energy family, analysis-tool E_ν/W/x/y, visible E, multiplicities, vertex) via `ndp data cache`; 844 / 43643 selected reproduced; the published-grid surrogate rebuilt from the new cache is bit-identical to the tracked one | `ndp/adapters/minerva_anatuple.py::build_cache`, `ndp data status` |
| surrogate per measurement | `ndp surrogate build --measurement` (binned + parametric) with closure printed; binned closure exact (max deviation 0) on all three example measurements; parametric total within 0.2 % of the selected signal | `ndp/surrogate/build.py`, `surrogates/minerva_me_cc_inclusive_ptpz/<measurement>/`, `tests/test_measurements.py` |
| folded-first pipeline | `ndp run --measurement`; report leads with the folded comparison, unfolded only where the grid was published; run dir carries `measurement.json`; regression: `reference_mc` on the published grid reproduces the 2026-09-04 scorecard exactly (folded 205.0/204, data/pred 0.938; unfolded 115.9 / +6.0 %) | `runs/2026-09-09_reference_mc_genie2126__minerva_me_cc_inclusive_ptpz` |

Forward-folded results on the user measurements (me1A slice, statistical errors only; surrogate = binned response from MC 110040):

| model | measurement | −2lnL/ndf | Pearson χ²/ndf (MC stat) | data/pred | run |
|---|---|---|---|---|---|
| reference MC (POT-normalised) | muon_p_theta (14×9) | 163.8/126 | 186.4/126 | 0.944 | `2026-09-09_reference_mc_genie2126__*__muon_p_theta` |
| GENIE 3.6.2 G18_02a (absolute) | muon_p_theta | 166.3/126 | 187.6/126 | 1.067 | `2026-09-09_genie_G18_02a_00_000__*__muon_p_theta` |
| GENIE 3.6.2 G18_02a (absolute) | enu_calorimetric (1D, 13) | 26.0/13 | 27.4/13 | 1.031 | `2026-09-09_genie_G18_02a_00_000__*__enu_calorimetric` |
| reference MC, MEC × 1.5 | enu_calorimetric | 54.9/13 | 54.6/13 | 0.924 | `2026-09-09_mec_x1p5__*__enu_calorimetric` |
| GENIE 3.6.2 G18_02a (absolute) | q2_calorimetric (1D, 14) | 15.7/14 | 15.9/14 | 1.046 | `2026-09-09_genie_G18_02a_00_000__*__q2_calorimetric` |

The two calorimetric measurements use `MasterAnaDev_recoil_E`; which recoil-energy branch is
MINERvA-canonical is an open question (carried as a note into every run), so they demonstrate the
mechanism rather than endorse the estimator.

First results (me1A slice, statistical errors only, see the run directories for every number):

| model | unfolded χ²_total/ndf (paper cov / +model stat) | norm offset | folded −2lnL/ndf | data/pred |
|---|---|---|---|---|
| MINERvA Tune v1 (shipped) | 33.03 / — | −10.5 % | — | — |
| reference MC (GENIE 2.12.6 CV, POT-normalised) | 115.9 / 5.39 | +6.0 % | 205.0/204 | 0.938 |
| reference MC, MEC × 1.5 | 132.7 / — | +8.2 % | 209.1/204 | 0.917 |
| reference MC, QE low-Q² suppression | 101.9 / — | +4.8 % | 204.5/204 | 0.952 |
| GENIE 3.6.2 G18_02a_00_000 (absolute) | 240.6 / 22.19 | −4.3 % | 215.6/203 | 1.063 |
| GENIE 3.6.2 G18_10a_02_11b (absolute, in-repo GENIE + own splines, 60k events) | 207.7 / 11.49 | −9.4 % | 230.0/203 | 1.132 |

Flux-averaged CC cross sections per nucleon from the spline files: G18_02a 4.49e-38 cm², G18_10a 4.32e-38 cm²
(the latter from splines generated in-repo with 100 knots to 100 GeV).

The paper-covariance-only χ² of the truth-sample models is inflated by their own MC statistics
(1/500 of the paper's exposure for the reference MC); the "+model stat" column is the like-for-like
number. The folded comparison has ~35 % statistical error per cell at this exposure.

## Next

1. **Ratify the defaults** in `docs/open_questions.md` (beam frame, Φ for POT-normalised MC,
   feed-in treatment, n_nucleons, target mix). Each is a one-line change in the channel YAML.
2. **Statistics.** Stream/download the full Playlist 1A (253 data + 41 MC files, ~36× the data,
   ~8× the MC) or run the adapters on the grid via the exploration repo's jobsub-lite skill; the
   adapters are already vectorised and cache per file.
3. **Systematics in the surrogate.** Alternative responses from the tuple's flux/GENIE/detector
   weight universes; a covariance for the folded prediction; a MC-stat band on the figures.
   The earlier `ndp-minerva-xsec` code (tag `ndp-minerva-xsec-final` in this repository's history)
   already builds these universes (flux, 56 GENIE knobs, muon energy scale, RPA, 2p2h, MINOS
   efficiency) and validates them against the release covariances — the natural source to port from.
4. **Learned surrogate.** A conditional normalising flow (or diffusion) trained on the same paired
   MC, implementing `SmearingSurrogate`'s interface (`fit`, `sample_reco`, `fold_events`), to carry
   the non-Gaussian tails and to smear *once* for every measurement (today each measurement gets its
   own binned response; an event-level surrogate conditioned on the full truth kinematics would serve
   all of them). No torch/sklearn in the environment yet — adding one is a `pixi.toml` change.
4b. **Measurements, richer.** Systematic variations of the reco observable (e.g. the recoil-energy
   family) as alternative surrogates; a surrogate built with the exploration repo's canonical
   recoil-E once that question is settled; hadronic truth observables (E_avail, q3) need the
   final-state list, which the reference MC has and generator samples carry through the adapters.
5. **Second channel.** Finish `minerva_me_lowrecoil_eavail_q3`: exact edges from Table II, a
   release-manifest kind for the 44-bin covariance, the reco E_avail estimator (exploration repo
   feasibility probe), then the hadronic surrogate.
6. **Theory → generator, richer.** GENIE Reweight knobs as a model kind; custom-tune runs through
   `gxmlpath` (supported, untested). NuWro, GiBUU 2025 and ACHILLES are built in the pixi
   environment and their outputs are readable (`external` model formats); what is missing are
   *runners* (`kind: nuwro | gibuu | achilles` that write the card, run on the channel flux and
   target, and normalise) mirroring `ndp/theory/generator.py`.
7. **Housekeeping.** Push the repository, `pixi install` when approved, pytest in CI, a `runs/`
   index generator.
