# Roadmap — where the vision stands

Vision: *a theorist provides a model; the platform converts it to a generator sample, smears the
generator truth into detector reconstruction with a surrogate model, and compares with data; theory
→ generator, channel, surrogate and data are all implemented from manifests.*

## Done (2026-09-04)

| stage | status | evidence |
|---|---|---|
| theory → generator | GENIE runs on the channel flux/target with absolute normalisation (180k-event G18_02a CC sample in 230 s on 6 cores); reweighting by formula or Python; external gst/npz samples; shipped curves | `runs/2026-09-04_genie_G18_02a_00_000__*`, `models/` |
| channel | `minerva_me_cc_inclusive_ptpz` fully specified with per-field status; low-recoil channel drafted | `channels/` |
| surrogate | binned response + parametric smearing learned from the me1A MC, exact closure, beam-frame truth | `surrogates/minerva_me_cc_inclusive_ptpz/`, tests |
| data | reco-level open data (844 selected / 835 in grid, 2.05e17 POT) and the 2106.16210 release + covariance through the certified benchmark engine (Tune v1 row 33.03 reproduced) | `runs/2026-09-04_MINERvA_Tune_v1_shipped__*_2` |
| orchestration | `ndp run` → manifest-backed run dir, report, figures; CLI; agent skill; 22 tests | `ndp/pipeline.py`, `tests/` |

First results (me1A slice, statistical errors only, see the run directories for every number):

| model | unfolded χ²_total/ndf (paper cov / +model stat) | norm offset | folded −2lnL/ndf | data/pred |
|---|---|---|---|---|
| MINERvA Tune v1 (shipped) | 33.03 / — | −10.5 % | — | — |
| reference MC (GENIE 2.12.6 CV, POT-normalised) | 115.9 / 5.39 | +6.0 % | 205.0/204 | 0.938 |
| reference MC, MEC × 1.5 | 132.7 / — | +8.2 % | 209.1/204 | 0.917 |
| reference MC, QE low-Q² suppression | 101.9 / — | +4.8 % | 204.5/204 | 0.952 |
| GENIE 3.6.2 G18_02a_00_000 (absolute) | 240.6 / 22.19 | −4.3 % | 215.6/203 | 1.063 |

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
   the non-Gaussian tails and to smear on arbitrary binnings.
5. **Second channel.** Finish `minerva_me_lowrecoil_eavail_q3`: exact edges from Table II, a
   release-manifest kind for the 44-bin covariance, the reco E_avail estimator (exploration repo
   feasibility probe), then the hadronic surrogate.
6. **Theory → generator, richer.** GENIE Reweight knobs as a model kind; custom-tune runs through
   `gxmlpath` (supported, untested); NuWro/GiBUU adapters (the user's nc1p workspace has both).
7. **Housekeeping.** Push the repository, `pixi install` when approved, pytest in CI, a `runs/`
   index generator.
