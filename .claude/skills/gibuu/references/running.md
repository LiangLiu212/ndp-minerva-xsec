# Part 2 — Running GiBUU here (job cards, fluxes, targets, ensembles, weights)

Verified on the FNAL EAF node on 2026-09-15 with the platform install (Release 2025 patch 5,
`external/gibuu/release2025/objects/GiBUU.x`, tables `external/gibuu/buuinput2025`). Every number below
was measured in this workspace; the commands are the ones that produced it. Project mechanics live in
`ndp/theory/gibuu.py` (card writer, local runner, merge), `resources/gibuu/minerva_me_numu_CC.job.tmpl`
(the card template), `models/gibuu_2025_me_fhc_c12.yaml` (a model spec of `kind: gibuu`),
`ndp/grid/gibuu_campaign.py` + `grid/gibuu_worker.sh` (grid campaigns; recipe in `grid/README.md`).

## 2.1 How a run is driven

`GiBUU.x < job.card` in an empty directory: GiBUU writes every output into the current directory
(`FinalEvents.dat`, the `neutrino_absorption_cross_section_*.dat` files, ~150 diagnostic `.dat`
files, `main.run`, `PYR.RG`). The card is a set of Fortran namelists; a namelist that is absent takes
GiBUU's defaults (the log says `namelist "xxx" missing in jobcard. Values-->Default!`). The platform
never writes cards by hand: `ndp.theory.gibuu.write_card` fills the template from the model spec
(`GibuuSpec.from_params`) and refuses to leave a `__PLACEHOLDER__` behind. `ndp gibuu smoke
<model.yaml> --channel <c> --ensembles N` runs one local job under
`runs/_generator_cache/gibuu_<fingerprint>/jobs/<name>/`.

The fingerprint of a sample = sha256 of the spec + the channel flux table path + the template file
(`GibuuSpec.fingerprint`), so a changed card or flux is a new cache directory.

## 2.2 The card the platform uses (from GiBUU's own MINERvA-ME card)

Template `resources/gibuu/minerva_me_numu_CC.job.tmpl`, trimmed from
`$GIBUU/testRun/jobCards/005_Neutrino_MINERvA_ME-nu.job`. Namelists and what the spec controls:

| namelist | keys | spec / default |
|---|---|---|
| `&neutrino_induced` | `process_ID` (2 CC), `flavor_ID` (2 numu), `nuXsectionMode = 16` (EXP_dSigmaMC, flux-averaged), `nuExp = 99` + `FileNameFlux`, `includeQE/DELTA/RES/1pi/DIS/2p2hQE/2pi`, `printAbsorptionXS = T`, `equalWeights_Mode`, `equalWeights_Max` | `process_ID`, `flavor_ID`, `include`, `card.equal_weights_*` |
| `&target` | `Z`, `A`, `densitySwitch_Static` (2), `fermiMotion` | `target`, `card.density_switch_static` |
| `&input` | `eventtype = 5`, `numEnsembles`, `numTimeSteps`, `delta_T`, `num_runs_SameEnergy`, `path_to_input`, `version = 2025` | `card.num_ensembles / num_time_steps / delta_T / num_runs`; the input path is filled at run time |
| `&width_Baryon` | `mediumSwitch`, `mediumSwitch_Delta`, `mediumSwitch_coll` | `card.medium_switch*` (GiBUU's card: T, T, F) |
| `&baryonPotential` | `EQS_Type` (5), `DeltaPot` (1) | `card.eqs_type`, `card.delta_pot` |
| `&neutrinoAnalysis` | `outputEvents = .true.` (writes `FinalEvents.dat`), detection thresholds 0, `applyCuts` (2: only unbound particles) | `card.apply_cuts` |
| `&initRandom` | `SEED` | the job seed (a 32-bit Fortran integer; job k of a campaign uses `seed + k`) |

Facts checked in this session:
- **`numEnsembles < 100` aborts** ("It is a bad idea to use #ensembles < 100 ! To enforce it anyway,
  give the number with a minus sign."). Smoke runs use 100.
- **Events per job are fewer than A × numEnsembles**: 1000 ensembles on C12 gave 5 759 events
  (0.48 × 12 000), 100 ensembles 588; the rest fall in forbidden phase space (GiBUU's card comment).
- **Timing / memory on the EAF node (1 core):** 100 ensembles 32 s, 1000 ensembles 88 s
  (~65 events/s), peak RSS 0.45 GB. A 4000-ensemble job is a few minutes.
- The in-medium width switches (`mediumSwitch`, `mediumSwitch_Delta`) do not change the QE cross
  section (identical σ_QE with them on and off at fixed energy); they change the resonance channels.

## 2.3 The flux (`nuExp = 99`)

`code/init/neutrino/esample.f90` reads the file named in `FileNameFlux`: lines starting with `#` are
comments (only at the top), then two columns `E [GeV]  flux` where `E` is the **bin centre** of an
**equidistant** grid; the energy is sampled by inverting the cumulative sum and drawing uniformly inside
the bin, so the second column is a bin-averaged density in arbitrary units. A name containing `/` is
taken as an absolute path; otherwise the file is looked up in `<path_to_input>/neutrino/`.

The platform writes the channel's flux table (MINERvA ME FHC, `resources/flux/arXiv2110.13372_supplemental.txt`
Table I, non-uniform bins) rebinned to uniform 0.5 GeV bins over 0–100 GeV (`ndp.theory.flux.rebin_uniform`,
integral conserved) as `runs/_generator_cache/gibuu_<fp>/flux_gibuu.dat`. GiBUU's own
`buuinput/neutrino/MINERvA_MEflux.dat` is a different digitisation of the ME flux (0.5 GeV bins,
0.25–94.75 GeV) and is not used, so that the flux-averaged cross section refers to the same flux and
range as the channel's `phi_per_pot_cm2`. Flux mean energy of the rebinned table: 5.95 GeV; the
σ-weighted mean energy of the generated events 6.15–6.56 GeV.

## 2.4 Weighted events and the equal-weights mode

In the default mode (`equalWeights_Mode = 0`) every FinalEvents row of an event carries the event's
`perweight` = σ(that test nucleon, sampled E)/(A × numEnsembles) in 10⁻³⁸ cm² per nucleon, and the
perweights of a run sum to the flux-averaged σ_CC per nucleon (checked against
`neutrino_absorption_cross_section_ALL.dat` column 2 to 5 × 10⁻⁶ relative on every run of this
session). Under the ME flux the weights are extremely non-uniform because the lepton kinematics are
sampled uniformly and QE is sharply forward-peaked: on the 1000-ensemble job (5 759 events, σ_CC =
4.078 × 10⁻³⁸ cm²/nucleon)

| channel | events | Σw (10⁻³⁸ cm²) | effective events Σw²/Σw² | max w / mean w |
|---|---|---|---|---|
| QE | 266 | 0.184 | 11 | 56 |
| RES | 1 073 | 0.682 | 454 | 9 |
| DIS | 3 420 | 3.082 | 1 390 | 18 |
| 2p2h | 1 000 | 0.130 | 260 | 13 |
| 1μ1p signal (all) | 328 | 0.240 | 21 | 53 |

so a weighted GiBUU sample is useless for a QE-dominated selection: 328 signal events carry the
statistical power of 21 (the QE part: 4). Per-nucleon weights (× A × numEnsembles): mean 8.5, median
5.1, 99 % 54, 99.9 % 160, max 463 in 5 759 events.

`equalWeights_Mode = 2` applies MC rejection against `equalWeights_Max` (per-nucleon σ units): an
event is kept with probability w/Max and then carries perweight Max/(A × numEnsembles), i.e. all
weights equal; the acceptance is ⟨w⟩/Max. **GiBUU aborts the run if an event's weight exceeds
`equalWeights_Max`** (`initNeutrino.f90`: "You have to increase 'equalWeights_Max'"), so the ceiling
must be set from a pilot in mode 1 (which prints the running maximum to unit 87 / `fort.87`) with a
safety margin; the cost is acceptance. Mode 2 measurements: see the "Gotchas" and the campaign record
`grid/campaigns/gibuu_*/campaign.json` (filled in when the campaign runs).

## 2.5 Grid campaigns

`grid/README.md` ("GiBUU generation campaigns"): `ndp gibuu plan` writes the card template with the
three run-time placeholders and the flux file into the cache directory, `grid/stage_gibuu_payload.sh`
assembles the payload (binary + `libgfortran.so.5`, `libquadmath.so.0`, `libgcc_s.so.1`,
`libbz2.so.1.0` from the pixi env per `ldd` + `buuinput` + cards) and checks it with a scrubbed
`ldd`, the jobsub-lite skill publishes it, `ndp gibuu submit-cmd` prints the submit line
(`gibuu_worker.sh -R payload -O pnfs_out -C card -S seed_base`, seed = base + PROCESS), and
`status / harvest / merge` bring the jobs back. The merge prescales weights by 1/K over K jobs and
refuses a job whose Σw disagrees with its absorption file by more than 10⁻³.

## Gotchas

- GiBUU writes to the cwd: never run two jobs in one directory.
- A card with `numEnsembles < 100` aborts unless the number is negative.
- The `__PLACEHOLDER__` convention: any `__NAME__` token left in a card is an error (the writer checks;
  a comment must not contain such a token).
- `equalWeights_Mode = 2` aborts on a weight above the ceiling — losing the job; pilot the maximum
  first and keep a margin.
- `nuXsectionMode = 6` (fixed energy) with only ~1000 test nucleons gives cross sections that scatter
  by factors of 2–5 between runs for QE (few events with large weights); do not read channel cross
  sections from small fixed-energy runs.
