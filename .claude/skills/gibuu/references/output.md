# Part 3 — GiBUU output: `FinalEvents.dat`, the cross-section files, the platform adapter

Verified on the FNAL EAF node on 2026-09-15 (Release 2025 patch 5, smoke and 1000-ensemble jobs of
`models/gibuu_2025_me_fhc_c12.yaml`). Adapter: `ndp/adapters/gibuu_finalevents.py::read_finalevents`;
merge: `ndp/theory/gibuu.py::read_job / merge_jobs`.

## 3.1 `FinalEvents.dat`

Written when `&neutrinoAnalysis outputEvents = .true.`; one row per outgoing particle, 15 columns
(header line `# 1:Run 2:Event 3:ID 4:Charge 5:perweight 6:position(1) 7:position(2) 8:position(3)
9:momentum(0) 10:momentum(1) 11:momentum(2) 12:momentum(3) 13:history 14:production_ID 15:enu`):

| column | meaning | unit |
|---|---|---|
| 1, 2 | run (`num_runs_SameEnergy` counter), event number within the run | – |
| 3, 4 | GiBUU particle ID and charge (`gibuu_pdg()` maps them to PDG: 1 → p/n by charge, 101 → π, 901–913 leptons, 32/33 Λ/Σ, 2 Δ) | – |
| 5 | `perweight` of the event, repeated on every row; 0 on the struck nucleon (and the incoming lepton with `outputEvents_lepIn`) | 10⁻³⁸ cm² per nucleon |
| 6–8 | production position | fm |
| 9–12 | (E, px, py, pz) | GeV; the neutrino is along +z, so the sample is already in the beam frame (`meta["frame"] = "beam"`) |
| 13 | history code | – |
| 14 | `production_ID` of the first interaction: 1 QE, 2–31 resonances, 32/33 π background, 34 DIS, 35/36 2p2h, 37 2π | mapped to `int_type` 1 QE / 2 RES / 3 DIS / 5 MEC (32/33/37 → DIS, the GENIE convention; the raw code is kept in `gibuu_production_id`) |
| 15 | E_ν | GeV |

`applyCuts = 2` writes only unbound particles (the default in the platform card). The outgoing
lepton row is identified by ID 901–913 with a non-zero weight; the adapter drops the struck nucleon and
the lepton from the `fs_*` list and fills `lep_*`, `E_nu`, `Q2`, `W`, `current`, `int_type`,
`target_Z/A`, `weight = perweight / n_runs`.

Sizes: 588 events → 0.64 MB, 5 759 events → ~6 MB; the grid worker gzips the file.

## 3.2 Cross-section files

`neutrino_absorption_cross_section_ALL.dat`: header `# 1:var 2:sum 3:QE 4:Delta 5:highRES 6:1pi 7:DIS
8:2p2h-QE 9:2p2h-Delta 10:2pi`, one row per run; for the flux-averaged mode `var` is the
σ-weighted mean neutrino energy and `sum` the flux-averaged σ_CC per nucleon in 10⁻³⁸ cm². **The
summed perweights of a run equal column 2** (5 × 10⁻⁶ relative on every run checked); the platform's
`read_job` refuses a job where they differ by more than 10⁻³. Per-channel files
`neutrino_absorption_cross_section_{QE,Delta,highRES,1pi,DIS,2p2h,2pi}.dat` and
`_numbers.dat` (event counts) sit next to it. `neutrino_initialized_energyFlux.dat` is the histogram
of the sampled neutrino energies (a check of the flux file: 1 200 entries for 100 ensembles × 12).

Measured with the MINERvA ME FHC flux on C12 (100 / 1000 ensembles): σ_CC = 4.16 / 4.08 × 10⁻³⁸
cm²/nucleon (GENIE G18_02a splines: 4.66 at 6 GeV, 4.49 flux-averaged on the tracker mix);
1μ1p signal fraction 5.9 % of σ_CC (GENIE truth on the official MC: 5.5 % of νμ CC events).
The channel split differs from GENIE's: GiBUU QE 0.18 vs GENIE 0.43, RES 0.68 vs 1.11, DIS 3.08 vs
2.96, 2p2h 0.13 vs 0.12 (10⁻³⁸ cm²/nucleon at ⟨E⟩ ≈ 6 GeV; the GiBUU QE number carries the weight
problem of Part 2.4 and is statistically poor).

## 3.3 Into the platform

- One job: `read_job(job_dir, spec)` → `TruthTable` with `norm = {kind: xsec_per_nucleon,
  xsec_per_unit_weight: 1e-38}`, `has_geometry: False`, `frame: beam`, `has_fs: True` (so the 1μ1p
  signal definition applies unchanged).
- K jobs: `merge_jobs` prescales each job's weights by 1/K and concatenates
  (`TruthTable.concatenate` keeps the first table's normalisation, which is why the prescale is
  needed), records per-job σ (mean/std/min/max) and saves `truth.npz` + `gibuu_run.json` in
  `runs/_generator_cache/gibuu_<fingerprint>/`.
- Model spec `kind: gibuu` (`ndp.theory.models._gibuu`): the cache when it exists; `local: {num_ensembles: N}`
  runs one local job; otherwise it raises with the `ndp gibuu` recipe (a grid campaign is never launched from
  `ndp run`).
- Absolute normalisation downstream: `N_true = σ_cell × n_nucleons × phi_per_pot_cm2 × POT_data`
  (`ndp/compare/folded.py::expected_true_cells`); for the 1μ1p channel n_nucleons = 3.23e30 (default,
  2026-09-15) and Φ = 6.32e-8 cm⁻² per POT (0–100 GeV, the range the flux file covers).

## Gotchas

- `FinalEvents.dat` is overwritten by every run in the same directory; `num_runs_SameEnergy > 1`
  appends runs with increasing run numbers (column 1), and the adapter divides the weights by the number
  of runs (`n_runs`, taken from the maximum run number unless given).
- The `.gz` produced by the worker is decompressed by `ndp.theory.gibuu._finalevents_path` on first use.
- Unit-weight samples (`equalWeights_Mode = 2`) still carry a common non-unit `perweight`
  (`equalWeights_Max/(A numEnsembles)`); nothing downstream assumes weights of 1.
