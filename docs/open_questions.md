# Open questions

Decisions the physicist owns and the platform has not settled. Format:
`- [ ] date | scope | question | where it came up`.

## Active

- [ ] 2026-09-04 | minerva channel | **Ratify the beam-frame truth default** (`phase_space.frame: beam`). Evidence in `docs/decisions.md`. Consequence for the exploration repo: runs 2026-06-09 and 2026-06-19 built migration/efficiency with detector-frame truth; re-running their unfolding in the beam frame would test whether the shape χ²/ndf ≈ 3 they report is partly this. | surrogate build, `tests/test_minerva_certification.py`
- [ ] 2026-09-04 | minerva channel | **Which Φ normalises a POT-based MC?** The reference MC normalised with the published ν-e-constrained Φ (6.32e-8) sits +6.0% above the published d²σ (unfolded) and +6.6% above the selected data (folded, data/pred 0.938) — the same −6% pull the exploration repo saw. The generated MC flux (unconstrained gen2thin/g4numiv6) is ~12% above the constrained one per arXiv:2110.13372 §2, and the CV MC carries no flux/tune weights. Decide: published Φ (current), the MC's own generated-flux integral, or apply flux CV weights. | `models/reference_mc.yaml` runs
- [ ] 2026-09-04 | minerva channel | **Feed-in treatment.** Selected signal whose truth is outside the phase space (1.2k of 43k selected MC events) is a POT-scaled additive background in the surrogate. For a model that changes the out-of-phase-space rate this term should scale with the model; a `has_geometry` sample could supply it directly. Acceptable as is? | `ndp/surrogate/binned.py`
- [ ] 2026-09-04 | minerva channel | **n_nucleons**: 3.23e30 (arXiv:2106.16210, adopted by the exploration repo) vs 3.115e30 (arXiv:2110.13372, −2% as-built mass, possibly different fiducial). Same tracker fiducial? | `channels/minerva_me_cc_inclusive_ptpz.yaml`
- [ ] 2026-09-04 | minerva channel | **Target composition for generator runs**: mass fractions from arXiv:2110.13372 §1 (C 88.51 %, H 8.18 %, O 2.5 %, Ti 0.47 %, Cl 0.2 %, Al 0.07 %, Si 0.07 %) are used as GENIE `-t` weights. Confirm this is the intended per-nucleon normalisation target (vs. pure CH). | `ndp/theory/generator.py`
- [ ] 2026-09-04 | surrogate | **Systematic uncertainties** are not in the surrogate (only MC statistics). Flux universes, GENIE knobs, muon-energy scale and MINOS efficiency would enter as alternative response matrices or as covariance on the folded prediction. Which set, and from which weights (the tuple carries `truth_genie_wgt_*` and `mc_wgt_*`)? | `ndp/compare/folded.py`
- [ ] 2026-09-04 | surrogate | **Tails.** The parametric surrogate carries Gaussian cores only; the 1/|q/p| momentum tails are in the binned response but not in the smearing model. Is a flow-type surrogate the next step, or is the binned response sufficient for the channels planned? | `ndp/surrogate/parametric.py`
- [ ] 2026-09-04 | low-recoil channel | `channels/minerva_me_lowrecoil_eavail_q3.yaml` is a draft: exact bin edges from the supplemental Table II, the reco E_avail estimator (exploration repo open question `[~]`), and a manifest kind for the 44-bin release are all needed before it can run. | `channels/`
- [ ] 2026-09-04 | data | The open-data MC generator tag (GENIE 2.12.6 + which tune) is not embedded in the tuple (exploration repo open question). Affects how the `reference_mc` model should be labelled. | `ndp/adapters/minerva_anatuple.py`

## Resolved

- [x] 2026-09-04 | genie | Are GENIE `-t` target-mix weights number or mass fractions? → mass fractions: `GMCJDriver` uses them as density-weighted path lengths and divides by A (`InteractionProbability`), consistent with the H₂O `[0.8888],[0.1111]` convention. Recorded in `docs/decisions.md`.
- [x] 2026-09-04 | flux | Can the MAT `MnvH1D` flux files be read here? → No (no MAT dictionaries; MakeProject needs a compiler the pixi ROOT lacks). The published 2110.13372 Table I is the flux source; integral 0.991 × the published Φ.
