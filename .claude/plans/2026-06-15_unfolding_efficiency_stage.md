# Plan: unfolding + efficiency stage (→ absolute d²σ/dp_T dp_∥, stat-only)

**Status: plan only — stepwise execution, each step runs on explicit user
go-ahead and ends with its gate evaluated. No implementation until triggered.**

## Context

The ingredients/cuts/weights stages are done and validated: the certified
selection, the truth signal/phase-space predicates (`xsec/signal.py`), the full
MnvTune v1 CV weight stack (`xsec/flux.py` + `xsec/weights.py`), and the paper
grid (`xsec/binning.py`). The weighted full-1A run reproduced the paper's
data/MnvTune v1 = 1.118 to within 0.3 %.

This stage completes the extraction chain
**bkg-subtract → unfold → efficiency-correct → normalize → d²σ** and validates
the absolute 2D cross section against `papers/.../2106.16210/anc/`.
**Stat-only** (systematics are the separate M4 stage; this stage builds the
hooks but propagates only data Poisson errors).

What exists: `data2d`, `mc2d_bkg`, a 224×224 in-grid migration, Φ_int
(`FluxCV.integral`), POT. **Missing** (this plan): a count-conserving migration
+ the efficiency denominator (Truth-tree loop), the D'Agostini unfolder, the
nucleon count, and the extraction orchestrator.

Reference implementations to port (exploration repo, frozen):
`exploring/dsigma_dpt.py:334` `dagostini_unfold(data, var, migration, prior, n_iter)`
and `:138` `tracker_n_nucleons()` (nPlanes=108 → 3.235e30). Chain order +
1e4 m²→cm² convention confirmed from `ExtractCrossSection.cpp:124-130`.

**Key design fact (from the reference):** with under/overflow included so the
migration is count-conserving, `eff_num[true] = migration.sum(reco axis)` comes
for free — only the **denominator** needs new streaming. Efficiency =
num/denom then corrects acceptance, exactly as the tutorial does.

---

## Step E0 — binning: count-conserving slots

`xsec/binning.py`: add flat-slot indexing that includes under/overflow so no
selected-signal event is dropped. p_T: 14 bins + overflow (>4.5); p_∥: under
(<1.5) + 16 bins + overflow (>60). Total `N_SLOTS` (e.g. 15×18=270); the 224
measurement cells are the in-grid subset, with an explicit `MEAS_SLOTS` index
map. New `slot_ids(pt, pl)` (never returns −1) alongside the existing
`cell_ids`. Standardize migration orientation as **M[reco_slot, true_slot]** to
match the unfolder (our current `plot_2d_ptpl` used [true, reco] — transpose).
*Gate: round-trip the 224 measurement slots ↔ GlobalIDs; every event lands in
exactly one slot; pytest green.*

## Step E1 — make_ingredients.py: the six ingredients (reco + Truth loops)

Finalized script (argparse + RunLog), streaming-only, `--weights cv` default,
`--playlist`. One pass streams, per file:
- **reco tree** (data + MC): data reco hist (slots), MC background (slots,
  weighted), migration M[reco_slot, true_slot] for reco-selected **signal**
  (weighted, count-conserving) — its reco-axis column-sum is the eff numerator;
- **MC Truth tree** (the new piece): efficiency **denominator** =
  signal ∧ phase-space (`signal.is_efficiency_denominator`), true-binned in
  slots, with the **truth-only** CV weight (`weights.truth_cv_weight` — no MINOS
  eff on the denominator, matching Model semantics). Handle the live Truth cycle
  (assert 544,600 on the golden file).

Output `ingredients.npz`: `data_reco`, `bkg`, `migration`, `eff_denom`,
`pot_data`, `pot_mc`, `weight_mode`, slot edges. Plus `summary.json`
(efficiency map, counts).

*Gates:*
- reco counts reproduce the weighted-1A run (data 357,547 in grid; selected
  signal weighted ≈ 1.44 M);
- Truth loop: unweighted signal∧phase-space count reproduces the frozen
  **65,712** (exploration repo anchor) before weighting;
- **average efficiency** = Σ eff_num / Σ eff_denom ≈ **0.64** (paper's stated
  64 % over the phase space) — the key physics gate;
- count conservation: migration total = selected-signal total (0 dropped).

## Step E2 — xsec/unfold.py: D'Agostini + MC closure

Port `dagostini_unfold` verbatim (dimension-agnostic, on flat slot vectors);
response = column-normalized migration; prior = MC true projection
(migration column-sum) or `eff_denom`. `n_iter` default **10** (paper).
*Gates:*
- synthetic 2-bin closure;
- **MC self-closure**: fold the MC-true slot vector through the response and
  unfold → recover MC-true to < 1e-6 on populated slots (count-conserving
  migration makes this exact);
- full-chain closure deferred to E4.

## Step E3 — xsec/targets.py: nucleon count

Port `tracker_n_nucleons` (nPlanes=108, hex area, MC mass fractions, PDG/
external atomic data per the project's no-hardcode rule where possible).
*Gate: N within 2 % of the published 3.23e30 (expect +0.16 %).*

## Step E4 — extract_xsec.py + compare_published.py: absolute d²σ

Orchestrate: bkg-subtract (`data_reco − (POT_data/POT_mc)·bkg`) → D'Agostini
unfold (10 iter) → divide by efficiency (num/denom) → divide by Φ_int → scale
`1/(N_nucleons·POT_data)`, ×1e4 (m²→cm²) → divide by 2D bin width (Δp_T·Δp_∥).
Output the 224-cell d²σ + the 1D projections. `compare_published.py` loads the
anc `data_result` table + covariances.
*Gates:*
- **MC self-closure through the full chain**: feed POT-scaled MC reco as
  pseudo-data → recover the MC truth cross section to < 1e-3 per populated bin;
- 1D p_T projection vs the frozen exploration result (+1.7 % unweighted; with
  CV weights expect the ~+12 % data/tune offset to move the absolute scale —
  cross-check direction);
- absolute 2D vs anc `data_result`: median per-bin ratio over the 205 reported
  cells recorded; integrated σ ratio. No hard external gate at stat-only +
  no-systematics, but the weighted result should sit within a few % given the
  1.121 event-rate match already shown.

---

## Records (docs/decisions.md when first touched)
under/overflow slot policy + migration orientation; efficiency asymmetry
(numerator = reco-selected signal, no phase space; denominator = signal ∧
phase space — tutorial Cutter semantics); truth-only weight on the denominator;
D'Agostini 10 iterations; nucleon nPlanes=108.

## Risks
1. Migration orientation/normalization bug → caught by MC self-closure gate.
2. Efficiency double-counting weights (MINOS-eff must NOT be on the denominator)
   → the 64 % average-efficiency gate catches gross errors.
3. Truth-cycle double count → assert 544,600; reuse the verified live-cycle read.
4. anc units/width factor → the 1D-projection + integrated-σ cross-checks catch
   any 1e4 / width / POT slip immediately.
5. Stat-only result won't match anc total error — expected; full validation
   waits on the M4 systematics stage (inputs already confirmed present).

## Verification
Each step ends with its gate. Stage deliverable: `extract_xsec.py` produces the
absolute 2D d²σ from a streamed `ingredients.npz`, MC-closes through the full
chain, and `compare_published.py` tabulates it against the paper's 205 reported
cells.
