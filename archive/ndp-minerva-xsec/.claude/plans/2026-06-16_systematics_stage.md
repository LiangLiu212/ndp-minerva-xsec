# Plan: systematics stage (M4) — total uncertainty + covariance vs anc

**Status: plan only — stepwise execution, each step on explicit go-ahead, ends
with its gate. No implementation until triggered.**

## Context

E0–E4 are done: the stat-only chain reproduces the published 2D cross section
to 0.13 % integrated (playlist 1A). The paper's *total* uncertainty includes
systematics; this stage propagates them and validates against the paper's
ancillary covariance files.

A systematic = an **ensemble of universes**, each a full re-run of the chain
with one input varied; the covariance is the **sample covariance of the
per-universe cross sections** (MnvVertErrorBand::CalcCovMx):
`Cov[i,k] = (1/ΣW) Σ_u W_u (σ_u[i]-σ̄[i])(σ_u[k]-σ̄[k])`; ±1σ pairs use
`((σ⁺−σ⁻)/2)`. **Total = Σ_groups Cov_group + Cov_stat** (MnvH1D::GetTotalErrorMatrix).

Two universe kinds:
- **vertical** (weight-only): same selected events/bins, only the MC weight
  changes — flux, GENIE knobs, 2p2h, RPA, MINOS-eff, Geant-hadron (~165 univ);
- **lateral**: the reco observable shifts, events re-bin/re-select — muon
  energy scale, beam angle, muon resolution (~10 univ).

Data is never varied (only the MC-derived bkg/migration/efficiency move). All
inputs are confirmed present in the MC tuples (GENIE 28 knobs as 7-element
arrays; hadron renorm tables in the tarball; etc. — see the systematics
analysis + docs/cv_reweight.md). The 1000 flux universes are already loaded in
`xsec.flux.FluxCV`.

**Validation targets** (papers/.../2106.16210/anc/):
`cov_ptpl_*_{stat, flux, energyscale, total}.txt`.

**Key efficiency:** the ~165 vertical universes ride along in ONE streaming
pass via a weight *matrix* W[event, universe] — per file, fill the per-universe
ingredient histograms, discard W. The ~10 lateral universes need a re-binning
pass but are cheap. Per-universe extraction (E4 chain) is in-memory matrix ops.

Inventory & per-universe formulas: the exploration repo's
`docs/minerva/systematics.md` (the "Yes" rows) and the MAT-MINERvA `universes/`
classes are the authority for which universes, their counts, and how each reads
its inputs.

---

## Step S1 — framework: universes, weight-matrix pass, covariance assembler

- `xsec/systematics.py`: a `Universe` abstraction + a covariance assembler
  `group_covariance(sigma_cv, sigma_universes, kind)` implementing the
  CalcCovMx sample covariance (many-universe outer product; ±1σ-pair variant);
  `total_covariance(group_covs, cov_stat)`.
- Extend the weighters/flux with a `universe_weights(arrays, names)` interface
  returning per-universe weight columns (CV × universe ratio); GENIE knobs read
  `truth_genie_wgt_<knob>` at idx 2/4 vs CV idx 3; 2p2h/RPA via their variation
  modes; flux via the 100 PPFX universe histograms.
- `make_ingredients.py --universes <spec>`: one streaming pass producing the
  ingredient stack `(n_univ, 288)` per ingredient (data_reco shared). Memory:
  ~110 MB for 165 vertical universes; per-file W discarded.
- A per-universe extraction runner (reuse `extract_xsec.extract` per column).
*Gates: the "cv" universe reproduces the E4 CV cross section bit-for-bit;
a 2-fake-universe toy gives the analytic covariance; covariance is symmetric
PSD; pytest green.*

## Step S2 — flux systematic (100 universes) → validate cov_flux

Vertical, biggest, and the inputs are in hand. Per-universe flux weight (each
PPFX universe) AND per-universe integrated flux Φ_u (each universe's histogram).
Run the chain per universe → Cov_flux (224×224).
*Gate: diagonal √Cov_flux per cell vs `cov_ptpl_*_flux.txt` — fractional flux
error in the well-populated cells within ~20 % (single-playlist + ν-e
constraint subtleties expected; document the constraint handling).*

## Step S3 — statistical covariance → validate cov_stat

Data-Poisson toys: throw N≈1000 Poisson replicas of `data_reco`, run the full
chain per toy, covariance from the spread (replaces RooUnfold's analytic
propagation, including the unfolding's iteration feedback that E4's approximate
`var` omits).
*Gate: diagonal √Cov_stat vs the E4 per-cell `dsigma_err` (same order) AND vs
`cov_ptpl_*_stat.txt` scaled by √(POT ratio) — full-dataset stat is ~√11.8×
smaller than 1A, so compare shapes/correlation structure, not absolute size.*

## Step S4 — muon energy scale (lateral) → validate cov_energyscale

±1σ shift of reco muon momentum (MuonUniverseMinerva + MuonUniverseMinos; the
MINOS one also re-evaluates the flux weight) → recompute reco p_T/p_∥, re-select,
re-bin → two shifted ingredient sets → σ⁺, σ⁻ → Cov_energyscale = outer product
of `(σ⁺−σ⁻)/2`. Shift magnitudes from NSFDefaults / the muon-systematics classes.
*Gate: diagonal vs `cov_ptpl_*_energyscale.txt` — the paper's muon-energy-scale
band (≈1 %, rising at high p_T) reproduced in shape.*

## Step S5 — remaining systematics + total → validate cov_total

Vertical (same weight-matrix pass as S2): GENIE 26-knob ×±1σ (52), 2p2h (3),
RPA HighQ²/LowQ² (4), MINOS-eff (2), Geant-hadron p/n (4). Lateral: beam angle
x/y (4), muon resolution (2). Each → its group covariance. Then
`Cov_total = Σ all groups + Cov_stat`.
*Gates: `Cov_total` vs `cov_ptpl_*_total.txt` — diagonal (per-cell total
fractional uncertainty) median ratio within ~30 % over reported cells; the
final d²σ with total error bars vs published `data_result` total column. (1A
stat dominates at single-playlist scale, so the total comparison is best done
after the full-dataset combine — note as a dependency.)*

---

## Records (docs/decisions.md)
universe inventory + counts actually enabled; covariance convention (sample cov
vs ±1σ pair); flux ν-e-constraint handling in universe space; lateral shift
magnitudes + provenance; toy count for stat.

## Risks
1. ν-e flux constraint in universe space — the constrained universes are a
   reweighting of the raw PPFX set; reproducing the constrained covariance (not
   just raw PPFX spread) is the subtle part → validate vs cov_flux early (S2).
2. Streaming cost / door stalls — the vertical mega-pass is one stream; add the
   xrootd timeout+retry hardening (flagged in E1) before launching it.
3. Lateral re-selection — events crossing the selection boundary under the shift
   must be re-cut, not just re-binned → covered by the per-universe full chain.
4. Single-playlist stat dominance — at 1A scale the stat band dwarfs most
   systematics; clean systematic validation (esp. cov_total) wants the
   full-dataset combine first → sequence accordingly or validate per-group
   shapes at 1A.
5. anc covariance ordering/units — GlobalID order + (cm²/(GeV/c)²)²,
   width-divided; cross-check against the data_result diagonal (= stat+sys).

## Verification
Each step's group covariance validated against its dedicated anc file
(stat/flux/energyscale); the assembled `Cov_total` against `cov_total` and the
final width-divided d²σ-with-total-errors against the published data_result.
Internal gates: cv-universe reproduces E4 exactly; covariances symmetric PSD;
toy stat matches analytic order.
