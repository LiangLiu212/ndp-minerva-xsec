# Open questions

Unresolved items the pipeline must settle before (or while) the relevant stage is
implemented. One entry per question; move to **Resolved** with the answer + evidence
when closed.

## Open

### 2026-06-12 — MINOS-efficiency weighter: which beam-intensity input, exactly?

Detail worth remembering for the weights stage: the per-batch `POT_Used_batchN`
leaves in the `Meta` tree are **not** what the MINOS-efficiency reweighter uses —
that one reads per-event beam-intensity branches from the reco tree. The `Meta`
tree is purely exposure bookkeeping.

To settle at Stage-2.3 (weighters): identify the exact reco-tree branch(es) behind
`CVUniverse::GetBatchPOT()` (MAT-MINERvA `MINOSEfficiencyReweighter.h` →
`MinosMuonEfficiencyCorrection::Get(...).GetCorrection(pmu, batchPOT, isFHC)`), and
confirm they are filled in the open-data MasterAnaDev tuples.

### 2026-06-12 — Flux usage 1 (CV): normalization integral Φ_int

Flux enters the cross-section denominator once, as a scalar:
σ_i = Σ_j U_ij (N_data,j − N_bkg,j) / (ε_i · **Φ_int** · T · POT_data · Δ_i), with
Φ_int = ∫₀¹⁰⁰ᴳᵉⱽ Φ_CV(E) dE of the ν-e-**constrained** reweighted flux
(`Integral(FindBin(0), FindBin(100), "width")`, broadcast into the
`reweightedflux_integrated` template; ExtractCrossSection.cpp:124-130 then applies
1/(T·POT), ×1e4 for m²→cm², ÷ bin width). Data is touched by flux ONLY here.
For me1A: file `MATFluxAndReweightFiles/flux/flux-gen2thin-pdg14-minervame1D_rearrangedUniverses.root`,
hist `flux_E_cvweighted` (me1A–1F all map to the minervame1D flux file,
FluxReweighter.cxx playlistString ~1953).

To settle at implementation:
- confirm the histogram's units are ν/m²/POT/GeV (from file metadata or by
  reproducing the paper integral);
- confirm the edge-bin convention (FindBin(0)→FindBin(100) integrates whole
  boundary bins, not partial);
- gate: our Φ_int must reproduce the paper's 6.32e-8 ν/cm²/POT (= 6.32e-4 ν/m²/POT)
  within ±5%.

### 2026-06-12 — Flux usage 2 (CV): per-event flux weight on all MC fills

w_flux(E_ν) = Φ_constrained(E_ν) / Φ_generated(E_ν), evaluated with
`TH1::Interpolate(Enu_GeV)` (Enu = mc_incomingE/1000); above 75 GeV plain bin
content (slope guard); weight = 1 if either histogram reads 0
(FluxReweighter::GetFluxCVWeight). Applied to ALL MC ingredients — background,
migration, efficiency numerator AND denominator; it does not cancel in
ε = num/denom (it reshapes the E_ν mix inside each (p_T, p_∥) bin).
Generated flux for me1A: `flux-g4numiv6-pdg14-minervame1D.root`, hist
`flux_E_unweighted`.

To settle at implementation:
- **the key unread piece**: `MnvHistoConstrainer::ConstrainHisto`
  (MAT/PlotUtils/MnvFluxConstraint.cxx) — the constrained CV is computed at load
  time as a weighted combination of the file's Flux-band universes with weights
  from `MParamFiles/data/FluxConstraints/sorted_NuEConstraint_FHC_RHC_IMD.txt`;
  reproduce that math in Python before any CV flux weight is trusted;
- confirm linear-interpolation edge behavior (first/last half-bins) matches
  TH1::Interpolate;
- sanity band: w_flux distribution over selected signal events should be smooth
  and O(1) (document observed range when first run).

### 2026-06-12 — Background (CV): definition, categories, subtraction

Source-confirmed semantics (runEventLoop.cpp:120-162, Cutter.cxx:133-151,
ExtractCrossSection.cpp:199-231, Variable.h:30-31):

- **Definition:** background = passes reco selection ∧ FAILS the truth signal
  definition (mc_incoming==14 && mc_current==1). **Truth phase space is NOT part
  of the signal/background split** — a true ν_μ-CC event outside the truth phase
  space that passes reco selection counts as signal (fills efficiency numerator +
  migration); the acceptance correction (efficiency denominator = signal-def ∧
  phase-space) maps data back to the fiducial. Two distinct predicates — keep them
  as two functions in `xsec/signal.py`, with tests.
- **Categories:** bkgd_ID 0 = "NC" (mc_current==2), bkgd_ID 1 = "Wrong_Sign"
  (everything else: ν̄_μ CC, ν_e CC, …). Histograms `<var>_background_<label>`,
  binned in RECO variables, full CV weight stack applied.
- **Subtraction:** first step of the extraction chain, before unfolding:
  N_sub,j = N_data,j − (POT_data/POT_mc)·Σ_cat N_bkg,cat,j. No sideband fit, no
  data-driven scale factor.
- **All subtracted background is MC prediction** (justified by 99.8% purity; paper:
  8655 / 4,105,696 = 0.2%). Nuances: (a) the MC is data-informed via MnvTunev1
  (2p2h fit, nonres-π reduction, ν-e flux constraint) but carries no
  background-specific data constraint; (b) sideband machinery exists
  (`isMCSelected(...).all()`) and is used by the CCQE-like companion analysis —
  not needed here; (c) non-beam backgrounds (rock muons, cosmics) are handled by
  the selection itself and never appear as an MC category.
- **Golden-pair reference:** 104 background of 43,643 selected MC (NC 47,
  Wrong_Sign 57); POT-scaled (×0.020521) ≈ 2.1 events vs 844 data.

To settle at implementation:
- confirm the mc_current convention (1=CC, 2=NC) holds in the open-data
  MasterAnaDev tuples (read values, don't assume);
- gate: reproduce the golden split exactly — 104 = 47 NC + 57 Wrong_Sign on the
  golden pair (unweighted CV);
- decide whether the Stage-3 ingredients file keeps the two categories separate
  (NC / Wrong_Sign) or stores one merged background histogram — separate is
  cheap and preserves the tutorial's breakdown plots.

## Resolved

(none yet)
