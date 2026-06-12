# The five CV reweights (MnvTune v1)

Every MC event in this analysis carries a central-value weight that is the
**product of five reweighters** — the MINERvA Tune v1 stack, instantiated in
`MINERvA-101-Cross-Section/runEventLoop.cpp:380-386` and multiplied by
`PlotUtils::Model::GetWeight` (`MAT-MINERvA/weighters/Model.h:53-61`). Data is
never weighted. Until this stack is implemented, all our MC products are
"unweighted CV" — the mapped ≈−9 % normalization + shape gap in
`docs/results/2026-06-12_playlist1A_2d_migration.md` is precisely what these
five weights correct.

Structural facts (source-confirmed):
- **Four weights are truth-only** (flux, GENIE non-res-π, 2p2h, RPA;
  `DependsReco() == false`) and apply to **every MC fill including the
  efficiency denominator** (Truth-tree loop). The **MINOS efficiency is the
  only reco-dependent one** and applies to reco-side fills only. `Model`
  caches the truth-only product per event (`fCVTruthOnlyWeight`).
- `LowQ2PiReweighter` is **not** part of Tune v1 (it belongs to later tunes) —
  the directory `MAT-MINERvA/weighters/` holds ~15 further reweighters used
  for alternative tunes, model curves, and systematics, none in the CV.

---

## 1. Flux CV weight — `FluxAndCVReweighter`

**Why.** The MC was generated with an old beam simulation (g4numi v6); the
best-estimate flux is the PPFX ("gen2thin") prediction further constrained by
MINERvA's ν-e elastic scattering measurement. The weight morphs the generated
flux into the constrained one, event by event in true E_ν.

**Formula** (`MAT/PlotUtils/FluxReweighter.cxx`, `GetFluxCVWeight`):
w = Φ_constrained.Interpolate(E_ν) / Φ_generated.Interpolate(E_ν), with plain
bin content above 75 GeV (slope guard) and w = 1 if either reads zero.
E_ν = `mc_incomingE`/1000 GeV.

**Inputs.** Branches: `mc_incomingE`, `mc_incoming`. Files (in our tarball,
`data/flux/`): `MATFluxAndReweightFiles/flux/flux-gen2thin-pdg14-<plist>_rearrangedUniverses.root`
(`flux_E_cvweighted`) and `flux-g4numiv6-pdg14-<plist>.root`
(`flux_E_unweighted`), with the playlist mapping me1A–1F → `minervame1D`,
1G/1L/1M → `minervame1M`, 1N/1O/1P → `minervame1N`; ν-e constraint weights:
`MParamFiles/data/FluxConstraints/sorted_NuEConstraint_FHC_RHC_IMD.txt`.

**Implementation note / open item.** The constrained CV is computed **at load
time** by `MnvHistoConstrainer::ConstrainHisto` (reweighting the file's Flux
error-band universes with the constraint weights) — that math must be
reproduced before any flux weight is trusted. Gate: Φ_int(0–100 GeV, width
integral) within ±5 % of the paper's 6.32×10⁻⁸ ν/cm²/POT.

**Expected size.** O(few %) normalization, E_ν-dependent — largest in the
falling flux edge above the focusing peak (the high-p_∥ data/MC rise we see).

## 2. Non-resonant pion reduction — `GENIEReweighter(true, false)`

**Why.** Reanalysis of deuterium bubble-chamber data (Rodrigues et al.) showed
GENIE's non-resonant single-pion rate is too high; Tune v1 scales those events
down to 43 %.

**Formula** (`GENIEReweighter.h:31-42`): w = 0.43 if the event is "non-res
single π", else 1. The tag (`GenieSystematics.cxx:368-373`):
`truth_genie_wgt_Rvn1pi[2] < 1 || truth_genie_wgt_Rvp1pi[2] < 1` — i.e. GENIE
itself marks sensitivity via its stored +1σ knob ratios. The constant:
`kNonResPiWeight = 0.43` (`GenieSystematics.h:14`). The second constructor
flag (deuterium MaRES pion tune) is **off** in Tune v1.

**Inputs.** Branches `truth_genie_wgt_Rvn1pi[2]`, `truth_genie_wgt_Rvp1pi[2]`
(present in both MC trees, verified). **No external file.**

**Expected size.** −57 % on a subset of RES/DIS-adjacent events; few-% net on
the inclusive sample.

## 3. Low-recoil 2p2h enhancement — `LowRecoil2p2hReweighter` (mode 0)

**Why.** MINERvA's own low-recoil data showed a large excess over
GENIE+Valencia 2p2h; Tune v1 adds an empirical 2D-Gaussian enhancement in
true (q0, q3) fitted to that data (≈+50 % 2p2h rate).

**Formula** (`MnvTuneSystematics.cxx:19-60`, `weight_2p2h.cxx:44-57`):
applies only to true 2p2h events (`mc_intType == 8`; CV mode skips QE) on
nuclei (`mc_targetZ ≥ 2` — no 2p2h on hydrogen);
w = 1 + N·exp(−½ z(q0,q3)/(1−ρ²)) with z the correlated 2D Gaussian argument;
fit parameters from
`MParamFiles/data/Reweight/fit-mec-2d-noScaleDown-penalty00300-best-fit` (in
our tarball; mode 0 = CV fit, modes 1/2/3 are nn+pp / np / QE variations used
as systematics).

**Inputs.** True (q0, q3) in GeV (`Getq0True`/`Getq3True` — from truth
kinematics branches), `mc_intType`, `mc_targetZ`, `mc_targetNucleon`.

**Expected size.** Large on the 2p2h channel (~×1.5), ~+2–4 % on the
inclusive sample, concentrated at low recoil / low p_T.

## 4. MINOS matching efficiency — `MINOSEfficiencyReweighter`

**Why.** The MINOS-match efficiency in data degrades with beam intensity
(dead time in MINOS) in a way the simulation does not model; a measured
data/MC correction vs (p_μ^MINOS, instantaneous intensity) fixes the
normalization of the *reco-selected* MC.

**Formula** (`MuonFunctions.h:163-168`): correction looked up from
`MinosMuonEfficiencyCorrection::Get(isFHC).GetCorrection(p_μ^MINOS [GeV],
batchPOT, isFHC)` — a measured table vs muon momentum × batch POT.
**Batch POT** (`MinervaUniverse.cxx:441-484`, resolves our open question):
`batch_pot = numi_pot / k` with k from the spill's batch structure —
structure 0 → k=6; 1 → k=4 (vertex batch < 3) else 8; 2 → k=5 (< 5) else 10;
3 or −1 → k=6.

**Inputs.** RECO-side branches, all four **verified present in both open-data
tuples** (2026-06-12): `MasterAnaDev_minos_trk_p` (via GetPmuMinos),
`numi_pot`, `batch_structure`, `reco_vertex_batch`. Efficiency table file:
under `MParamFiles/` in the tarball (exact path at the inventory step).

**Application.** The only weight NOT applied to the efficiency denominator
(it corrects reconstruction, which the denominator doesn't have).

**Expected size.** ~1–2 %, intensity-dependent.

## 5. RPA suppression — `RPAReweighter`

**Why.** The Valencia random-phase-approximation calculation screens the
weak response of the nucleus for QE at low Q² — GENIE's bare RFG lacks it.
Tune v1 multiplies true-QE events by the Valencia/GENIE ratio in (q0, q3).

**Formula** (`RPAReweighter.h` → `weightRPA`): w = ratio histogram lookup at
true (q0, q3) GeV, CV variation, for true QE on nuclei (target Z from truth;
no RPA on hydrogen); low-Q² region handled by the weight class's dedicated
parameterization. Ratio histograms from the RPA file under
`MParamFiles/data/Reweight/` (exact name at the inventory step).

**Inputs.** True (q0, q3), `mc_intType` (QE), `mc_targetZ`, ν PDG.

**Expected size.** Strong suppression (tens of %) of QE at low Q² — the main
driver of the low-p_T data/MC overshoot we observe.

---

## Status and plan

| weight | formula source read | inputs in tuples | external file located | implemented |
|---|---|---|---|---|
| Flux CV | ✓ | ✓ | ✓ (constraint math pending) | ✗ |
| Non-res π | ✓ | ✓ (36 genie_wgt families verified) | n/a | ✗ |
| 2p2h | ✓ | q0/q3 ingredients ✓ | ✓ (fit-mec file) | ✗ |
| MINOS eff | ✓ (incl. batchPOT) | ✓ (all four branches verified) | table file at inventory | ✗ |
| RPA | interface ✓ (class internals at implementation) | ✓ | file name at inventory | ✗ |

Implementation order (weights stage): unpack tarball → inventory →
`xsec/flux.py` (constrainer + CV weight + Φ_int gate) → `xsec/weights.py`
(non-res π, 2p2h, RPA, MINOS eff) → per-weighter hand-trace validation vs the
C++ formulas → weighted re-stream of playlist 1A → gates: normalization gap
moves from ≈−9 % toward |<3 %|, p_∥ tail flattens, low-p_T overshoot shrinks.
