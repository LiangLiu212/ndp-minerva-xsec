# Part 2 — The MasterAnaDev AnaTuple format

Verified 2026-09-13 by opening the two local files with uproot (`scripts/minerva_od.py inspect`) and
from the exploration repo's `docs/data_dictionary.md`, `docs/minerva/branches.md`,
`docs/minerva/analysis_recipe.md`, `docs/hadronic_system.md`, and ndp-platform's
`ndp/adapters/minerva_anatuple.py`. Official branch descriptions: `MAD_tuple_MainDoc.csv` in
github.com/MinervaExpt/Tuple-Documentation (257 rows, incomplete).

## 2.1 Trees

| tree | data file | MC file | content |
|---|---|---|---|
| `MasterAnaDev` | 6,304 entries, 3,687 branches (run 10066) | 186,205 entries, 4,190 branches (run 110040) | one entry per reconstructed event (pre-selected muon or electron candidate); MC entries also carry the matched truth (`mc_*`, `truth_*`) and the GENIE weight branches |
| `Truth` | — | 544,600 entries, 560 branches | every generated event (efficiency denominators); MC only |
| `Meta` | 1 entry | 1 entry | `POT_Used`, `POT_Total` (+ MC: `Total_Reco_Entries`, `Total_Truth_Entries`) |

**Truth key cycles:** the MC file holds `Truth;285` (live, 544,600) and `Truth;284` (stale, 543,576).
`TFile::Get("Truth")` and uproot's `f["Truth"]` return the highest cycle; anything iterating
`GetListOfKeys()` / `f.keys()` must dedupe by name.

POT: sum `Meta.POT_Used` over all files of a sample (`POT_Used ≤ POT_Total`; e.g. data run 6040 has
2.58403e16 used of 2.59901e16 total). Local pair: data `POT_Used = POT_Total = 2.04977e17`; MC `POT_Used = 9.98880e18`, `POT_Total = 9.98921e18`;
data/MC POT scale = 0.0205207 (the exploration repo's 5+5-file slice has 0.049825 — do not mix).
Unique event key: `(ev_run, ev_subrun, ev_gate, physEvtNum)`.

## 2.2 Units and key branches

All energies/momenta are **MeV**, angles radians, positions mm, `mc_Q2` in MeV².

| purpose | branch(es) |
|---|---|
| muon 4-vector | `MasterAnaDev_leptonE[4]` (px,py,pz,E in MeV); `MasterAnaDev_muon_E`, `MasterAnaDev_muon_theta` |
| muon angles | `muon_theta`, `muon_thetaX`, `muon_thetaY` — **beam frame** (see 2.5); `*_allNodes` variants |
| MINOS match | `isMinosMatchTrack`, `MasterAnaDev_minos_trk_p`, `MasterAnaDev_minos_trk_qp` (charge sign) |
| vertex | `vtx[4]`; dead-time `phys_n_dead_discr_pair_upstream_prim_track_proj` |
| recoil / hadronic | `MasterAnaDev_recoil_E`, `MasterAnaDev_hadron_recoil`, `_CCInc`, `_default`, `_two_track`, `MasterAnaDev_recoil_passivecorrected`, `blob_recoil_E`, `recoil_E_polylinecorrected`, `recoil_summed_energy[…]` (2.6) |
| calorimetric kinematics | `MasterAnaDev_E`, `_W`, `_x`, `_y` (`MasterAnaDev_Q2` is empty on the selected sample) |
| beam monitors | `numi_pot`, `numi_horn_curr`, `numi_bpm1..6`, `numi_tor101`, `numi_tortgt`, `numi_x/y(_width)` |
| truth (MC, both trees) | `mc_incoming` (PDG), `mc_current` (1 = CC), `mc_intType` (1 QE, 2 RES, 3 DIS, 4 COH, 8 MEC), `mc_primFSLepton[4]` (**detector frame**), `mc_vtx[4]`, `mc_Q2`, `mc_w`, `mc_targetZ/A`, `mc_FSPart*` (final-state particles); `truth_is_fiducial`, `truth_muon_E`, `truth_muon_theta` |
| GENIE weights (MC) | `genie_wgt_n_shifts`, `truth_genie_wgt_<knob>[n]` (vertical universes; knobs seen: `MaCCQE`, `MaCCQEshape`, `NormCCQE`, `MaRES`, `MvRES`, `NormCCRES`, `NormDISCC`, `Rvn1pi/2pi`, `Rvp1pi/2pi`, `AhtBY`, `BhtBY`, `CV1uBY`, `CV2uBY`, `AGKYxF1pi`, `Theta_Delta2Npi`, `RDecBR1gamma`, `VecFFCCQEshape`, `CCQEPauliSupViaKF`, `EtaNCEL`, `MaNCEL`, `NormNCRES`, FSI: `FrAbs_N/pi`, `FrCEx_N/pi`, `FrElas_N/pi`, `FrInel_N/pi`, `FrPiProd_N/pi`, `MFP_N/pi`); `truth_genie_wgt_shifts` |
| detector-systematic shifts | `MasterAnaDev_muon_theta_biasUp/Down`, `MasterAnaDev_proton_E_{BetheBloch,MEU}_biasUp/Down`, `_Birks_bias`, … |

Flux universes are **not** in the tuples (the 2026 update promises "increase in and fixes for flux
information saved in the tuples"); MAT's `FluxReweighter` takes them from `MATFluxAndReweightFiles`.

## 2.3 The CC-inclusive selection (MINERvA-101 recipe, as reproduced here)

Cuts in order (`docs/minerva/analysis_recipe.md:25-30`; constants in `ndp/adapters/minerva_anatuple.py:40-46`):

1. ZRange: 5980 ≤ vtx z ≤ 8422 mm (tracker fiducial)
2. Apothem: hexagonal apothem 850 mm (slope −1/√3)
3. MaxMuonAngle: θμ < 20°
4. HasMINOSMatch: `isMinosMatchTrack == 1`
5. NoDeadtime: `phys_n_dead_discr_pair_upstream_prim_track_proj ≤ 1`
6. IsNeutrino: `MasterAnaDev_minos_trk_qp < 0` (μ⁻)

Signal: `mc_incoming == 14 && mc_current == 1`; phase space / efficiency denominator: true vertex in
the same fiducial, θ_lep ≤ 20°, true p_z ≥ 1500 MeV/c. On the local pair this gives 844 selected data
events, 43,643 selected MC (43,266 in-grid signal), 65,041 truth signal-in-phase-space (`tests/test_minerva_certification.py`).

## 2.4 Data/MC normalisation and cross section

`dσ/dX = (N_data − (POT_data/POT_MC)·N_bkg) / (ε · Φ · N_nucleons · ΔX · POT_data)` with
`Scale(1e4)` for m² → cm² (`analysis_recipe.md:137-152`). MC histograms are scaled by the POT ratio;
the flux integral Φ comes from MAT's `flux_reweighter()` in the official chain (νμ, 0–100 GeV) and
from the arXiv:2110.13372 Table I flux here (part 3).

## 2.5 Frames (audited finding, 2026-09-04)

Reco `muon_thetaX/Y` (and MAT's `GetThetalepTrue`) are relative to the **NuMI beam axis**;
`mc_primFSLepton` is in the **detector frame**. Rotating the truth lepton about x by
`MinervaUnits::numi_beam_angle_rad = −0.05887` shrinks the reco−true pT core width from 0.278 to
0.061 GeV (same-cell migration fraction 0.10 → 0.40, pT bias −54 → +2 MeV) on MC run 110040. The
platform channel defaults to `frame: beam`; the exploration repo's audited runs used the detector frame.

## 2.6 Recoil energies (exploration `docs/hadronic_system.md`, marked "under review")

Seven branches, three families on the selected sample: ≈700 MeV `recoil_summed_energy[0]` (its array
axis is a 20×35 time-window scan, *not* universes); ≈844–885 MeV `blob_recoil_E ≡ recoil_E_nopolyline`
and `MasterAnaDev_recoil_passivecorrected ≡ MasterAnaDev_hadron_recoil_default`; ≈1281–1347 MeV
`MasterAnaDev_recoil_E ≡ _wide_window ≡ MasterAnaDev_hadron_recoil_CCInc` and
`recoil_E_polylinecorrected`. No canonical choice has been made; every calorimetric result carries this note.

## 2.7 Gotchas

- Everything is MeV; convert once (`MEV = 1e-3`, `MEV2 = 1e-6` in the platform adapter).
- `Truth` has two key cycles in the MC file (2.1).
- The generator/tune of StandardMC is not recorded in the file.
- Data files carry no `Truth` tree and no `mc_*` values; MC `MasterAnaDev` entries are reco-level
  and their `mc_*` branches are the matched truth, while `Truth` is the full generated sample.
- `MasterAnaDev_Q2` is unpopulated on the selected CC-inclusive sample; compute Q² from the lepton.
- Reading with PyROOT (exploration repo) and uproot (platform) agree row by row on the data file
  (`parity_vs_tool`, 6,304 rows, 0 mismatches).
