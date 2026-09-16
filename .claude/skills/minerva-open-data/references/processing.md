# Part 3 — Processing MINERvA data and MC in this workspace

Verified 2026-09-13 (survey of `/exp/dune/data/users/liangliu/minerva-open-data/`,
`ndp-minerva-data-release-exploration`, `ndp-platform`).

## 3.1 Two routes

| route | what | status here |
|---|---|---|
| **Official (C++)** | MAT + MAT-MINERvA + UnfoldUtils + MParamFiles + MATFluxAndReweightFiles; `MINERvA-101-Cross-Section`: `runEventLoop <data.txt> <mc.txt>` → `runEventLoop.root` (MnvH1D with all universes) → `ExtractCrossSection` | not built; the ROOT in the pixi env (6.40) satisfies the wiki's requirement, MAT itself would be a new build (skill part 4, unwritten) |
| **This workspace (Python)** | exploration repo (PyROOT tools, docs, 26 run directories incl. the audited 1D/2D cross sections) + ndp-platform (uproot adapter → npz caches → channel/measurement → surrogate → model comparison) | in use; certified against the audited numbers |

Consequences of the Python route: MAT `MnvH1D` flux/reweight files cannot be opened (no MAT
dictionaries), so the flux comes from the published arXiv:2110.13372 supplemental Table I
(`ndp-platform/resources/flux/arXiv2110.13372_supplemental.txt`; its integral 6.265e-8 = 0.991 × the
published Φ); flux/detector systematic universes are not yet propagated (the archived
`ndp-platform/archive/ndp-minerva-xsec/` code built flux/GENIE/detector universes and is the starting
point); the GENIE vertical universes *are* in the tuples (part 2).

## 3.2 Getting files

- Skill helper (stdlib + optional uproot):
  ```bash
  S=.claude/skills/minerva-open-data/scripts/minerva_od.py
  python3 $S filelists                              # scrape the site, download the 116 lists, summarise
  python3 $S filelists --playlist 1A --kind data    # print one list's URLs
  python3 $S fetch --playlist 1A --kind data --runs 6040 --dest xxx     # xrdcp with size check, skip existing
  python3 $S inspect /exp/dune/data/users/liangliu/minerva-open-data/*.root
  python3 $S status
  ```
- Exploration repo: `scripts/download_data.sh` (1+1 quick set = data 10066 + MC 110040; `--full` adds
  data 10062–10065 and MC 110036–110039; writes `shortData.txt`/`shortMC.txt`); files live outside
  the repo with a `data/ →` symlink and sha256 in `data/.hashes.json` (`.claude/data.md`).
- Budget: ~21.5 GB per StandardMC file, 41 of them for FHC 1A alone; `/exp/dune/data` is the place
  (6 PB filesystem), never the home directory. Transfers ran at ~85 MB/s on 2026-09-04.

## 3.3 What is on disk (2026-09-13)

`/exp/dune/data/users/liangliu/minerva-open-data/` (= `ndp.yaml: data_dir`, env `NDP_DATA_DIR`):

| file | size | content |
|---|---|---|
| `MasterAnaDev_data_AnaTuple_run00010066_Playlist.root` | 196,118,266 B | ME FHC 1A data, 6,304 reco entries, POT 2.04977e17 |
| `MasterAnaDev_mc_AnaTuple_run00110040_Playlist.root` | 21,594,091,620 B | ME FHC 1A StandardMC, 186,205 reco / 544,600 truth, POT 9.98880e18 |
| `cache/reco_data10066.npz`, `cache/reco_mc110040.npz`, `cache/reco_mc110040_truthcols.npz`, `cache/truth_mc110040.npz` | 0.95 MB, 27.8 MB, 14.4 MB, 197 MB | ndp-platform caches (`cache_version: 2`, rebuilt 2026-09-09); POT only on the truth cache (`__meta__.norm.pot`) |
| `xrdcp_mc_110040.log`, `cache/cache_build.log` | 0 B, 4 lines | stale leftovers, not manifests |

## 3.4 ndp-platform commands

```bash
cd /exp/dune/data/users/liangliu/ndp-dev/ndp-platform
python3 -m ndp data status [--channel minerva_me_cc_inclusive_ptpz]     # files present? caches fresh?
python3 -m ndp data cache --channel minerva_me_cc_inclusive_ptpz [--which data,mc] [--reco-only] [--entry-stop N]   # ~55 s
python3 -m ndp channels; python3 -m ndp measurements --channel minerva_me_cc_inclusive_ptpz
python3 -m ndp surrogate build --channel minerva_me_cc_inclusive_ptpz [--measurement <name>] --kind all
python3 -m ndp flux --channel minerva_me_cc_inclusive_ptpz [--out flux.root]
python3 -m ndp run models/<model>.yaml --channel minerva_me_cc_inclusive_ptpz [--measurement <name>]
python3 tests/run_tests.py            # includes tests/test_minerva_certification.py
```

- Adapter `ndp/adapters/minerva_anatuple.py`: uproot, vectorised; `RECO_BRANCHES`/`RECO_EXTRA_BRANCHES`
  (branch → cache column, e.g. `MasterAnaDev_recoil_E → reco_recoil_E`), `TRUTH_SCALARS/VECTORS/FS`,
  `read_pot` from `Meta`, cutflow with explicit NaN semantics, `MINERVA_INT_TYPE` 1/2/3/4/8 → NDP
  QE/RES/DIS/COH/MEC, `RECO_CACHE_VERSION = 2` (stale caches are refused with the rebuild command),
  `parity_vs_tool` re-runs the exploration repo's `tools/cc_inclusive_selector.py`.
- Reco cache columns: `passed, failing_cut, reco_p, reco_theta, reco_pT, reco_pz, reco_E_mu,
  reco_thetaX, reco_thetaY, reco_minos_qp, reco_minos_p, reco_vtx_{x,y,z}, reco_recoil_E,
  reco_recoil_E_passive, reco_recoil_E_calo, reco_recoil_E_polyline, reco_anatool_{E_nu,W,x,y},
  reco_visible_E, reco_n_prongs, reco_n_hadron_tracks`; truth cache: `nu_pdg, E_nu, lep_pdg,
  lep_{px,py,pz,E}, current, int_type, target_Z, target_A, Q2, W, weight, vtx_{x,y,z},
  generator_int_type` (+ `fs_*` final-state arrays).
- Channel `channels/minerva_me_cc_inclusive_ptpz.yaml`: `frame: beam`; fiducial z 5980–8422 mm,
  apothem 850 mm, θ ≤ 20°, p_z ≥ 1.5 GeV; pT edges `[0, 0.075, 0.15, 0.25, 0.325, 0.4, 0.475, 0.55,
  0.7, 0.85, 1.0, 1.25, 1.5, 2.5, 4.5]`, 16 p∥ bins 1.5–60 GeV; selection `minerva_cc_inclusive_v1`;
  POT "read from the Meta tree at run time, never typed in"; normalisation `phi_per_pot_cm2 = 6.32e-8`
  (±3.9 %), `n_nucleons = 3.23e30` (±1.4 %), tracker mass fractions C 0.8851 / H 0.0818 / O 0.0250 /
  Ti 0.0047 / Cl 0.0020 / Al 0.0007 / Si 0.0007 (also the GENIE `target_mix` of the platform's productions).
- Measurements (`measurements/<channel>/*.yaml`): user-defined truth↔reco observable pairs
  (`muon_p_theta`, `enu_calorimetric`, `q2_calorimetric` shipped as examples); reports are
  folded-first, unfolded only on the published grid.

## 3.5 Certified numbers (`tests/test_minerva_certification.py`, from the audited 2026-06-19 run)

| quantity | value |
|---|---|
| data selected (all 6,304 rows) | 844 (835 in-grid) |
| MC selected | 43,643 (43,361 in-grid; 43,266 in-grid signal; 43,175 reco&true in-grid) |
| truth signal in phase space | 65,041 (65,003 in-grid) |
| efficiency numerator | 41,948 (41,922 in-grid); pop-weighted ε = 0.64492 |
| purity | 0.9978; background raw 95 (43 NC + 52 other), POT-scaled 1.949 |
| Φ, N, POT_data | 6.32e-8 cm²/POT, 3.23e30 nucleons, 2.04977e17 |
| σ_tot, mean ratio to arXiv:2106.16210 | 3.1679e-38 cm²/nucleon, 1.0317 |

The certified counts are in the **detector** frame (the audited run's convention); the platform's
beam-frame default changes the migration, and the test pins both.

## 3.6 Known caveats of the audited runs (exploration `docs/decisions.md`)

Unweighted CV MC (≈ −6 % normalisation pull: no MnvTune/flux CV weights applied), published
ν-e-constrained Φ against unconstrained MC, NumPy D'Agostini with 10 iterations, detector frame with
the NuMI tilt neglected, negative-bin clamp, N_nucleons from TargetUtils with nPlanes = 108. The
exploration repo's `docs/notes/README.md` has no notes; `docs/open_questions.md` still lists the
recoil-energy choice and the generator tag.

## 3.7 Adding another playlist or sample

1. `python3 $S filelists --playlist <X> --kind <data|mc|Extended2p2h|…> --beam <FHC|RHC>` and
   `fetch … --dest /exp/dune/data/users/liangliu/minerva-open-data/<subdir>` (check free space first).
2. `python3 $S inspect <files>` — confirm the trees, POT and (MC) the `Truth` entry count.
3. Point a channel's `data:` section (or a new channel YAML) at the files; keep POT read from `Meta`.
4. `python3 -m ndp data cache --channel <channel>` then `python3 -m ndp data status`; rebuild
   surrogates. Never mix FHC with RHC, or LE with ME, in one sample.
5. Special MC samples need their POT-ratio weighting rules from part 1.4 (e.g. Extended 2p2h replaces
   the standard 2p2h events); nothing in the platform implements that yet.

## 3.8 Gotchas

- `Meta` POT must be summed over *all* files of a sample; the local pair's data/MC ratio (0.0205)
  differs from the 5+5 slice's (0.0498) and from full Playlist 1A (0.22027).
- The MC file is 21.5 GB: reading `Truth` takes ~11 s and reco ~26 s with the platform adapter;
  use `--entry-stop` while developing.
- uproot cannot stream xrootd in the default pixi python (no fsspec-xrootd); `xrdcp` first, or use the
  platform pixi env which has `xrootd` + `fsspec-xrootd`.
- MAT `MnvH1D` objects (flux files, `runEventLoop.root`) need the MAT libraries; plain ROOT/uproot see
  only the TH1 base part at best.
- Frame, Truth cycles, recoil families: part 2.
