# Plan: pipeline scaffold + Stage 1 (inputs) + Stage 2 (cuts)

**Status: approved 2026-06-12. Execution stepwise — each step runs only on explicit
user go-ahead and ends with its gate evaluated. Steps 0 done at commit time of this file.**

## Context

This repo is the cross-section pipeline; the exploration repo
`../ndp-minerva-data-release-exploration` stays as reference (certified selector +
golden counts, frozen 1D result, paper answer key in `papers/minerva/2106.16210/anc/`).

Pipeline = MINERvA-101 ingredients taxonomy in four stages:
**1. input files → 2. cuts (reco selection + signal definition) → 3. raw ingredients
(data hist, background, migration, efficiency num/denom) → 4. unfolding +
normalization (flux, POT, targets T) → cross section.**

**Scope of THIS plan: Step 0 scaffold + Stage 1 (input-file structure, all external
data) + Stage 2 (cuts code). Stages 3–4 outlined only; separate plans later.**

Decisions captured:
- Own git repo (this one); GitHub remote added later by user.
- Run metadata via `runlog_tools` (git@github.com:LiangLiu212/runlog_tools.git,
  private): exports `RunLog`, `make_parser`, `args_to_inputs`, `add_label`,
  `add_outdir`, `default_outdir`; JSON logs to `~/log/<project>/<ts>.log`.
- **All physics constants from external databases, never hardcoded**: PDG IDs,
  masses, widths, charges via the official `pdg` PyPI package (v2026.0, local DB).
  MINERvA-specific constants (fiducial z, apothem, nPlanes…) in `config/constants.py`
  with provenance comments (NSFDefaults.h / TargetUtils / paper).
- Environment: parent pixi workspace `/home/feanor/ndp-genesis-agent/pixi.toml`,
  auto-discovered (this repo has no own pixi.toml).

## Step 0 — Scaffold + git init + commit plan  ✅ (this commit)

README.md (4-stage workflow), .gitignore, this plan, `config/`, `xsec/__init__.py`,
`tests/`, `results/`. `git init` branch main. *Gate: clean status, plan committed.*

## Step 1 — Environment additions (runlog_tools, pdg)

1. Clone runlog_tools to `~/runlog_tools`; `pixi run pip install -e ~/runlog_tools`.
2. Add `pdg` to parent `pixi.toml` pypi-dependencies (locked) + `pixi install`.
3. Smoke test: `from runlog_tools import RunLog`; `pdg.connect()` returns muon
   mass/charge and ν_μ/μ⁻ MC IDs offline.
*Gate: both imports + printed PDG lookup (m_μ, IDs 13/14 fetched via API).* STOP.

## Step 2 — Stage 1: input-file structure + all external data

Every input declared in a spec; one finalized script materializes and verifies.
1. `config/datasets/me1A_single_pair.json` — data run 10066 (196 MB) + MC run 110040
   (21.6 GB): role, playlist (minervame1A), run, xrootd URL
   (`root://fndcadoor.fnal.gov:1095//pnfs/fnal.gov/usr/minerva/persistent/OpenData/MediumEnergy_FHC/{Data,MC/StandardMC}/Playlist1A/…`),
   expected sha256 + size from the frozen exploration-repo manifest
   (`runs/2026-06-09_dsigma_dpt/manifest.json`), POT recorded after first read.
   Spec format scales to playlists 1B–1P without code changes.
2. `config/datasets/aux_flux_reweight.json` —
   `FluxAndReweightFiles_Tarred_Feb_20_2026_1145_FNALTime.tgz`
   (`root://fndca1.fnal.gov:1095//pnfs/fnal.gov/usr/minerva/persistent/OpenData/FluxAndReweightFiles/`),
   size via `xrdfs stat`; unpack target recorded.
3. `config/published.json` — pointer to the vendored anc/ answer key (referenced,
   not duplicated).
4. `fetch_data.py` (argparse + RunLog): `--spec <json> --data-root <dir>` (default
   `/home/feanor/ndp-genesis-agent/data/`, shared outside the repo); `xrdfs ls`
   check → resumable `xrdcp` of missing files → sha256 verify → tarball unpack →
   summary; RunLog records every file.
5. Run for both specs (21.6 GB MC in background).
*Gates: sha256 match frozen manifest; tarball unpacked + contents listed; re-run is
an idempotent no-op.* STOP.

## Step 3 — Stage 2: cuts code (reco selection + signal definition)

Vectorized uproot/awkward, arrays in / boolean masks out, no I/O inside:
1. `xsec/kinematics.py` — θ3d(θx,θy); p_T, p_∥ incl. NuMI beam rotation; m_μ and
   PDG IDs from the `pdg` API (no literals).
2. `xsec/cuts.py` — 6 named mask functions + `reco_selection(arrays)`:
   ZRange [5980,8422] mm, Apothem 850 mm, θ_μ<20°, isMinosMatchTrack==1,
   deadtime≤1, minos_trk_qp<0 (MINOS-match guard ordering preserved);
   exports `REQUIRED_BRANCHES`. Constants from `config/constants.py`.
3. `xsec/signal.py` — truth signal (mc_incoming==ν_μ && mc_current==CC) and
   phase space (z, apothem, θ≤20°, p_z≥1.5 GeV) for the efficiency denominator.
4. `tests/` — synthetic per-cut boundary tests + golden parity gate on staged
   files: **exactly 844 selected data / 43643 selected MC (43539 signal, 104
   background)**; skips cleanly when data not staged.
*Gates: pytest green; parity exact; branch list documented.* STOP.

## Later stages (outline only)

- Stage 3: `make_ingredients.py` — event loop → `.npz` with the six ingredients,
  binning-agnostic flat slots (1D pT and 2D pT×p_∥ as config).
- Stage 4: `extract_xsec.py` — bkg subtraction → NumPy D'Agostini (10 iter) →
  efficiency → Φ, T (TargetUtils port, ±2% gate vs 3.23e30), POT, widths → dσ;
  `compare_published.py` vs anc tables.
- Weights (MnvTunev1) + scale-up: per M2–M4 of the 2026-06-11 plan in the
  exploration repo; execution home is now this repo.

## Risks

1. Private runlog_tools needs SSH — verified working.
2. `pdg` offline DB behavior — smoke-tested in Step 1 first.
3. 21.6 GB xrootd egress slow — background + resumable.
4. Vectorized-cuts drift vs certified selector — blocking parity gate in Step 3.
