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
- **Streaming-only policy (2026-06-12): AnaTuples are never copied to local disk.**
  All reads go through xrootd streaming; uproot fetches only the requested
  branches/trees (verified offsite: Meta-tree POT reads of both golden files).
  File identity is checked by streamed fingerprints — POT_Used/POT_Total + tree
  entry counts vs the frozen 2026-06-09 manifest — not sha256 (which would require
  reading whole files). The only local-fetch input in the entire pipeline is the
  FluxAndReweightFiles tarball (a .tgz cannot be streamed), deferred to the weights
  stage that consumes it.

## Step 0 — Scaffold + git init + commit plan  ✅ (this commit)

README.md (4-stage workflow), .gitignore, this plan, `config/`, `xsec/__init__.py`,
`tests/`, `results/`. `git init` branch main. *Gate: clean status, plan committed.*

## Step 1 — Environment additions (runlog_tools, pdg)

1. Clone runlog_tools to `~/runlog_tools`; `pixi run pip install -e ~/runlog_tools`.
2. Add `pdg` to parent `pixi.toml` pypi-dependencies (locked) + `pixi install`.
3. Smoke test: `from runlog_tools import RunLog`; `pdg.connect()` returns muon
   mass/charge and ν_μ/μ⁻ MC IDs offline.
*Gate: both imports + printed PDG lookup (m_μ, IDs 13/14 fetched via API).* STOP.

## Step 2 — Stage 1: input-file structure (streaming-only) + input summary

Every input declared in a spec; verification and summary happen **by streaming** —
no AnaTuple is ever downloaded. (Already done ahead of this step, commit bdb6066:
official playlist file lists ingested to `config/playlists/` with 12 pinned input
tests — specs derive per-file xrootd URLs from those lists, never hand-written.)

1. `config/datasets/me1A_single_pair.json` — data run 10066 + MC run 110040:
   role, playlist (minervame1A), run, xrootd URL (taken from
   `config/playlists/*.txt`), and **expected streamed fingerprints** — Meta
   POT_Used/POT_Total and reco/Truth entry counts — referenced against the frozen
   exploration-repo manifest (`runs/2026-06-09_dsigma_dpt/manifest.json`:
   POT_Used 2.049772e17 data / 9.988797e18 MC, 6304 reco data, 186205 reco +
   544600 Truth MC). Remote size from `xrdfs stat` is informational.
   Spec format scales to playlists 1B–1P without code changes.
2. `config/datasets/aux_flux_reweight.json` —
   `FluxAndReweightFiles_Tarred_Feb_20_2026_1145_FNALTime.tgz`
   (`root://fndca1.fnal.gov:1095//pnfs/fnal.gov/usr/minerva/persistent/OpenData/FluxAndReweightFiles/`).
   **The one exception to streaming-only** (a .tgz cannot be streamed): spec written
   now; the actual fetch+unpack happens at the weights stage that consumes it.
3. `config/published.json` — pointer to the vendored anc/ answer key (referenced,
   not duplicated).
4. **`config/branches.json` — curated analysis branch catalog** (single source of
   truth for "which branches this analysis uses and why"): branch name → tree
   (reco/Truth/Meta), dtype, units, one-line definition, and role(s) ∈
   {`measurement` (muon p_T / p_∥ ingredients: lepton momentum + angles),
   `reco_selection` (the 6 cuts' branches), `signal_definition`
   (mc_incoming, mc_current), `phase_space` (truth vertex + truth lepton
   kinematics), `pot` (Meta/POT_Used, POT_Total)}. Definitions sourced from the
   exploration repo `docs/minerva/branches.md` + vendored Tuple-Documentation —
   cited per entry. Step 3's `REQUIRED_BRANCHES` is **imported from this catalog**
   (docs and code cannot drift).
5. **`summarize_inputs.py` (argparse + RunLog) — the Stage-1 summary feature.**
   Given a dataset spec (+ branch catalog), streams via xrootd and produces
   `results/<ts>__summarize_inputs/{summary.md,summary.json}`:
   - per file: run, role, playlist, remote size (`xrdfs stat`), Meta
     POT_Used/POT_Total, reco-tree entries (+ Truth entries for MC), fingerprint
     status (streamed values vs spec expectations);
   - per data file: **data-taking time span** from the event GPS-time branches
     (exact branch names verified against the tuple at implementation; candidates
     `ev_gps_time_sec`/`ev_gps_time_usec`), cross-checked against the published
     run-period window for the playlist (ME1A: 12-Sep-2013 → 14-Jan-2014); MC
     marked "simulated" (timestamps not physical);
   - dataset level: total POT, file count, run range, date range;
   - **branch contract check**: every catalogued branch exists in the expected tree
     with the expected dtype (data-tree vs MC-tree presence semantics respected),
     reported per role group — the input is "valid for this analysis" iff the
     contract passes.

## Step 3 — Stage 2: cuts code (reco selection + signal definition)

Vectorized uproot/awkward, arrays in / boolean masks out, no I/O inside:
1. `xsec/kinematics.py` — θ3d(θx,θy); p_T, p_∥ incl. NuMI beam rotation; m_μ and
   PDG IDs from the `pdg` API (no literals).
2. `xsec/cuts.py` — 6 named mask functions + `reco_selection(arrays)`:
   ZRange [5980,8422] mm, Apothem 850 mm, θ_μ<20°, isMinosMatchTrack==1,
   deadtime≤1, minos_trk_qp<0 (MINOS-match guard ordering preserved);
   `REQUIRED_BRANCHES` imported from `config/branches.json` (role
   `reco_selection`). Constants from `config/constants.py`.
3. `xsec/signal.py` — truth signal (mc_incoming==ν_μ && mc_current==CC) and
   phase space (z, apothem, θ≤20°, p_z≥1.5 GeV) for the efficiency denominator.
4. `tests/` — synthetic per-cut boundary tests + golden parity gate **streaming the
   golden pair via xrootd** (network-marked, skips cleanly offline; reads only the
   catalogued branches ≈ 1–2% of file bytes): **exactly 844 selected data / 43643
   selected MC (43539 signal, 104 background)**.
*Gates: pytest green; parity exact; branch list documented.* STOP.

## Later stages (outline only)

- Stage 3: `make_ingredients.py` — event loop **streaming source tuples via xrootd**
  → `.npz` with the six ingredients (the intermediates are the only local
  artifacts), binning-agnostic flat slots (1D pT and 2D pT×p_∥ as config).
- Stage 4: `extract_xsec.py` — bkg subtraction → NumPy D'Agostini (10 iter) →
  efficiency → Φ, T (TargetUtils port, ±2% gate vs 3.23e30), POT, widths → dσ;
  `compare_published.py` vs anc tables.
- Weights (MnvTunev1) + scale-up: per M2–M4 of the 2026-06-11 plan in the
  exploration repo; execution home is now this repo.

## Risks

1. Private runlog_tools needs SSH — verified working.
2. `pdg` offline DB behavior — smoke-tested in Step 1 first.
3. Streaming latency/repeat-read cost (no local copies by policy) — uproot reads
   only catalogued branches (~1–2% of bytes); slow full-loop passes are network-
   marked tests or RunLog-tracked runs; Stage-3 intermediates make downstream
   stages local.
4. Vectorized-cuts drift vs certified selector — blocking parity gate in Step 3.
