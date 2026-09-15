---
name: gibuu
description: "GiBUU (Giessen Boltzmann-Uehling-Uhlenbeck transport model) in this project, built up part by part. Part 1 = installation: what is installed where, verifying it, downloading from hepforge (bot wall), building from source, standalone installs, optional extras, updating. Use when the user says 'install GiBUU', 'build GiBUU', 'is GiBUU installed', 'which GiBUU version', 'update GiBUU', 'download GiBUU', 'buuinput', 'GiBUU.x', or when a task needs a working GiBUU. Part 2 = running: job cards, the platform card template and model spec, flux files (nuExp 99), ensembles, weighted vs equal-weights events, grid campaigns (ndp gibuu). Part 3 = output: FinalEvents.dat, cross-section files, the adapter and merge, normalisation. Triggers - gibuu, GiBUU.x, buuinput, hepforge, release2025, gibuu install, gibuu build, gibuu version, gibuu update, gibuu job card, FinalEvents.dat, perweight, equalWeights, gibuu flux, ndp gibuu."
argument-hint: "[install|check|update] [--gibuu <dir>]"
---

# gibuu — GiBUU in the NDP platform

Skill root: `<ndp-platform>/.claude/skills/gibuu/` (also reachable as `ndp-dev/.claude/skills/gibuu`,
a symlink). Every fact in `references/` was verified on the FNAL EAF node and carries the date it
was checked; re-verify anything that looks stale before acting on it.

## Parts

| part | covers | file | status |
|---|---|---|---|
| 1 Installation | what is installed, verify, download, build, standalone install, extras, update | `references/install.md` | written 2026-09-13 |
| 2 Running | the platform's card template and model spec, the `nuExp = 99` flux-file contract, ensembles/timing/memory, weighted vs equal-weights events, grid campaigns | `references/running.md` | written 2026-09-15 |
| 3 Output | `FinalEvents.dat` columns and units, the cross-section and flux-check files, the adapter and merge, normalisation downstream | `references/output.md` | written 2026-09-15 |

A part marked *not written yet* does not exist. Do not improvise it from general GiBUU knowledge
as if it were verified here: say so, help from GiBUU's own material (`$GIBUU/testRun/jobCards/`,
`$GIBUU/namelists.pdf`, https://gibuu.hepforge.org) with that caveat, and offer to write the part.

## Parts 2 and 3 in short (details and every number: `references/running.md`, `references/output.md`)

- **Run a sample:** a model spec of `kind: gibuu` (`models/gibuu_2025_me_fhc_c12.yaml`) + the channel
  flux → `python -m ndp gibuu smoke <model> --channel <c> --ensembles 100` (local, ~30 s) or a grid
  campaign (`ndp gibuu plan / submit-cmd / status / harvest / merge`, `grid/README.md`). Cards come
  from `resources/gibuu/minerva_me_numu_CC.job.tmpl` (GiBUU's own MINERvA-ME card) via
  `ndp/theory/gibuu.py`; the flux file is the channel table rebinned to equidistant 0.5 GeV bins.
- **Measured (C12, ME FHC flux, EAF node):** 1000 ensembles → 5 759 events in 88 s, 0.45 GB;
  σ_CC = 4.1 × 10⁻³⁸ cm²/nucleon; `numEnsembles < 100` aborts.
- **Weights:** default events are weighted (perweight = σ/(A numEnsembles), summing to σ_CC); under
  the ME flux QE weights are so uneven that 328 signal events carry 21 effective ones — use
  `equalWeights_Mode = 2` with a ceiling from a pilot (GiBUU aborts above it).
- **Output → platform:** `FinalEvents.dat` (15 columns, GeV, neutrino along +z = beam frame) is read
  by `ndp/adapters/gibuu_finalevents.py`; `neutrino_absorption_cross_section_ALL.dat` column 2 must
  equal the summed weights (the merge enforces it).

## Part 1 in short (details and every number: `references/install.md`)

- **Installed:** Release 2025, patch 5 (April 24, 2026) at `external/gibuu/release2025/objects/GiBUU.x`,
  tables at `external/gibuu/buuinput`. `activate.sh` exports `GIBUU`, `GIBUU_INPUT` and puts
  `$GIBUU/objects` on PATH under `pixi run` / `pixi shell`.
- **Verify:** `bash .claude/skills/gibuu/scripts/gibuu_check.sh [--smoke]` — checks the tree, the
  binary's shared libraries, the tables, compares version and tarball sha1s with what hepforge
  serves now, and with `--smoke` runs the numu-C12 smoke job (~25 s, a few hundred events).
- **Build (platform):** `pixi run build-gibuu` (= `scripts/build_gibuu.sh`: fetch, untar,
  `make FORT=gfortran -j8`); it skips any step whose result already exists, so an update needs the
  old tree moved aside first. Smoke test for all generators: `pixi run test-generators`.
- **Download:** hepforge sits behind an Anubis bot wall. Plain `curl` (default User-Agent) against
  `https://gibuu.hepforge.org/downloader?f=<file>` works; browser-like agents and WebFetch do not.
  Always verify with `?f=sha1sum.txt` and read `?f=version.txt` for the patch level.
- **Standalone install elsewhere:** section 1.5 of the reference, with `xxx` for the paths the
  user must choose.
- **Optional extras (not installed):** in-medium tables, RootTuple (ROOT output), HEPMC3event
  (HepMC3 output). Section 1.6 gives GiBUU's own procedure; it has not been exercised here.

## Conventions for growing this skill

- One `references/<part>.md` per part; this file stays a router (parts table + short summaries).
  When a part is added: fill its row above, extend the `description` triggers, keep numbering.
- A reference states only measured or verified values, each with its date and source (command,
  file, URL). Unknowns that the user must supply are written `xxx`, never a plausible-looking
  value (global CLAUDE.md rule).
- Project mechanics live in the repo (`scripts/build_gibuu.sh`, `pixi.toml`, `activate.sh`,
  `resources/gibuu/`, `ndp/adapters/gibuu_finalevents.py`); references point at them instead of
  copying them, and get updated when they change.
- `scripts/` here are bash or stdlib-Python only, default to the platform install and accept
  `--gibuu` / `--input` overrides so they also serve a standalone install.
- Gotchas found while working go into the part they belong to, in its final "Gotchas" section.
