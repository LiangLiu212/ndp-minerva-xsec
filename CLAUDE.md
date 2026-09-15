# CLAUDE.md — NDP platform

You are working in the Neutrino Discovery Platform: a theorist's model goes in, a comparison with
neutrino-scattering data comes out, with every step recorded. The physicist owns every physics
choice (signal definition, phase space, binning, normalisation constants, which surrogate, whether
a result is believable); you run the machinery and report faithfully.

## Operating rules

- **Physics choices live in manifests, not code.** Channel YAMLs (`channels/`), measurement YAMLs
  (`measurements/<channel>/`) and model YAMLs (`models/`) carry every value that is a physics
  decision, each with a `status` (decided / default / open; measurements: decided / user / example).
  Never change a `decided` value; propose changes to `default` values with evidence; never fill in
  an `open` value silently — ask, or record it in `docs/open_questions.md`. Example measurements
  are illustrations — an analyst who adopts one owns its edges.
- **Forward folding is the workflow.** Data = what was recorded; official MC = the detector's
  performance, used only to learn the surrogate; the theorist's model = what is tested, pushed
  forward through the surrogate into reco space. The platform never unfolds data.
- **A number without a run directory does not exist.** Every comparison goes through
  `ndp.pipeline.run_model` so it lands in `runs/<id>/` with `manifest.json`, `scorecard.json`,
  `report.md` and figures. Quote from `scorecard.json`, not from memory.
- **Report both comparison modes side by side where both exist** (folded: surrogate-smeared
  prediction vs reconstructed data — always; unfolded: published d²σ + covariance — only on the
  channel's `published` measurement). Do not combine χ² values into a verdict; read the folded
  −2lnL together with the data/pred ratio, and the shape χ² together with its α.
- **Surrogates are learned from paired MC and certified by closure.** One per (channel,
  measurement); a rebuilt surrogate must fold the training MC's truth back onto its own reco counts
  exactly (`ndp surrogate build` prints the closure; `tests/test_measurements.py`,
  `tests/test_minerva_certification.py`). Say which surrogate a run used.
- **The environment is `pixi.toml`.** `pixi install` builds it; `pixi run build-pythia6` /
  `build-genie` / `snapshot-genie-env` produce the in-repo GENIE under `external/`. Do not
  `pip install` into other environments or edit `external/genie/Generator` sources by hand —
  patches belong in `scripts/build_genie.sh` so a rebuild reproduces them. Outside pixi the default
  Python has no pytest; use `python tests/run_tests.py`.
- **Upstream stays upstream.** The MINERvA exploration repo (`ndp.yaml: minerva_repo`) is read,
  imported and cited, never edited from here. Findings that concern it (e.g. the truth-frame
  finding in `docs/decisions.md`) are reported to the user, who owns that repo.

## Layout

```
ndp/            package (theory/ adapters/ channels/ surrogate/ compare/ pipeline.py cli.py)
channels/       channel manifests            models/      example model specs
measurements/   observable pairs + binnings per channel (the published grid is implicit)
surrogates/     trained detector surrogates  resources/   flux tables etc.
runs/           run directories (see runs/README.md)
grid/           grid worker, payload staging, campaigns (see grid/README.md)
docs/           architecture, decisions, open_questions, roadmap
tests/          pytest-style tests + run_tests.py fallback runner
.claude/skills/ndp-model   the agent workflow for "test my model"
```

## Everyday commands

```bash
python -m ndp channels                       # what can be tested (channels + their measurements)
python -m ndp measurements --channel minerva_me_cc_inclusive_ptpz   # published grid + user observables, surrogate status
python -m ndp models                         # example model specs (validated)
python -m ndp run models/<m>.yaml --channel minerva_me_cc_inclusive_ptpz [--measurement <name>]
python -m ndp surrogate build --channel minerva_me_cc_inclusive_ptpz [--measurement <name>] --kind all
python -m ndp data status                    # are the AnaTuples / caches present (and current)
python -m ndp data cache --channel minerva_me_cc_inclusive_ptpz     # (re)build the truth/reco caches
python -m ndp signal --channel minerva_me_ccqelike_1mu1p           # truth-level signal definition on the cached MC -> diagnostics run
python -m ndp selection --channel minerva_me_ccqelike_1mu1p        # reco selection on cached data + MC: cutflow, purity/efficiency, data-vs-MC figures
python -m ndp data cache --channel minerva_me_ccqelike_1mu1p --url root://... --kind mc --out DIR   # one streamed file, grid-style outputs
python -m ndp grid plan <campaign> --channel <c> --beams FHC   # worklists; then grid/README.md: publish, submit, status, resubmit, harvest
python -m ndp data merge --beam FHC --playlist 1A               # per-file products -> playlist products (+ pot_1A.json)
python -m ndp surrogate build --channel minerva_me_ccqelike_1mu1p --measurement all --kind binned   # every grid's binned response in one pass over the playlists
python -m ndp efficiency run --channel minerva_me_ccqelike_1mu1p    # efficiency maps + background by category + factorised-ansatz closure -> runs/<date>_efficiency_<channel>/
python -m ndp efficiency apply --channel <c> --run runs/<eff run> --sample truth.npz --out weights.npz   # weight any truth sample with the maps
python -m ndp run models/gibuu_2025_me_fhc_c12.yaml --channel minerva_me_ccqelike_1mu1p --measurement all --modes folded --efficiency-run runs/<eff run>   # data vs model signal + MC background on every grid
python -m ndp gibuu smoke models/gibuu_2025_me_fhc_c12.yaml --channel minerva_me_ccqelike_1mu1p --ensembles 100   # one local GiBUU job; campaigns: grid/README.md (ndp gibuu plan/submit-cmd/status/harvest/merge)
python tests/run_tests.py                    # or, inside pixi: pixi run test
pixi run build-genie && pixi run snapshot-genie-env   # (re)build the in-repo GENIE
```

Site paths (data, MINERvA repo, GENIE environment, splines) are in `ndp.yaml`; see
`ndp/config.py` for the environment-variable overrides. `docs/architecture.md` explains the data
contracts; `docs/roadmap.md` says what is done and what is next.
