# CLAUDE.md — NDP platform

You are working in the Neutrino Discovery Platform: a theorist's model goes in, a comparison with
neutrino-scattering data comes out, with every step recorded. The physicist owns every physics
choice (signal definition, phase space, binning, normalisation constants, which surrogate, whether
a result is believable); you run the machinery and report faithfully.

## Operating rules

- **Physics choices live in manifests, not code.** Channel YAMLs (`channels/`) and model YAMLs
  (`models/`) carry every value that is a physics decision, each with a `status` (decided /
  default / open). Never change a `decided` value; propose changes to `default` values with
  evidence; never fill in an `open` value silently — ask, or record it in `docs/open_questions.md`.
- **A number without a run directory does not exist.** Every comparison goes through
  `ndp.pipeline.run_model` so it lands in `runs/<id>/` with `manifest.json`, `scorecard.json`,
  `report.md` and figures. Quote from `scorecard.json`, not from memory.
- **Report both comparison modes side by side** (unfolded: published d²σ + covariance; folded:
  surrogate-smeared prediction vs reconstructed data). Do not combine χ² values into a verdict;
  read the shape χ² together with its α, and the folded −2lnL together with the data/pred ratio.
- **Surrogates are learned from paired MC and certified by closure.** A rebuilt surrogate must
  fold the training MC's truth back onto its own reco counts exactly
  (`tests/test_minerva_certification.py`). Say which surrogate a run used.
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
surrogates/     trained detector surrogates  resources/   flux tables etc.
runs/           run directories (see runs/README.md)
docs/           architecture, decisions, open_questions, roadmap
tests/          pytest-style tests + run_tests.py fallback runner
.claude/skills/ndp-model   the agent workflow for "test my model"
```

## Everyday commands

```bash
python -m ndp channels                       # what can be tested
python -m ndp models                         # example model specs (validated)
python -m ndp run models/<m>.yaml --channel minerva_me_cc_inclusive_ptpz
python -m ndp surrogate build --channel minerva_me_cc_inclusive_ptpz --source mc --kind all
python -m ndp data status                    # are the AnaTuples / caches present
python tests/run_tests.py                    # or, inside pixi: pixi run test
pixi run build-genie && pixi run snapshot-genie-env   # (re)build the in-repo GENIE
```

Site paths (data, MINERvA repo, GENIE environment, splines) are in `ndp.yaml`; see
`ndp/config.py` for the environment-variable overrides. `docs/architecture.md` explains the data
contracts; `docs/roadmap.md` says what is done and what is next.
