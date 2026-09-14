# runs/

One directory per model test, created by `python -m ndp run <model.yaml> --channel <channel>`:
`<YYYY-MM-DD>_<model>__<channel>[_N]/` with `manifest.json` (inputs, fingerprints, versions, git
state, timings), `model.yaml`, `channel.json`, `scorecard.json`, `report.md`, `figs/`.

`_generator_cache/` holds GENIE productions keyed by the spec fingerprint (reused across runs) and
`_generator_cache/spline_cache/` the parsed total cross sections. Run directories are not committed;
the manifest is what makes a run reproducible.

`<YYYY-MM-DD>_signal_<channel>[_N]/` are diagnostics runs from `python -m ndp signal --channel <c>`:
the channel's truth-level signal definition applied to the cached official MC (`manifest.json`,
`channel.json`, `cutflow.json`, `summary.json`, `report.md`, `figs/`). They certify a signal
definition; they compare nothing with data.

`<YYYY-MM-DD>_selection_<channel>[_N]/` are diagnostics runs from `python -m ndp selection --channel <c>`:
the channel's reconstruction-level selection applied to the cached data and official MC
(`cutflow.json` with purity/efficiency and the proton-score scan, `summary.json` with the MC
composition and per-grid counts, `report.md`, `figs/data_vs_mc_<measurement>.png`). POT-scaled,
unweighted CV MC; a data/MC shape check, not a model test.
