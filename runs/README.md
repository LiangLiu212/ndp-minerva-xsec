# runs/

One directory per model test, created by `python -m ndp run <model.yaml> --channel <channel>`:
`<YYYY-MM-DD>_<model>__<channel>[_N]/` with `manifest.json` (inputs, fingerprints, versions, git
state, timings), `model.yaml`, `channel.json`, `scorecard.json`, `report.md`, `figs/`.

`_generator_cache/` holds GENIE productions keyed by the spec fingerprint (reused across runs) and
`_generator_cache/spline_cache/` the parsed total cross sections. Run directories are not committed;
the manifest is what makes a run reproducible.
