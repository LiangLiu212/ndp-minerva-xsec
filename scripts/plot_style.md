# scripts/plot_style.py

The single entry point for plot colors, markers, and style preferences.

## Purpose

Reads a user-owned preferences YAML (located as described below) and exposes
**role-based** accessors. This module is intended to be the **only** way a color
or marker enters a plot — plotting code asks `color_for("data")`, never writes a
hex string.

### Where preferences come from (no config required)

The module ships a **built-in default**, `scripts/style_preferences.default.yaml`,
so plotting works on a fresh checkout without first running any skill and
**without this repo's `docs/` layout**. Precedence:

1. an explicit `use_prefs_file(path)` (tests / alternate configs) — a missing
   path here raises rather than falling back, so a typo is not silently ignored;
2. `$NDP_STYLE_PREFS`, an absolute path to a prefs YAML — the hook for a user
   whose project is laid out differently from this repo;
3. `docs/style_preferences.yaml` if it exists (what the `/ndp-plot-prefs` and
   `/ndp-execute` skills materialize on first use);
4. the bundled `scripts/style_preferences.default.yaml`.

`use_prefs_file(None)` restores the normal order. `default_prefs_path()` reports
what 2–4 resolve to.

**What the bundled default actually is:** one project's seeded choices — the
`solid-dark` role palette, step histograms, magma sequential map — chosen to
reproduce this repo's pre-refactor figures. It is *not* the Okabe-Ito
colorblind-safe scheme; that ships alongside as the `okabe_ito` swatch table,
available to name a colorblind-safe color from an `override=`.

## Roles

Series **roles** drive every lookup, not call sites:

| Role | Used for |
|---|---|
| `data` | data points (black) |
| `mc_signal` | stacked MC signal; cutflow MC bar |
| `mc_bkg` | stacked MC background |
| `efficiency` | efficiency-vs-variable points |
| `purity` | purity bars |
| `ratio` | data/MC ratio points |
| `reference` | y=1 guide lines, ratio unity line + band |
| `syst_band` | systematic error band fill (same hue as `data`, shown via alpha) |

## API

| Call | Returns |
|---|---|
| `color_for(role, override=None)` | color string for a role |
| `marker_for(role, override=None)` | point marker (`data`/`efficiency`/`ratio`) |
| `palette(override=None)` | full `{role: color}` dict |
| `okabe_ito(override=None)` | named Okabe-Ito swatch table |
| `hist_style(override=None)` | `"bar"` \| `"step"` \| `"step_filled"` |
| `sequential_palette(override=None)` | colormap name for 2D plots |
| `style(key, override=None)` | a structural knob (markersize, capsize, linestyle, band bounds) |
| `response_prefs()` | the `response` block |
| `apply_rcparams()` | apply the `rcparams` block to matplotlib |
| `draw_hist(ax, centres, widths, heights, *, color, ...)` | draw one series honoring `hist_style` |
| `reload_prefs()` | drop the cache after the YAML is rewritten |

## Override path (one-off vs. persistent)

The YAML is **global and persistent** — editing it changes every plot. Every
accessor also takes `override=`:

- `override is None` → the persistent prefs lookup (the default path).
- a non-`None` `override` → used verbatim **for that call only**, touching
  neither the preferences YAML nor any mirror of it.

This is the explicit, auditable channel for a one-off "just this figure" tweak
and keeps the invariant intact: a literal color/marker still never appears at a
call site — it enters only through a named `override=` argument.

**Decision rule:**
- a persistent request ("from now on…") → run `/ndp-plot-prefs`.
- a one-off request ("just this plot…") → pass `override=` at the call site,
  change nothing on disk.

## The no-literal-colors invariant (audited)

No plotted color may remain as a string literal in a consuming script.
Audit a refactored plotting module with:

```bash
grep -nE "color=" scripts/diagnostic_plots.py | grep -v "ps\.color_for"   # → empty
grep -nE "#[0-9A-Fa-f]{6}" scripts/diagnostic_plots.py                     # → empty
grep -nE "fmt="  scripts/diagnostic_plots.py                               # → empty
```

`fmt=` (which bundles a marker letter) is replaced by
`marker=marker_for(role), linestyle="none"`. Non-color structural params
(markersize, capsize, alpha, linewidth) are sourced via `style(...)` where
shared, or kept as call-site literals where they legitimately vary per plot.

## Inputs

A preferences YAML (`schema_version` `1.0`) resolved as described above; the
bundled default is used when nothing else is present.

## Dependencies

- Python 3.10+
- PyYAML
- matplotlib (only for `apply_rcparams` / `draw_hist`)
