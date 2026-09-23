"""Single entry point for plot colors, markers, and style preferences.

Reads a user-owned preferences YAML (see `default_prefs_path()` for how it is
located) and exposes role-based accessors. This module
is intended to be the ONLY way a color or marker enters a plot: plotting code
asks for `color_for("data")`, never writes `"#000000"`.

Series ROLES drive every lookup, not call sites:
    data, mc_signal, mc_bkg, efficiency, purity, ratio, reference

Per-call override:
    Every accessor takes `override=`. `override is None` -> the persistent
    prefs lookup (default). A non-None override wins for that one call and
    touches neither the YAML nor the mirror — the explicit, auditable
    channel for a one-off "just this figure" tweak. Persistent changes go
    through `/ndp-plot-prefs`.

Usage:
    from scripts import plot_style as ps
    ax.errorbar(x, y, color=ps.color_for("data"), marker=ps.marker_for("data"),
                linestyle="none")
    ax.bar(x, h, color=ps.color_for("mc_signal", override="#999999"))  # one-off
"""

import os
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]
PREFS_PATH = REPO / "docs" / "style_preferences.yaml"

#: Env var naming an explicit prefs YAML. Highest precedence after an explicit
#: `use_prefs_file()` call, so this module works for a user whose project does
#: not have this repo's `docs/` layout at all.
PREFS_PATH_ENV = "NDP_STYLE_PREFS"

# Built-in default bundled with the tool. Loaded when neither the env var nor a
# project docs/style_preferences.yaml is present, so plotting works on a fresh
# checkout without first running a skill. NOTE: these are one project's seeded
# choices (solid-dark palette, step histograms, magma sequential) chosen to
# reproduce this repo's pre-refactor figures -- NOT the Okabe-Ito colorblind-safe
# set, which is carried separately as the `okabe_ito` swatch table.
DEFAULT_PREFS_PATH = Path(__file__).resolve().parent / "style_preferences.default.yaml"

# Active prefs file (swappable for tests / alternate configs via
# use_prefs_file) and a tiny path-keyed cache.
_prefs_path = None
_cache = {}


def default_prefs_path():
    """Where prefs come from absent an explicit `use_prefs_file()`.

    Resolution order: ``$NDP_STYLE_PREFS`` → the project's
    ``docs/style_preferences.yaml`` → the bundled default shipped next to this
    module.
    """
    env = os.environ.get(PREFS_PATH_ENV)
    if env:
        return Path(env)
    if PREFS_PATH.exists():
        return PREFS_PATH
    return DEFAULT_PREFS_PATH


def _active_prefs_path():
    """The file to actually read.

    An explicit `use_prefs_file()` path wins and is returned as-is even when
    missing, so `open()` raises a clear error rather than silently falling back.
    Otherwise resolve through `default_prefs_path()`.
    """
    if _prefs_path is not None:
        return _prefs_path
    return default_prefs_path()


def _prefs():
    """Load and cache the active YAML prefs (falling back to the bundled default
    when docs/style_preferences.yaml is not present)."""
    path = _active_prefs_path()
    key = str(path)
    if key not in _cache:
        with open(path) as fh:
            _cache[key] = yaml.safe_load(fh)
    return _cache[key]


def reload_prefs():
    """Drop the cached prefs (call after the harness rewrites the YAML)."""
    _cache.clear()


def use_prefs_file(path):
    """Point the module at a different prefs YAML (tests / alternate configs).

    Pass `None` (or `plot_style.PREFS_PATH`) to restore the normal resolution
    order — `default_prefs_path()`: env var, then project docs/, then the
    bundled default. Clears the cache.
    """
    global _prefs_path
    _prefs_path = None if path is None or Path(path) == PREFS_PATH else Path(path)
    _cache.clear()


# --------------------------------------------------------------------------
# Color / marker accessors (role-based, override-aware)
# --------------------------------------------------------------------------

def color_for(role, override=None):
    """Color for a series role. `override` (a color) wins for this call only."""
    if override is not None:
        return override
    roles = _prefs()["palette"]["roles"]
    try:
        return roles[role]
    except KeyError:
        raise KeyError(
            f"no color for role {role!r}; known roles: {sorted(roles)}"
        )


def marker_for(role, override=None):
    """Point marker for a role. `override` wins for this call only.

    Bar roles (mc_signal, mc_bkg, purity) have no marker; asking for one
    is a programming error unless overridden.
    """
    if override is not None:
        return override
    markers = _prefs().get("markers", {})
    try:
        return markers[role]
    except KeyError:
        raise KeyError(
            f"no marker for role {role!r}; markered roles: {sorted(markers)}"
        )


def palette(override=None):
    """The full {role: color} mapping (override replaces it wholesale)."""
    if override is not None:
        return override
    return dict(_prefs()["palette"]["roles"])


def okabe_ito(override=None):
    """The named Okabe-Ito swatch table (orange/sky/green/...)."""
    if override is not None:
        return override
    return dict(_prefs()["palette"]["okabe_ito"])


# --------------------------------------------------------------------------
# Structural style accessors
# --------------------------------------------------------------------------

def hist_style(override=None):
    """Histogram drawing style: 'bar' | 'step' | 'step_filled'."""
    if override is not None:
        return override
    return _prefs().get("hist_style", "bar")


def sequential_palette(override=None):
    """Colormap name for sequential / 2D plots."""
    if override is not None:
        return override
    return _prefs().get("sequential_palette", "viridis")


def style(key, override=None):
    """A non-color structural knob from `error_bars` / `reference_line` /
    `ratio_band` (e.g. 'efficiency_markersize', 'guide_linestyle', 'lo').

    Looked up across those sections; first hit wins. `override` short-circuits.
    """
    if override is not None:
        return override
    prefs = _prefs()
    for section in ("error_bars", "reference_line", "ratio_band", "syst_band"):
        block = prefs.get(section, {})
        if key in block:
            return block[key]
    raise KeyError(f"unknown style key {key!r}")


def ratio_panel(override=None):
    """The `ratio_panel` block (height ratios, ylabels, ylim, model lw).

    `override` (a dict) is merged over the persistent values for this call.
    """
    base = dict(_prefs().get("ratio_panel", {}))
    if override:
        base.update(override)
    return base


#: Fallbacks used when a prefs file omits the `syst_band` block entirely.
_SYST_BAND_DEFAULTS = {"alpha": 0.25, "ratio_alpha": 0.20, "edgecolor": "none"}


def syst_band(override=None):
    """The `syst_band` block (fill alphas, edgecolor) — section-scoped.

    Use this rather than `style("alpha")`: that helper searches several
    sections and returns the first hit, which for the key `alpha` is
    `ratio_band`'s value, not this band's. `override` (a dict) is merged over
    the persistent values for this call.
    """
    base = dict(_SYST_BAND_DEFAULTS)
    base.update(_prefs().get("syst_band", {}))
    if override:
        base.update(override)
    return base


def response_prefs():
    """The response-style block ({results, summary_verbosity, ...})."""
    return dict(_prefs().get("response", {}))


def apply_rcparams():
    """Apply the `rcparams` block to matplotlib's global rcParams.

    Empty block (the seeded default) leaves matplotlib untouched, preserving
    current output. The harness can populate font/figure defaults here.
    """
    rc = _prefs().get("rcparams") or {}
    if not rc:
        return
    import matplotlib
    matplotlib.rcParams.update(rc)


# --------------------------------------------------------------------------
# Histogram drawing helper — honors hist_style so every consumer renders
# the user's chosen style without re-deciding bar-vs-step at the call site.
# --------------------------------------------------------------------------

def draw_hist(ax, centres, widths, heights, *, color, alpha=1.0, label=None,
              width_scale=1.0, style_override=None):
    """Draw one histogram series on `ax` in the user's hist_style.

    'bar' -> filled matplotlib bars (current default). 'step' / 'step_filled'
    -> outline / filled step histogram via edges reconstructed from
    centres+widths. Color comes in already resolved (call color_for first).
    """
    import numpy as np
    centres = np.asarray(centres, dtype=float)
    widths = np.asarray(widths, dtype=float)
    heights = np.asarray(heights, dtype=float)
    hs = hist_style(override=style_override)

    if hs == "bar":
        return ax.bar(centres, heights, width=widths * width_scale,
                      color=color, alpha=alpha, label=label)

    # Reconstruct contiguous bin edges from centres/widths for step drawing.
    edges = np.empty(len(centres) + 1)
    edges[:-1] = centres - widths / 2
    edges[-1] = centres[-1] + widths[-1] / 2
    fill = hs == "step_filled"
    return ax.stairs(heights, edges, color=color, alpha=alpha, label=label,
                     fill=fill)
