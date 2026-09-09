---
name: ndp-model
description: "Take a theorist's model description (a GENIE tune, a reweighting formula, an external event sample, or a shipped generator curve), express it as an NDP model YAML, run it through generator -> detector surrogate -> data comparison for a channel, and report the scorecard with its caveats. Use when the user says 'test my model', 'how does X compare to MINERvA', 'run this tune against the data', 'what if MEC were 30% larger'."
argument-hint: "[model description or models/<file>.yaml] [--channel <name>]"
---

# ndp-model — test a theory model against data

## 1. Understand what the theorist is handing you

Map the description onto one model `kind` (see `ndp/theory/models.py`):

| they give you | kind | you write |
|---|---|---|
| a GENIE tune name (+ optional custom tune XML dir) | `genie` | `tune`, `generator_list`, `n_events`, `n_jobs`, `seed`, optional `gxmlpath` |
| "scale / suppress / enhance <process> as a function of <kinematics>" | `reweight` | `base: reference_mc`, `weight_expr` in the numpy namespace (columns, observables, QE/RES/DIS/COH/MEC) |
| a Python function of the truth table | `reweight` | `weight_module: path.py:function` |
| their own event file (NDP npz or GENIE gst) | `external` | `path`, `format`, `sigma_per_nucleon_cm2` if they know the flux-averaged total |
| a generator curve already in the paper release | `shipped_curve` | `curve` |

Ask only if the mapping is ambiguous (e.g. weights that could be per nucleon or per nucleus).
Never invent a normalisation: if the sample's absolute cross section is unknown, say the folded
comparison will be shape-only-invalid and the unfolded one impossible, and ask.

## 2. Pick the channel and the measurement

`python -m ndp channels`. Use the channel the data belong to; check its `status` and the
`decided/default/open` flags on the physics fields. If a needed value is `open`, stop and ask.

The comparison is made on a *measurement* (`python -m ndp measurements --channel <channel>`):
the channel's `published` grid (default; folded + unfolded) or a user-defined observable pair
(folded only). If the user describes an observable that is not there ("compare in muon energy",
"in reconstructed Q²"), write `measurements/<channel>/<slug>.yaml` with `x`/`y` (`observable` =
truth name or expression, `reco` = reco name or expression over the cache columns, edges, units),
`status: user`, then `python -m ndp surrogate build --channel <channel> --measurement <slug>` and
check the printed closure is exact. Edges are the user's choice — ask if they gave none, or
propose some and label them clearly as a proposal.

## 3. Write the model YAML and validate

Put it in `models/<slug>.yaml` with `name`, `kind`, `description`, `author`, then
`python -m ndp validate models/<slug>.yaml --channel <channel>`.

## 4. Run

`python -m ndp run models/<slug>.yaml --channel <channel> [--measurement <name>]`. GENIE productions
are cached by spec fingerprint under `runs/_generator_cache/`; a 180k-event ME sample takes ~4 minutes
on 6 cores. If the folded comparison is skipped for lack of a surrogate, build one
(`python -m ndp surrogate build --channel <channel> --measurement <name>`) — it needs the caches
(`python -m ndp data cache --channel <channel>`; `python -m ndp data status` says if they are stale).

## 5. Report

Open `runs/<id>/report.md` and relay, in this order: (1) what the model is, which measurement
(observables, grid) and how it was normalised; (2) the folded result — data/pred ratio, −2lnL/ndf,
Pearson χ²/ndf, which surrogate and its closure; (3) on the published grid, the unfolded row(s) —
χ²_total/ndf with N cells, χ²_shape/ndf **with α**, norm offset — and where it lands among the
shipped curves; (4) the figures; (5) every warning in the manifest (measurement notes included).
State the reference frame and the normalisation constants used (they are in `channel.json`). Do
not declare a model "better" — the physicist decides; you give them both numbers and the plots.
