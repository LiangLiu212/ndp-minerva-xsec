# Environment snapshot

Copy of the **live** pixi workspace manifest + lock from
`/home/feanor/ndp-genesis-agent/` (the parent directory of this repo), where the
environment actually lives. This repo deliberately has no `pixi.toml` at its root so
that pixi auto-discovers the parent workspace.

Snapshot taken 2026-06-12, after plan Step 1 (added `pip` and the `pdg` pypi
dependency; PDG edition 2026). On top of this env, `runlog_tools` is installed
editably: `pixi run pip install -e ~/runlog_tools`
(from `git@github.com:LiangLiu212/runlog_tools.git`).

To recreate elsewhere: place `pixi.toml` + `pixi.lock` in a workspace root,
`pixi install`, then the editable runlog_tools install above.

If the parent manifest changes, re-copy both files here in the same commit as the
change that needed them.
