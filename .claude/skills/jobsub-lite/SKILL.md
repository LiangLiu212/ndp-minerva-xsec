---
name: jobsub-lite
description: Submit, monitor, and harvest FNAL FIFE grid jobs with jobsub_lite, including RCDS/CVMFS tarball publishing. Use when the user wants to run grid jobs, submit to HTCondor at Fermilab, ship a code payload via RCDS/dropbox tarball, publish to CVMFS fifeuser repos, check jobsub job status, fetch grid logs, or pull outputs from PNFS with ifdh. Triggers - jobsub, jobsub_lite, jobsub_submit, FNAL grid, FIFE, grid campaign, RCDS, dropbox tarball, CVMFS publish, fifeuser, ifdh, PNFS outputs, condor grid jobs.
---

# jobsub-lite — FNAL grid submission with RCDS publishing

Bundled CLI: `scripts/jobsub.py` (stdlib-only). **Always run it with system
`python3` and NEVER through pixi/conda** — jobsub_lite is its own Python venv
and inherited `PYTHONPATH`/`CONDA_*`/`PIXI_*` break it (the lib scrubs those
before every jobsub/ifdh call, but don't tempt fate by wrapping the CLI itself).

```bash
SKILL=~/.claude/skills/jobsub-lite       # or <project>/.claude/skills/jobsub-lite
python3 $SKILL/scripts/jobsub.py <subcommand> ...
```

## Adoption model

- The skill directory is **read-only and copyable**: to pin it to a project or
  share via a repo, copy the whole `jobsub-lite/` dir into
  `<project>/.claude/skills/`.
- All mutable state lives in `<project>/.jobsub/` (config.json, catalog.json,
  tarballs/, runs/) — created by `init`, self-gitignored, never inside the skill.
- Commands find `.jobsub/` by walking up from CWD, so run them from anywhere
  inside the project (or set `$JOBSUB_STATE_DIR`).

## Workflow

1. **First use in a project**: `jobsub.py init --group <exp>` (e.g. `dune`).
   Probes the jobsub_lite binaries, writes `.jobsub/config.json`, reports auth.
   Valid kerberos (`kinit`) is required; the bearer token is auto-fetched at
   submit time. Ask the user for the group if it isn't obvious — never guess.

2. **Decide how the payload travels** (details: `references/rcds.md`):
   - a few small files → `-f /pnfs/...` or `-f file:///abs/path` per file;
   - a software install / toolchain reused across jobs → **tarball + RCDS
     publish** (lands *unpacked* on CVMFS; workers run it in place — no
     per-job transfer or untar);
   - one-shot big tarball, no reuse → `submit --tar-file X.tar` (dropbox
     without cataloging).

3. **Build the payload tarball** (cached on a full-tree content key — path,
   mtime, size + include/exclude rules; rebuilds are free no-ops when nothing
   changed):
   ```bash
   python3 $SKILL/scripts/jobsub.py tarball build \
       --build-dir /path/to/project --include env bin data \
       --exclude-component .git --name-prefix myproj
   ```
   Multi-GB payloads: add `--background` and watch the `.log` next to the tar.

4. **Publish to RCDS + catalog it** (runs ONE real grid job — the sentinel —
   to discover the CVMFS path; typically 5–20 min):
   ```bash
   python3 $SKILL/scripts/jobsub.py publish --tarball <tar> --label my_payload
   python3 $SKILL/scripts/jobsub.py verify --label my_payload   # before reuse
   ```
   RCDS garbage-collects after ~30 days — `verify` warns at 21 d and says
   `republish` at 28 d or when the dir is gone.

5. **Write the worker**: copy `templates/worker_skeleton.sh`, fill the two
   PROJECT-SPECIFIC sections (env + run), keep the conventions — seed is
   `CLUSTER * 100000 + PROCESS` (64-bit; unique across submissions), outputs
   go to `<pnfs-out>/<%04d PROCESS>/` via `ifdh cp -D`, zero copied outputs
   exits 4, and the last line on success is a bare `DONE` (completion
   counting greps for it). Patterns + pitfalls: `references/worker-patterns.md`.

6. **Smoke test with `-N 1` first — always** — then scale:
   ```bash
   python3 $SKILL/scripts/jobsub.py submit --worker w.sh -N 1 \
       --tar-label my_payload --expected-lifetime 2h --memory 2000MB \
       --output-suffix .root -- -R @TAR_DIR@ -O @PNFS_OUT@
   ```
   `@TAR_DIR@` → the published CVMFS dir; `@PNFS_OUT@` → the PNFS output dir
   (default `<pnfs_scratch_base>/$USER/jobsub-lite/<runtype>/<stem>`, override
   with `--pnfs-out`). `--dry-run` builds+records the full command with
   `--no_submit` and submits nothing. Flag reference: `references/cheatsheet.md`.

7. **Track / debug / harvest**:
   ```bash
   python3 $SKILL/scripts/jobsub.py status <jobid>     # re-polls jobsub_q
   python3 $SKILL/scripts/jobsub.py list --active
   python3 $SKILL/scripts/jobsub.py fetchlog <jobid>   # worker logs -> record dir
   python3 $SKILL/scripts/jobsub.py pull <jobid> --suffix .root
   python3 $SKILL/scripts/jobsub.py cancel <jobid>
   ```
   Every submission writes `<stem>.{command.json,submit.log,gridlog}` under
   `.jobsub/runs/<runtype>-<date>/`; the jobid printed at submit is the handle
   for everything else.

## Hard rules

- **Scale needs consent**: submitting real jobs costs shared grid resources.
  A 1-node smoke test after the user asked for grid work is fine; a campaign
  (`-N` large) needs the user to have named or approved the scale.
- **`--no_submit` must precede `file://`** in a raw jobsub_submit command —
  after it, everything is worker args (the CLI's `--dry-run` handles this).
- Don't hand-edit `.jobsub/catalog.json`; use `publish`/`adopt`/`verify`.
- Arch-tuned binaries (`-march=native` builds) SIGILL on older workers — either
  rebuild generic or set `append_condor_requirements` to pin a floor, e.g.
  `'(TARGET.Microarch>="x86_64-v3")'` (see `references/worker-patterns.md`).
- If jobsub commands fail with odd Python errors, check for a poisoned env
  (`PYTHONPATH`, conda/pixi vars) — run from a clean shell.

## Files

- `scripts/jobsub.py` — CLI (init | tarball | publish | verify | labels |
  adopt | submit | status | list | cancel | fetchlog | pull)
- `scripts/lib/` — the machinery (submit records, RCDS publish, monitors)
- `templates/worker_skeleton.sh` — starting point for every new worker
- `templates/publish_only.sh` — RCDS sentinel (used by `publish`; don't edit)
- `references/cheatsheet.md` — jobsub_lite flags, worker env vars, auth
- `references/rcds.md` — how RCDS works, publish flow, GC, payload decision
- `references/worker-patterns.md` — worker anatomy, seeds, ifdh, portability
