# jobsub_lite cheatsheet

jobsub_lite is FNAL's thin wrapper over HTCondor (`/opt/jobsub_lite/bin/*`).
One `jobsub_submit` call creates one *cluster* of `-N` *processes*; each process
runs the worker script with `$PROCESS` = 0..N-1.

## Auth

- Kerberos first: `kinit user@FNAL.GOV` (check: `klist -s`).
- Bearer tokens (WLCG/SciTokens) are fetched automatically by jobsub_lite via
  `htgettoken` on submit; stored at `/run/user/$UID/bt_u$UID` (or
  `$BEARER_TOKEN_FILE`). X.509 proxies are legacy — don't build new flows on them.
- Group/VO comes from `-G <group>` (dune, uboone, ...); `--role Analysis` is
  the normal analysis role.

## jobsub_submit — flags that matter

| Flag | Meaning |
|---|---|
| `-G <group>` | experiment/VO (required) |
| `--role Analysis` | VOMS role |
| `-N <n>` | number of processes in the cluster |
| `--memory 2000MB` | RAM request; job is held if exceeded |
| `--disk 20GB` | scratch disk request |
| `--expected-lifetime 8h` | max walltime class; shorter = more slots |
| `-f <uri>` | per-job input file; `/pnfs/...` (server-side ifdh) or `file:///abs` — lands in `$CONDOR_DIR_INPUT` |
| `--tar_file_name dropbox:///abs/x.tar` | payload tarball via RCDS (see rcds.md) |
| `-e VAR=VAL` (or `-e VAR`) | export env var into the job |
| `--append_condor_requirements '<classad>'` | extra worker constraints, e.g. `'(TARGET.Microarch>="x86_64-v3")'` |
| `--singularity-image <path>` | override the default apptainer image (e.g. an EL9 image under `/cvmfs/singularity.opensciencegrid.org/fermilab/`) |
| `--no_submit` | build everything, submit nothing (dry run) |
| `file:///abs/worker.sh arg1 arg2` | the executable; **everything after it is worker args** |

**Gotcha — flag order:** jobsub_submit stops parsing its own options at the
`file://` executable. `--no_submit` (or any jobsub flag) placed after it is
silently passed to the worker and the job submits for real.

**Output parsing:** a successful submit prints the cluster id
`<cluster>.<proc>@<schedd>.fnal.gov` — that id is what `jobsub_q`,
`jobsub_rm`, and `jobsub_fetchlog` take via `--jobid`.

## Worker-side environment (set by the glidein)

| Var | Meaning |
|---|---|
| `$PROCESS` | 0..N-1 within the cluster (also `$PROCID`) |
| `$CLUSTER` | HTCondor cluster number (also `$CLUSTERID`) |
| `$CONDOR_DIR_INPUT` | where `-f` files land |
| `$INPUT_TAR_FILE` | path *into the unpacked tarball dir* when `--tar_file_name` was used (its dirname — or itself, if a dir — is the payload root) |
| `$INPUT_TAR_DIR_LOCAL` | worker-local copy of the unpacked tarball dir |
| `$_CONDOR_SCRATCH_DIR` | per-job scratch = CWD at start; wiped afterwards |
| `$GRID_USER` | submitting user |

Jobs run inside an apptainer container (AL9 default) on the OSG; assume no
network beyond CVMFS + ifdh-mediated transfer, and no `/exp`, no home dirs.

## Monitor / control

```bash
jobsub_q       -G dune --jobid <cluster>@<schedd>          # queue state
jobsub_q       -G dune --jobid <id> --long                 # full classads
jobsub_rm      -G dune --jobid <id>                        # kill
jobsub_fetchlog -G dune --jobid <id> --unzipdir DIR        # worker stdout/err
```

HTCondor `JobStatus` codes in `--long` output: 1 idle, 2 running, 3 removed,
4 completed, 5 held, 6 transferring.

Held jobs: almost always resource overrun (memory/disk/lifetime) — check
`jobsub_q --long` for `HoldReason`, bump the request, resubmit.

## Data movement (ifdh)

- Workers: `ifdh cp -D <files> <dir>/` copies into a directory; `ifdh mkdir_p`
  first. Never write directly to `/pnfs` paths with POSIX tools from a worker.
- Login nodes: same commands work for pulling outputs back; reads of
  `/pnfs/.../scratch` are fine, but dCache scratch is subject to eviction —
  move keeper outputs to persistent space or tape-backed areas.
- The skill CLI wraps the common cases: `pull <jobid> --suffix .root`.

## Environment poisoning (why the scrub exists)

jobsub_lite and ifdh are Python programs in their own venvs. Inherited
`PYTHONPATH`, `PYTHONHOME`, `CONDA_*`, `PIXI_*`, `MAMBA_*`, `VIRTUAL_ENV`
make them import the wrong interpreter/libs and crash in confusing ways.
`scripts/lib/submit_env.py` strips exactly those vars (auth vars pass through)
for every subprocess call the skill makes.
