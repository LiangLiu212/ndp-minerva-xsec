#!/bin/bash
# worker_skeleton.sh — template HTCondor worker for the jobsub-lite skill.
#
# Copy this into your project, fill every `xxx` in the PROJECT-SPECIFIC
# sections, chmod +x, then submit:
#
#   python3 <skill>/scripts/jobsub.py submit --worker myworker.sh -N 10 \
#       --tar-label my_payload --expected-lifetime 8h --memory 2000MB \
#       -- -R @TAR_DIR@ -O @PNFS_OUT@
#
# Conventions this skeleton implements (see references/worker-patterns.md):
#   -R <dir>   payload dir (RCDS-published CVMFS path, substituted for @TAR_DIR@);
#              falls back to $INPUT_TAR_FILE when the payload came via dropbox://
#   -O <dir>   PNFS output base (required); each process writes to <dir>/<%04d>/
#   seed       CLUSTER*100000 + PROCESS — unique per process across submissions,
#              reproducible; needs a 64-bit seed consumer (see worker-patterns.md)
#   `DONE`     the final line on success (only after ≥1 output copied);
#              completion counting greps for it
set -e
date; hostname; uname -r

PAYLOAD_OVERRIDE=""
PNFS_OUTPUT_DIR=""
while getopts "R:O:" opt; do          # append your own flags to this getopts string
  case $opt in
    R) PAYLOAD_OVERRIDE=$OPTARG ;;
    O) PNFS_OUTPUT_DIR=$OPTARG ;;
    *) echo "bad flag: $opt" >&2; exit 2 ;;
  esac
done
if [ -z "${PNFS_OUTPUT_DIR}" ]; then
  echo "-O <pnfs_output_dir> is required" >&2
  exit 2
fi

PROCESS_STR=$(printf "%04d" "${PROCESS:-0}")
# CLUSTER*100000 + PROCESS: plain CLUSTER+PROCESS collides across submissions
# whenever Δcluster < N. Seeds land near 1e12–1e13 (cluster ids are ~1e7–1e8),
# so the consumer must take 64-bit seeds — int32 programs truncate silently.
SEED=$((CLUSTER * 100000 + PROCESS))
echo "CLUSTER=${CLUSTER} PROCESS=${PROCESS} SEED=${SEED}"

# ── ifdh (data transfer) — the standard FIFE pattern ──────────────────────────
source /cvmfs/larsoft.opensciencegrid.org/setup-env.sh
spack load --first ifdhc

# ── Locate the payload ────────────────────────────────────────────────────────
if [ -n "${PAYLOAD_OVERRIDE}" ]; then
  PAYLOAD="${PAYLOAD_OVERRIDE%/}"
elif [ -n "${INPUT_TAR_FILE}" ]; then
  if [ -d "${INPUT_TAR_FILE}" ]; then
    PAYLOAD="${INPUT_TAR_FILE%/}"
  else
    PAYLOAD=$(dirname "${INPUT_TAR_FILE}")
  fi
else
  echo "no -R override and no INPUT_TAR_FILE — payload unavailable" >&2
  exit 3
fi
echo "PAYLOAD=${PAYLOAD}"

# ── PROJECT-SPECIFIC: environment ─────────────────────────────────────────────
# Set up your software from $PAYLOAD (read-only; run in $PWD = condor scratch).
# For a relocatable conda/pixi payload, point LD_LIBRARY_PATH at its lib dirs —
# baked-in RPATHs won't resolve on the worker and the loader falls through to
# LD_LIBRARY_PATH. Example shape:
#   export PATH="${PAYLOAD}/xxx/bin:${PATH}"
#   export LD_LIBRARY_PATH="${PAYLOAD}/xxx/lib:${LD_LIBRARY_PATH}"
xxx

# ── PROJECT-SPECIFIC: run ─────────────────────────────────────────────────────
# Run in $PWD (per-job scratch, wiped after). Use $SEED for reproducible RNG and
# $PROCESS_STR in output names. Keep logs terse — they are fetched per process.
#   xxx --seed "${SEED}" -o "out_${PROCESS_STR}.root"
xxx

ls -alh

# ── Copy outputs to PNFS ──────────────────────────────────────────────────────
DEST="${PNFS_OUTPUT_DIR%/}/${PROCESS_STR}"
echo "ifdh mkdir_p ${DEST}"
ifdh mkdir_p "${DEST}"
# xxx: adjust the glob(s) to your real outputs. Zero matches exits 4 so a
# wrong glob can't count as success (DONE would otherwise still print).
N_COPIED=0
for f in *.root; do
  [ -e "$f" ] || continue
  ifdh cp -D "$f" "${DEST}/"
  N_COPIED=$((N_COPIED + 1))
done
if [ "${N_COPIED}" -eq 0 ]; then
  echo "ERROR: no files matched the output glob(s) — nothing copied to ${DEST}" >&2
  exit 4
fi
echo "copied ${N_COPIED} file(s) to ${DEST}"

date
echo "DONE"
