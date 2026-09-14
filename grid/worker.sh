#!/bin/bash
# worker.sh — HTCondor worker for the ndp-platform MINERvA Open Data campaign (jobsub-lite skill
# skeleton conventions: -R payload, -O pnfs out, %04d process dirs, exit 4 if nothing copied, DONE).
#
#   -R <dir>    payload dir on CVMFS (RCDS-published; @TAR_DIR@), falls back to $INPUT_TAR_FILE
#   -O <dir>    PNFS output base (required); this process writes to <dir>/<%04d PROCESS>/
#   -W <file>   worklist name (one root:// URL per line) — looked up in $CONDOR_DIR_INPUT, then $PAYLOAD/worklists
#   -K mc|data  kind
#   -C <name>   channel manifest name (channels/<name>.yaml in the payload)
#   -n <int>    files per process: this process handles worklist lines [PROCESS*n, (PROCESS+1)*n)
#   -t <sec>    XRootD request timeout (default 300)
#
# Each file is processed by `python -m ndp.grid.process_file` from the payload's slim pixi env;
# a failed file leaves a manifest_<tag>.json with status "failed" and the process continues.
set -e
date; hostname; uname -r

PAYLOAD_OVERRIDE=""; PNFS_OUTPUT_DIR=""; WORKLIST=""; KIND=""; CHANNEL=""; NPER=1; TIMEOUT=300
while getopts "R:O:W:K:C:n:t:" opt; do
  case $opt in
    R) PAYLOAD_OVERRIDE=$OPTARG ;;
    O) PNFS_OUTPUT_DIR=$OPTARG ;;
    W) WORKLIST=$OPTARG ;;
    K) KIND=$OPTARG ;;
    C) CHANNEL=$OPTARG ;;
    n) NPER=$OPTARG ;;
    t) TIMEOUT=$OPTARG ;;
    *) echo "bad flag: $opt" >&2; exit 2 ;;
  esac
done
for v in PNFS_OUTPUT_DIR WORKLIST KIND CHANNEL; do
  if [ -z "${!v}" ]; then echo "missing required flag for $v" >&2; exit 2; fi
done

PROCESS_STR=$(printf "%04d" "${PROCESS:-0}")
SEED=$((${CLUSTER:-0} * 100000 + ${PROCESS:-0}))
echo "CLUSTER=${CLUSTER} PROCESS=${PROCESS} SEED=${SEED} KIND=${KIND} CHANNEL=${CHANNEL} NPER=${NPER}"

# ── ifdh (data transfer) — the standard FIFE pattern ──────────────────────────
source /cvmfs/larsoft.opensciencegrid.org/setup-env.sh
spack load --first ifdhc

# ── Locate the payload ────────────────────────────────────────────────────────
if [ -n "${PAYLOAD_OVERRIDE}" ]; then
  PAYLOAD="${PAYLOAD_OVERRIDE%/}"
elif [ -n "${INPUT_TAR_FILE}" ]; then
  if [ -d "${INPUT_TAR_FILE}" ]; then PAYLOAD="${INPUT_TAR_FILE%/}"; else PAYLOAD=$(dirname "${INPUT_TAR_FILE}"); fi
else
  echo "no -R override and no INPUT_TAR_FILE — payload unavailable" >&2; exit 3
fi
echo "PAYLOAD=${PAYLOAD}"

# ── Environment: the payload's slim pixi env + the ndp package (read-only CVMFS) ──────────────
export PATH="${PAYLOAD}/env/bin:${PATH}"
export LD_LIBRARY_PATH="${PAYLOAD}/env/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export PYTHONPATH="${PAYLOAD}"
export PYTHONNOUSERSITE=1
export NDP_CONFIG="${PAYLOAD}/grid/ndp_grid.yaml"
export NDP_PAYLOAD="${PAYLOAD}"
export XRD_REQUESTTIMEOUT="${TIMEOUT}"
export XRD_STREAMTIMEOUT="${TIMEOUT}"
export MPLBACKEND=Agg
if [ -f "${PAYLOAD}/ndp_version.json" ]; then
  export NDP_GIT_SHA=$(sed -n 's/.*"sha": *"\([0-9a-f]*\)".*/\1/p' "${PAYLOAD}/ndp_version.json" | head -1)
fi
PY="${PAYLOAD}/env/bin/python3"
"$PY" -c "import sys, numpy, uproot, awkward, XRootD; print('python', sys.version.split()[0], 'numpy', numpy.__version__, 'uproot', uproot.__version__)"

# ── Worklist slice for this process ────────────────────────────────────────────
if [ -f "${CONDOR_DIR_INPUT:-/nonexistent}/${WORKLIST}" ]; then
  WL="${CONDOR_DIR_INPUT}/${WORKLIST}"
elif [ -f "${PAYLOAD}/worklists/${WORKLIST}" ]; then
  WL="${PAYLOAD}/worklists/${WORKLIST}"
elif [ -f "${WORKLIST}" ]; then
  WL="${WORKLIST}"
else
  echo "worklist ${WORKLIST} not found" >&2; exit 2
fi
START=$(( ${PROCESS:-0} * NPER + 1 ))
END=$(( START + NPER - 1 ))
mapfile -t URLS < <(sed -n "${START},${END}p" "${WL}" | grep -v '^\s*$')
echo "worklist ${WL}: lines ${START}-${END} -> ${#URLS[@]} file(s)"
if [ "${#URLS[@]}" -eq 0 ]; then echo "nothing to do for PROCESS=${PROCESS}"; echo "DONE"; exit 0; fi

# ── Process each file in scratch ($PWD) ────────────────────────────────────────
OUT="${PWD}/out"; mkdir -p "${OUT}"
N_FAIL=0
for url in "${URLS[@]}"; do
  echo "=== $(date '+%H:%M:%S') processing ${url}"
  set +e
  "$PY" -m ndp.grid.process_file --url "${url}" --channel "${CHANNEL}" --kind "${KIND}" --out "${OUT}" --timeout "${TIMEOUT}"
  rc=$?
  set -e
  if [ "$rc" -ne 0 ]; then echo "process_file exit ${rc} for ${url}"; N_FAIL=$((N_FAIL + 1)); fi
done
ls -alh "${OUT}"

# ── Copy outputs to PNFS ──────────────────────────────────────────────────────
DEST="${PNFS_OUTPUT_DIR%/}/${PROCESS_STR}"
echo "ifdh mkdir_p ${DEST}"
ifdh mkdir_p "${DEST}"
N_COPIED=0
for f in "${OUT}"/*.npz "${OUT}"/*.json; do
  [ -e "$f" ] || continue
  ifdh cp -D "$f" "${DEST}/"
  N_COPIED=$((N_COPIED + 1))
done
if [ "${N_COPIED}" -eq 0 ]; then
  echo "ERROR: no output files — nothing copied to ${DEST}" >&2
  exit 4
fi
echo "copied ${N_COPIED} file(s) to ${DEST}; ${N_FAIL} file(s) failed processing"
date
echo "DONE"
