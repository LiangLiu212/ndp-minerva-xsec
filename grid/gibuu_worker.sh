#!/bin/bash
# gibuu_worker.sh — HTCondor worker for a GiBUU generation campaign (jobsub-lite skill skeleton
# conventions: -R payload, -O pnfs out, %04d process dirs, exit 4 if nothing copied, DONE).
#
#   -R <dir>    payload dir on CVMFS (RCDS-published; @TAR_DIR@), falls back to $INPUT_TAR_FILE
#   -O <dir>    PNFS output base (required); this process writes to <dir>/<%04d PROCESS>/
#   -C <stem>   card template stem: $PAYLOAD/cards/<stem>.job.tmpl (placeholders __PATH_TO_INPUT__,
#               __FLUX_FILE__, __SEED__ left by ndp.theory.gibuu.write_card; everything else filled)
#   -S <int>    base seed; this process runs with seed = S + PROCESS (a 32-bit Fortran integer)
#   -F <file>   flux file name in $PAYLOAD/cards (default flux_gibuu.dat)
#
# The payload holds GiBUU.x, the shared libraries it needs (lib/), the buuinput tables and the card.
set -e
date; hostname; uname -r

PAYLOAD_OVERRIDE=""; PNFS_OUTPUT_DIR=""; CARD=""; SEED_BASE=""; FLUX="flux_gibuu.dat"
while getopts "R:O:C:S:F:" opt; do
  case $opt in
    R) PAYLOAD_OVERRIDE=$OPTARG ;;
    O) PNFS_OUTPUT_DIR=$OPTARG ;;
    C) CARD=$OPTARG ;;
    S) SEED_BASE=$OPTARG ;;
    F) FLUX=$OPTARG ;;
    *) echo "bad flag: $opt" >&2; exit 2 ;;
  esac
done
for v in PNFS_OUTPUT_DIR CARD SEED_BASE; do
  if [ -z "${!v}" ]; then echo "missing required flag for $v" >&2; exit 2; fi
done
PROCESS_STR=$(printf "%04d" "${PROCESS:-0}")
SEED=$(( SEED_BASE + ${PROCESS:-0} ))
echo "CLUSTER=${CLUSTER} PROCESS=${PROCESS} SEED=${SEED} CARD=${CARD}"

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
export LD_LIBRARY_PATH="${PAYLOAD}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
GIBUU_X="${PAYLOAD}/GiBUU.x"
[ -x "${GIBUU_X}" ] || { echo "no GiBUU.x in the payload" >&2; exit 3; }
ldd "${GIBUU_X}" | grep -i "not found" && { echo "missing shared libraries" >&2; exit 3; } || true

# ── The card for this process (run in $PWD: GiBUU writes everything to the cwd) ────────────────
WORK="${PWD}/gibuu"; mkdir -p "${WORK}"; cd "${WORK}"
sed -e "s|__PATH_TO_INPUT__|${PAYLOAD}/buuinput|" -e "s|__FLUX_FILE__|${PAYLOAD}/cards/${FLUX}|" -e "s|__SEED__|${SEED}|" \
    "${PAYLOAD}/cards/${CARD}.job.tmpl" > job.card
if grep -q "__[A-Z_]*__" job.card; then echo "unfilled placeholder in job.card" >&2; grep "__[A-Z_]*__" job.card; exit 2; fi
T0=$(date +%s)
set +e
"${GIBUU_X}" < job.card > gibuu.log 2>&1
RC=$?
set -e
WALL=$(( $(date +%s) - T0 ))
N_ROWS=0; N_EVENTS=0
if [ -s FinalEvents.dat ]; then
  N_ROWS=$(grep -vc '^#' FinalEvents.dat || true)
  N_EVENTS=$(awk '!/^#/ {print $1":"$2}' FinalEvents.dat | sort -u | wc -l)
fi
STATUS="ok"; [ "${RC}" -eq 0 ] && [ "${N_EVENTS}" -gt 0 ] || STATUS="failed"
echo "GiBUU rc=${RC} wall=${WALL}s rows=${N_ROWS} events=${N_EVENTS} status=${STATUS}"
tail -3 gibuu.log

# ── Outputs: the event file (gzipped), the cross-section and flux checks, the card, the log, a sidecar ─
OUT="${PWD}/out"; mkdir -p "${OUT}"
[ -s FinalEvents.dat ] && gzip -c FinalEvents.dat > "${OUT}/FinalEvents.dat.gz"
for f in neutrino_absorption_cross_section_ALL.dat neutrino_initialized_energyFlux.dat neutrino_absorption_cross_section_numbers.dat; do
  [ -f "$f" ] && cp "$f" "${OUT}/"
done
cp job.card "${OUT}/job.card"; gzip -c gibuu.log > "${OUT}/gibuu.log.gz"
CARD_SHA=$(sha256sum job.card | cut -c1-64)
XSEC=$(awk '!/^#/ && NF>1 {v=$2} END {print v}' neutrino_absorption_cross_section_ALL.dat 2>/dev/null || echo null)
cat > "${OUT}/manifest_${PROCESS_STR}.json" <<EOF
{"process": "${PROCESS:-0}", "cluster": "${CLUSTER:-}", "seed": ${SEED}, "card": "${CARD}", "card_sha256": "${CARD_SHA}",
 "payload": "${PAYLOAD}", "status": "${STATUS}", "returncode": ${RC}, "wall_s": ${WALL}, "n_rows": ${N_ROWS}, "n_events": ${N_EVENTS},
 "xsec_file_1e-38cm2": ${XSEC:-null}, "host": "$(hostname)", "finished": "$(date -u +%Y-%m-%dT%H:%M:%SZ)"}
EOF
ls -alh "${OUT}"
DEST="${PNFS_OUTPUT_DIR%/}/${PROCESS_STR}"
echo "ifdh mkdir_p ${DEST}"
ifdh mkdir_p "${DEST}"
N_COPIED=0
for f in "${OUT}"/*; do
  [ -e "$f" ] || continue
  ifdh cp -D "$f" "${DEST}/"
  N_COPIED=$((N_COPIED + 1))
done
if [ "${N_COPIED}" -eq 0 ]; then echo "ERROR: nothing copied to ${DEST}" >&2; exit 4; fi
echo "copied ${N_COPIED} file(s) to ${DEST}"
date
echo "DONE"
[ "${STATUS}" = "ok" ] || exit 5
