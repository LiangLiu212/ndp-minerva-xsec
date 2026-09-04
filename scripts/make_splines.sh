#!/usr/bin/env bash
# Generate GENIE cross-section splines for a tune with the in-repo GENIE, one gmkspl per target
# nuclide in parallel, merged with gspladd. A theorist's own tune (via GXMLPATH) needs exactly this.
#
#   pixi run make-splines -- --tune G18_10a_02_11b [--targets 1000060120,1000010010,...] \
#        [--probe 14] [--emax 100] [--knots 100] [--genlist CC] [--out external/splines/<name>.xml]
set -euo pipefail
: "${GENIE:?run via pixi}"
TUNE=G18_02a_00_000; PROBE=14; EMAX=100; KNOTS=100; GENLIST=CC
TARGETS=1000060120,1000010010,1000080160,1000220480,1000170350,1000130270,1000140280   # MINERvA tracker nuclides
OUT=""
while [ $# -gt 0 ]; do case "$1" in
  --tune) TUNE=$2; shift 2;; --probe) PROBE=$2; shift 2;; --emax) EMAX=$2; shift 2;; --knots) KNOTS=$2; shift 2;;
  --genlist) GENLIST=$2; shift 2;; --targets) TARGETS=$2; shift 2;; --out) OUT=$2; shift 2;;
  *) echo "unknown option $1" >&2; exit 2;; esac; done
OUT="${OUT:-$NDP_EXTERNAL/splines/${TUNE}_${GENLIST}_nu${PROBE}_e${EMAX}_n${KNOTS}.xml}"
W="$(dirname "$OUT")/work_${TUNE}_${GENLIST}"; mkdir -p "$W"
echo "tune $TUNE genlist $GENLIST probe $PROBE targets $TARGETS Emax $EMAX knots $KNOTS -> $OUT"
pids=(); parts=()
IFS=, read -ra TGT <<< "$TARGETS"
for t in "${TGT[@]}"; do
  part="$W/spl_${t}.xml"; parts+=("$part")
  if [ -s "$part" ]; then echo "have $part"; continue; fi
  ( gmkspl -p "$PROBE" -t "$t" -e "$EMAX" -n "$KNOTS" --tune "$TUNE" --event-generator-list "$GENLIST" \
           --message-thresholds Messenger_laconic.xml -o "$part" > "$W/gmkspl_${t}.log" 2>&1 \
      || { echo "gmkspl failed for $t (see $W/gmkspl_${t}.log)" >&2; exit 1; } ) &
  pids+=($!)
done
for p in "${pids[@]}"; do wait "$p"; done
gspladd -f "$(IFS=,; echo "${parts[*]}")" -o "$OUT" > "$W/gspladd.log" 2>&1
echo "wrote $OUT ($(grep -c '<spline ' "$OUT") splines)"
