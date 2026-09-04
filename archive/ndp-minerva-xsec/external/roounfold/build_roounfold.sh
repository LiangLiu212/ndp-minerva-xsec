#!/usr/bin/env bash
# Build a minimal RooUnfold shared library (kBayes / D'Agostini path only) from
# the vendored MINERvA UnfoldUtils sources, for the PyROOT-based 2D-unfolding
# cross-check (mirrors MnvUnfold::UnfoldHisto2D).
#
# No environment install: uses the conda ROOT compiler already in the pixi env
# (root-config --cxx). Run inside the pixi env:
#     pixi run bash external/roounfold/build_roounfold.sh
# Reproducible: re-run to rebuild. Outputs libRooUnfoldMin.so (+ _rdict.pcm) here.
#
# -DNOTUNFOLD drops the ROOT-TUnfold and (absent) Dagostini algorithms, which the
# Bayes path does not use. The other RooUnfold::New() algorithm references
# (SVD/Invert/BinByBin, +vendored TSVDUnfold) ARE compiled in: ROOT resolves
# symbols at load time, so the base RooUnfold.cxx (which names them in New())
# needs them defined even though we only ever construct RooUnfoldBayes.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$(cd "${ROOUNFOLD_SRC:-$HERE/../../../MinervaExpt/UnfoldUtils/RooUnfold}" && pwd)"
echo "RooUnfold sources: $SRC"

CXX="$(root-config --cxx)"
CFLAGS="$(root-config --cflags)"
LIBS="$(root-config --libs)"
DEFS="-DNOTUNFOLD"
SRCS=(RooUnfold RooUnfoldResponse RooUnfoldBayes RooUnfoldErrors
      RooUnfoldInvert RooUnfoldBinByBin RooUnfoldSvd TSVDUnfold matrix_mult_mt)

echo "1/2 generating dictionary (rootcling)"
rm -f "$HERE/roounfold_dict.cxx" "$HERE/roounfold_dict_rdict.pcm"
rootcling -f "$HERE/roounfold_dict.cxx" $DEFS -I"$SRC" \
    RooUnfold.h RooUnfoldResponse.h RooUnfoldBayes.h RooUnfoldErrors.h \
    RooUnfoldInvert.h RooUnfoldBinByBin.h RooUnfoldSvd.h \
    "$HERE/RooUnfoldMinLinkDef.h"

echo "2/2 compiling libRooUnfoldMin.so ($CXX)"
CXXSRCS=()
for s in "${SRCS[@]}"; do CXXSRCS+=("$SRC/$s.cxx"); done
# shellcheck disable=SC2086
"$CXX" -fPIC -shared -w $CFLAGS $DEFS -I"$SRC" \
    "${CXXSRCS[@]}" "$HERE/roounfold_dict.cxx" \
    -o "$HERE/libRooUnfoldMin.so" $LIBS

echo "OK -> $HERE/libRooUnfoldMin.so"
