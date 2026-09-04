# Sourced by pixi on `pixi shell` / `pixi run` ([activation] in pixi.toml).
# Exports the in-repository GENIE build and its dependencies. Paths are relative to the
# project so a clone elsewhere works unchanged; nothing here is machine-specific except the
# optional CVMFS spline default.

export NDP_EXTERNAL="$PIXI_PROJECT_ROOT/external"

# Pythia6 + ROOT's TPythia6 interface (removed from ROOT in 6.30; conda ROOT has neither)
export ROOTEGPythia6_ROOT="$NDP_EXTERNAL/ROOTEGPythia6/install"
export PYTHIA6="$ROOTEGPythia6_ROOT/lib"          # libPythia6.so + libEGPythia6.so (GENIE's --with-pythia6-lib)

# GENIE Generator, built in place
export GENIE_VERSION="${GENIE_VERSION:-R-3_06_02}"
export GENIE="$NDP_EXTERNAL/genie/Generator"

# GENIE's C++ dependencies from the conda environment
export LHAPDF6_INC="$CONDA_PREFIX/include"
export LHAPDF6_LIB="$CONDA_PREFIX/lib"
export LHAPATH="$CONDA_PREFIX/share/LHAPDF"      # PDF-set directory (GENIE's default GRV98LO is internal; sets optional)
export LOG4CPP_INC="$CONDA_PREFIX/include"
export LOG4CPP_LIB="$CONDA_PREFIX/lib"
export LIBXML2_INC="$CONDA_PREFIX/include/libxml2"
export LIBXML2_LIB="$CONDA_PREFIX/lib"
export ROOTSYS="${ROOTSYS:-$CONDA_PREFIX}"

# The unprefixed conda-forge compiler packages (gcc/gxx/gfortran) put plain gcc/g++/gfortran on
# PATH but export no CC/CXX/FC; the build scripts and CMake read these.
export CC="${CC:-$(command -v x86_64-conda-linux-gnu-gcc 2>/dev/null || command -v gcc)}"
export CXX="${CXX:-$(command -v x86_64-conda-linux-gnu-g++ 2>/dev/null || command -v g++)}"
export FC="${FC:-$(command -v x86_64-conda-linux-gnu-gfortran 2>/dev/null || command -v gfortran)}"

# Pre-computed cross-section splines: default to the CVMFS G18_02a set when it is mounted
# (same tune family the platform channel uses); override GENIEXSECFILE to use another.
_cvmfs_spl=/cvmfs/larsoft.opensciencegrid.org/products/genie_xsec/v3_06_00/NULL/G1802a00000-k250-e1000/data/gxspl-NUsmall.xml
if [ -z "${GENIEXSECFILE:-}" ] && [ -f "$_cvmfs_spl" ]; then
  export GENIEXSECFILE="$_cvmfs_spl"
fi
unset _cvmfs_spl

# NuWro (in-place build; data/ found through $NUWRO; event1.so lives in bin/)
export NUWRO="$NDP_EXTERNAL/nuwro"
# GiBUU release 2025 (binary objects/GiBUU.x; job cards set path_to_input = $GIBUU_INPUT)
export GIBUU="$NDP_EXTERNAL/gibuu/release2025"
export GIBUU_INPUT="$NDP_EXTERNAL/gibuu/buuinput"
# ACHILLES (CMake install prefix; data in share/Achilles)
export ACHILLES_VERSION="${ACHILLES_VERSION:-v0.3.1}"
export ACHILLES="$NDP_EXTERNAL/achilles/install"
export ACHILLES_SRC="$NDP_EXTERNAL/achilles/Achilles"

export PATH="$GENIE/bin:$ROOTEGPythia6_ROOT/bin:$NUWRO/bin:$GIBUU/objects:$ACHILLES/bin:$PATH"
export LD_LIBRARY_PATH="$GENIE/lib:$ROOTEGPythia6_ROOT/lib:$NUWRO/bin:$ACHILLES/lib:$ACHILLES/lib64:$CONDA_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
