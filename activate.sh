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

# GENIE Reweight (in-place build beside the Generator; NUISANCE/nusystematics find it via GENIE_REWEIGHT)
export GENIE_REWEIGHT_VERSION="${GENIE_REWEIGHT_VERSION:-R-1_04_02}"
export GENIE_REWEIGHT="$NDP_EXTERNAL/genie/Reweight"
# nusystematics (CMake install prefix; systematicstools + the standalone fhicl-cpp suite share it)
export NUSYST_VERSION="${NUSYST_VERSION:-v02_00_07}"
export NUSYST="$NDP_EXTERNAL/nusystematics/install"
export NUSYST_SRC="$NDP_EXTERNAL/nusystematics/nusystematics"
export nusystematics_ROOT="$NUSYST" systematicstools_ROOT="$NUSYST" fhicl_cpp_standalone_ROOT="$NUSYST"
# fhicl lookups (DumpConfiguredTweaksNuSyst etc.) search FHICL_FILE_PATH: the working directory, the
# installed systematicstools examples and the nusystematics source fcl/ (not installed by its CMake).
export FHICL_FILE_PATH=".:$NUSYST/fcl/fcl:$NUSYST_SRC/fcl${FHICL_FILE_PATH:+:$FHICL_FILE_PATH}"
# nusystematics' GENIE tools resolve their tune from GENIE_XSEC_TUNE (the FNAL genie_xsec UPS convention);
# default to the platform's G18_02a tune (matches the CVMFS spline default above), override per sample.
export GENIE_XSEC_TUNE="${GENIE_XSEC_TUNE:-G18_02a_00_000}"
# NUISANCE (CMake install prefix; its apps find data/ and parameters/ through $NUISANCE)
export NUISANCE="$NDP_EXTERNAL/nuisance/install"
export NUISANCE_SRC="$NDP_EXTERNAL/nuisance/nuisance"

export PATH="$GENIE/bin:$GENIE_REWEIGHT/bin:$ROOTEGPythia6_ROOT/bin:$NUWRO/bin:$GIBUU/objects:$ACHILLES/bin:$NUSYST/bin:$NUISANCE/bin:$PATH"
# NUISANCE and ACHILLES each install a libnuhepmc_cpputils.so (different NuHepMC cpputils versions, same
# soname): NUISANCE's must come first; ACHILLES keeps its own through the RPATH baked into its binaries.
export LD_LIBRARY_PATH="$GENIE/lib:$GENIE_REWEIGHT/lib:$ROOTEGPythia6_ROOT/lib:$NUWRO/bin:$NUSYST/lib:$NUISANCE/lib:$ACHILLES/lib:$ACHILLES/lib64:$CONDA_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
