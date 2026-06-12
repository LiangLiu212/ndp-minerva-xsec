"""Analysis constants with provenance.

Two kinds live here, both deliberately outside the analysis logic:

1. Particle properties — fetched from the PDG API at import (standing project
   rule: physics constants come from external databases, never typed in).
2. MINERvA/tutorial constants — experiment-specific numbers PDG cannot
   provide; each carries its provenance. Changing any of these is a physics
   decision (log it in docs/decisions.md when that file exists).

All angles in radians, lengths in mm, momenta in MeV/c unless noted.
"""
import math

import pdg

from . import catalog

# --- particle properties (PDG API, edition pinned by the pixi lock) --------
_api = pdg.connect()
NU_MU_PDG = int(_api.get_particle_by_name("nu_mu").mcid)    # 14
MUON_PDG = int(_api.get_particle_by_name("mu-").mcid)       # 13
MUON_MASS_GEV = float(_api.get_particle_by_name("mu-").mass)

# --- tuple conventions (config/branches.json _meta.conventions) ------------
CC_CURRENT = int(catalog.conventions()["mc_current"]["CC"])  # 1
NC_CURRENT = int(catalog.conventions()["mc_current"]["NC"])  # 2

# --- tracker fiducial volume (the "on hydrocarbon" region) -----------------
# Provenance: certified selector tools/cc_inclusive_selector.py (exploration
# repo) == MINERvA-101 runEventLoop.cpp:363 (minZ, maxZ, apothem).
Z_MIN_MM = 5980.0
Z_MAX_MM = 8422.0
APOTHEM_MM = 850.0
# Flat-top hexagon edges for the apothem test (CCInclusiveCuts.h:177,185-186).
APOTHEM_SLOPE = -1.0 / math.sqrt(3.0)
APOTHEM_INTERCEPT_MM = 2.0 * APOTHEM_MM / math.sqrt(3.0)

# --- muon kinematic limits --------------------------------------------------
# Reco cut is strict <, truth phase-space cut is <= (CCInclusiveSignal.h:68);
# no float input lands exactly on the boundary, so the distinction is
# documentation, not behavior (see selector's equivalent-mutant note).
MAX_MU_THETA_RAD = 20.0 * math.pi / 180.0
PZ_MIN_TRUE_MEV = 1500.0  # truth::PZMu >= 1500 MeV/c (runEventLoop.cpp:377)

# --- beam geometry ----------------------------------------------------------
# NuMI beam points downward in detector coordinates; true kinematics are
# evaluated in beam coordinates via RotateX(NUMI_BEAM_ANGLE_RAD).
# Provenance: dsigma_dpt.py:92 (frozen 1D chain) == TruthFunctions.h:21-27.
NUMI_BEAM_ANGLE_RAD = -0.05887
