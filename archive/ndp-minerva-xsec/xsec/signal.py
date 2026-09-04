"""Truth signal definition and phase space for the CC-inclusive measurement.

Two deliberately distinct predicates
------------------------------------
``is_signal``      — WHAT process we measure: a muon neutrino charged-current
                     interaction. Two branches only: mc_incoming == nu_mu and
                     mc_current == CC. Nothing about where it happened or
                     where the muon went.

``in_phase_space`` — WHERE the measurement is defined: the kinematic and
                     geometric window, in TRUE variables, inside which the
                     detector has usable acceptance and the cross section is
                     reported. True vertex in the tracker fiducial volume
                     (z in [5980, 8422] mm, hexagonal apothem 850 mm — this
                     is also what makes the target "hydrocarbon": the volume
                     selects the scintillator, and the nucleon count T is
                     computed for the same volume), true muon angle <= 20 deg
                     to the NuMI beam, true muon p_z >= 1500 MeV/c.

Where each applies (source-confirmed, MINERvA-101 + MAT Cutter):
- signal/background split of reco-selected MC (numerator, migration,
  background categories):           is_signal ONLY        (Cutter.cxx:133-137,
                                    runEventLoop.cpp:134)
- efficiency denominator:           is_signal AND in_phase_space
                                    (Cutter.cxx:145-151 isEfficiencyDenom)

Keeping the phase space OUT of the signal split is load-bearing: a true
nu_mu-CC event outside the true fiducial that passes the reco selection is
treated as acceptance (corrected by the efficiency denominator), NOT as
background to subtract (ExtractCrossSection.cpp:197). The golden counts
(43539 signal / 104 background; 397604 truth-signal rows of which 65712 in
phase space on MC run 110040) are only reproducible with this split.

Provenance: predicates mirror exploring/dsigma_dpt.py:300-330 (frozen 1D
chain) == CCInclusiveSignal.h / runEventLoop.cpp:371-377; constants live in
xsec/constants.py; branch lists come from config/branches.json (roles
``signal_definition`` and ``phase_space``).

All functions are vectorized: NumPy arrays in, boolean masks out, no I/O.
"""
import numpy as np

from . import catalog
from .constants import (CC_CURRENT, MAX_MU_THETA_RAD, NU_MU_PDG,
                        PZ_MIN_TRUE_MEV, Z_MAX_MM, Z_MIN_MM)
from .kinematics import in_hex_apothem, true_theta_p

# Branch lists from the catalog (single source of truth; code cannot drift
# from config/branches.json).
SIGNAL_BRANCHES = catalog.branches_for_role("signal_definition")
PHASE_SPACE_BRANCHES = catalog.branches_for_role("phase_space")


def is_signal(mc_incoming, mc_current):
    """Truth signal definition: nu_mu charged-current, nothing else.

    Parameters are the per-event arrays of the two SIGNAL_BRANCHES (from the
    MC reco tree when classifying selected events, or from the Truth tree
    when building the efficiency denominator).
    """
    return ((np.asarray(mc_incoming) == NU_MU_PDG)
            & (np.asarray(mc_current) == CC_CURRENT))


def in_phase_space(mc_vtx, mc_primFSLepton):
    """Truth phase space: the window the cross section is defined in.

    Parameters
    ----------
    mc_vtx : (n, 4) array — true vertex (x, y, z, t) in mm; only x, y, z used.
    mc_primFSLepton : (n, 4) array — true primary FS lepton 4-momentum
        (px, py, pz, E) in MeV; only the 3-momentum is used.

    Cuts (boundary semantics match the frozen chain exactly):
      z inclusive on both ends; hexagon strict <; theta <= 20 deg
      (CCInclusiveSignal.h:68); p_z(beam) >= 1500 MeV/c, where
      p_z(beam) = p * cos(theta) with theta measured to the beam axis.
    """
    vtx = np.asarray(mc_vtx, dtype=np.float64)
    lep = np.asarray(mc_primFSLepton, dtype=np.float64)

    z_ok = (Z_MIN_MM <= vtx[:, 2]) & (vtx[:, 2] <= Z_MAX_MM)
    hex_ok = in_hex_apothem(vtx[:, 0], vtx[:, 1])

    theta, p = true_theta_p(lep[:, 0], lep[:, 1], lep[:, 2])
    theta_ok = theta <= MAX_MU_THETA_RAD
    pz_ok = p * np.cos(theta) >= PZ_MIN_TRUE_MEV

    return z_ok & hex_ok & theta_ok & pz_ok


def is_efficiency_denominator(mc_incoming, mc_current, mc_vtx, mc_primFSLepton):
    """Efficiency-denominator predicate: is_signal AND in_phase_space.

    The ONLY place the two predicates combine (Cutter::isEfficiencyDenom).
    Run over the Truth tree.
    """
    return (is_signal(mc_incoming, mc_current)
            & in_phase_space(mc_vtx, mc_primFSLepton))
