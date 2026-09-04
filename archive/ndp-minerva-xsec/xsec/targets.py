"""Tracker fiducial nucleon count (normalization target T).

Port of `PlotUtils::TargetUtils::GetTrackerNNucleons(5980, 8422, isMC=true,
apothem=850)` via the frozen exploration-repo `tracker_n_nucleons`
(dsigma_dpt.py:138), which reproduces the published N = 3.23e30 to +0.16 %.

Constants are pinned to MINERvA's TargetUtils values (Avogadro, molar masses,
areal density), NOT CODATA-latest — the cross-section normalization must match
the count MINERvA actually used, so reproduction fidelity is the requirement,
not the most recent external value. (Particle properties elsewhere in the
package still come from the pdg API; these are detector/target constants.)

  mass   = nPlanes * hexArea(apothem) * arealDensity
  N      = Σ_el  (mass * massFrac_el * N_A / molarMass_el) * (Z_el + N_el)
  with N_el = molarMass_el - Z_el (TargetUtils convention), so the per-element
  term reduces to mass * massFrac_el * N_A and N ≈ mass * N_A.
"""
import math

# --- detector / target geometry (provenance in comments) ---------------------
# nPlanes for z in [5980, 8422] mm: modules 27-80, 2 planes each
# (NSFDefaults.h:90-91); verified vs published N = 3.23e30 ± 1.4 %.
N_PLANES = 108
APOTHEM_MM = 850.0
AREAL_DENSITY_G_PER_MM2 = 1.9872 / 100.0   # 1.9872 g/cm^2 (TargetUtils.h:31)
AVOGADRO = 6.0221412927e23                 # TargetUtils.h:26 (pinned)

# isMC=true path -> NXMassFractionMC (TargetUtils.h:68-78)
MASS_FRACTION_MC = {
    "C": 0.8896, "H": 0.07533, "O": 0.02432, "Ti": 0.00697,
    "Al": 0.001613, "Si": 0.001613, "Cl": 0.000062, "N": 0.000651,
}
MOLAR_MASS = {  # g/mol (TargetUtils.h:140-153)
    "C": 12.0107, "H": 1.00794, "O": 15.9994, "Ti": 47.867,
    "Al": 26.982, "Si": 28.0855, "Cl": 35.453, "N": 14.007,
}
PROTONS_PER_ATOM = {  # Z (TargetUtils.h:172-186)
    "C": 6, "H": 1, "O": 8, "Ti": 22, "Al": 13, "Si": 14, "Cl": 17, "N": 7,
}

N_NUCLEONS_PUBLISHED = 3.23e30   # arXiv:2106.16210
N_NUCLEONS_GATE = 0.02           # |ours/published - 1| must be <= 2 %


def hex_area_mm2(apothem_mm=APOTHEM_MM):
    """Flat-top regular-hexagon area from the apothem (TargetUtils.cxx:28)."""
    return apothem_mm ** 2 * 6.0 / math.sqrt(3.0)


def tracker_fiducial_mass_g(n_planes=N_PLANES, apothem_mm=APOTHEM_MM):
    return n_planes * hex_area_mm2(apothem_mm) * AREAL_DENSITY_G_PER_MM2


def tracker_n_nucleons(n_planes=N_PLANES, apothem_mm=APOTHEM_MM):
    """Return (N_nucleons, fiducial_mass_g) for the tracker fiducial volume."""
    mass_g = tracker_fiducial_mass_g(n_planes, apothem_mm)
    total = 0.0
    for el, frac in MASS_FRACTION_MC.items():
        atoms = mass_g * frac * AVOGADRO / MOLAR_MASS[el]
        protons = atoms * PROTONS_PER_ATOM[el]
        neutrons = atoms * (MOLAR_MASS[el] - PROTONS_PER_ATOM[el])
        total += protons + neutrons
    return total, mass_g
