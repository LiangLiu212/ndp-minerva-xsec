"""Tests for xsec.targets: nucleon count must match the published value."""
import math

from xsec import targets


def test_nucleon_count_within_gate():
    n, mass_g = targets.tracker_n_nucleons()
    rel = n / targets.N_NUCLEONS_PUBLISHED - 1.0
    assert abs(rel) <= targets.N_NUCLEONS_GATE, f"N={n:.4e} rel={rel:+.3%}"
    # frozen exploration result was +0.16 %
    assert abs(rel - 0.0016) < 5e-4, f"rel {rel:+.4%} != frozen +0.16%"


def test_fiducial_mass_scale():
    # ~5.4 tonne tracker fiducial (paper quotes 5.48 t for the reported volume)
    mass_g = targets.tracker_fiducial_mass_g()
    assert 5.0e6 < mass_g < 5.8e6


def test_hex_area_formula():
    # flat-top hexagon: area = 2*sqrt(3)*a^2 = a^2 * 6/sqrt(3)
    a = 850.0
    assert math.isclose(targets.hex_area_mm2(a), 2 * math.sqrt(3) * a * a)


def test_nucleons_scale_with_planes_and_area():
    n1, _ = targets.tracker_n_nucleons(n_planes=108)
    n2, _ = targets.tracker_n_nucleons(n_planes=54)
    assert math.isclose(n2, n1 / 2, rel_tol=1e-12)
