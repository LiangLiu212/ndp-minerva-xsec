"""Native-frame handling: generator samples (neutrino along +z) are already in the NuMI beam frame and
must not be rotated again; detector-frame tables (MINERvA caches, toys without meta) keep rotating."""
import numpy as np

from ndp.channels import observables as obs
from ndp.channels.signal import leading_proton, rotate_to_frame
from ndp.events import TruthTable
from _helpers import toy_truth


def _beam_native(n=500):
    t = toy_truth(n, seed=7, with_fs=True)
    t.meta.update({"generator": "GiBUU", "has_geometry": False})   # no explicit frame: inferred as beam
    return t


def test_native_frame_inference():
    assert obs.native_frame(toy_truth(50)) == "detector"                      # toy: neither key
    t = _beam_native()
    assert obs.native_frame(t) == "beam"
    t.meta["frame"] = "detector"                                               # an explicit frame wins
    assert obs.native_frame(t) == "detector"
    assert obs.frame_rotation_angle("detector", "beam") == obs.NUMI_BEAM_ANGLE_RAD
    assert obs.frame_rotation_angle("beam", "detector") == -obs.NUMI_BEAM_ANGLE_RAD
    assert obs.frame_rotation_angle("beam", "beam") == 0.0


def test_beam_native_sample_is_not_rotated_in_beam_frame():
    t = _beam_native()
    px, py, pz = obs._lep_p3(t, "beam")
    assert np.array_equal(py, t["lep_py"]) and np.array_equal(pz, t["lep_pz"])
    th_beam = obs.lep_theta(t, "beam")
    th_raw = np.arctan2(np.hypot(t["lep_px"], t["lep_py"]), t["lep_pz"])
    assert np.allclose(th_beam, th_raw)
    # asking for the detector frame rotates by the inverse angle: a forward muon tilts by 3.37 deg
    th_det = obs.lep_theta(t, "detector")
    fwd = np.argmin(th_raw)
    assert abs(np.degrees(th_det[fwd] - th_raw[fwd])) < 3.4 and abs(np.degrees(th_det[fwd])) > 3.0
    # protons follow the same rule
    lp_beam = leading_proton(t, "beam", {"theta_max_deg": 180.0})
    lp_raw = leading_proton(TruthTable(dict(t.columns), {"source": "toy-detector"}), "detector", {"theta_max_deg": 180.0})
    ok = lp_beam["index"] >= 0
    assert np.allclose(lp_beam["theta"][ok], lp_raw["theta"][ok])


def test_detector_native_sample_keeps_rotating():
    t = toy_truth(300, seed=3)
    px, py, pz = obs._lep_p3(t, "beam")
    a = obs.NUMI_BEAM_ANGLE_RAD
    assert np.allclose(py, t["lep_py"] * np.cos(a) - t["lep_pz"] * np.sin(a))
    assert np.allclose(pz, t["lep_py"] * np.sin(a) + t["lep_pz"] * np.cos(a))
    x, y, z = rotate_to_frame(np.zeros(3), np.zeros(3), np.ones(3), "beam")            # native defaults to detector
    assert np.allclose(np.degrees(np.arctan2(np.hypot(x, y), z)), 3.373, atol=0.01)
    x, y, z = rotate_to_frame(np.zeros(3), np.zeros(3), np.ones(3), "beam", native="beam")
    assert np.allclose(z, 1.0) and np.allclose(y, 0.0)
