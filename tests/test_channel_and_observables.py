import numpy as np
from _helpers import toy_truth
from ndp.channels import load_channel, list_channels
from ndp.channels import observables as obs


def test_channels_load():
    names = list_channels()
    assert "minerva_me_cc_inclusive_ptpz" in names
    ch = load_channel("minerva_me_cc_inclusive_ptpz")
    assert ch.binning.n_cells == 224 and ch.binning.formula == "ix*n_y + iy"


def test_phase_space_and_signal_masks():
    ch = load_channel("minerva_me_cc_inclusive_ptpz")
    t = toy_truth(5000)
    sig = ch.is_signal(t); ps = ch.in_phase_space(t)
    assert sig.all()
    th = obs.lep_theta(t, "beam"); pz = obs.lep_pz(t, "beam")
    inside_kin = (th <= np.deg2rad(20)) & (pz >= 1.5)
    assert np.all(ps <= inside_kin)  # phase space also needs the fiducial vertex


def test_beam_rotation_is_a_rotation():
    t = toy_truth(1000)
    assert np.allclose(obs.lep_p(t, "beam"), obs.lep_p(t, "detector"))
    assert not np.allclose(obs.lep_pT(t, "beam"), obs.lep_pT(t, "detector"))


def test_e_avail_definition():
    t = toy_truth(400, with_fs=True)
    ea = obs.E_avail(t)
    pdg, E = t["fs_pdg"], t["fs_E"]
    manual = np.zeros(400); idx = t.fs_event_index()
    for i, (p, e) in enumerate(zip(pdg, E)):
        if p == 2212: manual[idx[i]] += e - 0.93827208816
        elif abs(p) == 211: manual[idx[i]] += e - 0.13957039
        elif p == 111: manual[idx[i]] += e
        # neutrons contribute nothing
    assert np.allclose(ea, manual)
    assert np.allclose(obs.q3(t) ** 2, t["Q2"] + obs.q0(t) ** 2)


def test_vertex_cut_skipped_for_samples_without_geometry():
    ch = load_channel("minerva_me_cc_inclusive_ptpz")
    t = toy_truth(3000)
    with_geom = ch.in_phase_space(t)
    t.meta["has_geometry"] = False
    without = ch.in_phase_space(t)
    assert without.sum() > with_geom.sum()            # the fiducial-vertex requirement no longer applies
    th = obs.lep_theta(t, "beam"); pz = obs.lep_pz(t, "beam")
    assert np.array_equal(without, (th <= np.deg2rad(20)) & (pz >= 1.5))
