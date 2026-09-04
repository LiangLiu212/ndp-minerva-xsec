# archive/

Earlier work preserved verbatim, with its git history (subtree-merged, so `git log --follow`
reaches the original commits).

- `ndp-minerva-xsec/` — pure-Python (uproot/awkward/numpy) reproduction of the MINERvA ME FHC
  inclusive CC νμ d²σ/dpT dp∥ (arXiv:2106.16210) from the Open Data AnaTuples, with the full
  systematic ladder (flux universes, GENIE universes, muon energy scale, RPA, 2p2h, MINOS
  efficiency), a RooUnfold cross-check of the D'Agostini stage, and Fig-8/Fig-13 style results.
  Formerly the GitHub repository `LiangLiu212/ndp-minerva-xsec` (renamed to `ndp-platform`).
  Its own `README.md` describes its four-stage workflow; nothing in the platform imports it yet
  (see `docs/roadmap.md` — its systematics machinery is the natural source for surrogate
  uncertainty bands).
