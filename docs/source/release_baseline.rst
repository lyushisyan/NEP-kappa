NEP-kappa 2.0 scientific baseline
===================================

The 2.0 refactor is protected by a compact, machine-readable scientific
local baseline in ``benchmarks/scientific-baselines.yaml``. It freezes selected
published-workflow observables instead of committing large displacement or
transport directories.

The initial baseline covers:

- 3C-SiC DFT, material-specific NEP, and MACE 300 K RTA conductivity versus
  q-point mesh;
- 3C-SiC NEP QHA properties at representative temperatures;
- the completed 3C-SiC 300 K Phonopy-SSCHA smoke result;
- BAs DFT-versus-NEP FC2, FC3, and FC4 comparison metrics;
- BAs DFT and NEP 300 K three-plus-four-phonon conductivity.

The frozen manifest records a pending BAs MACE transport case. This is a
historical baseline status, not a live scheduler report. Pending entries are
not treated as reference results.

These files and ``tests/`` are ignored by Git and must be available locally
to run the baseline checks. Original BAs source summaries are retained under
``benchmarks/data/BAs/sources/`` so their hashes
can be checked without private ``calculations/`` data. Snapshot checks do not
recompute production observables; they verify the archived reference data.

When both local directories are available, run the baseline tests with::

   pytest -q tests/test_scientific_baselines.py

The baseline files record both strict snapshot tolerances and looser scientific
tolerances. Strict tolerances detect accidental data or parsing changes. The
scientific tolerances are the acceptance limits to use when an intentional
numerical-library or algorithm change is being evaluated.
