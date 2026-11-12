NEP-kappa
========================

NEP-kappa uses a unified YAML input to connect atomic-force calculations,
force-constant generation, phonon properties, lattice thermal transport and
nonequilibrium phonon dynamics.

- **Calculation backends**: NEP, VASP, MACE and other ASE-compatible calculators.
- **Phonons and transport**: integrates Phonopy, phono3py and FourPhonon.
- **Nonequilibrium dynamics**: supports spatially homogeneous and periodic
  one-dimensional time-dependent BTE calculations.

New users should first complete the two examples in Chapters 3 and 4, then
consult the parameter reference and topic-specific chapters as needed.

.. toctree::
   :maxdepth: 1
   :numbered: 2
   :caption: User guide

   introduction
   installation
   starting_phonons
   starting
   workflow
   input_files
   temperature
   thermal_transport
   tdbte
   tdbte_spatial

.. toctree::
   :maxdepth: 1
   :caption: Project and support

   reference
   contact
   Releases <https://github.com/lyushisyan/NEP-kappa/releases>

`GitHub source <https://github.com/lyushisyan/NEP-kappa>`_ ·
`Online English documentation <https://nep-kappa.readthedocs.io/en/latest/>`_
