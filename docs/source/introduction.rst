Capabilities
==============

Each calculation uses a YAML input to select the force calculator,
force-constant method, and transport solver.

Calculators and force constants
---------------------------------

- NEP through Calorine CPUNEP.
- VASP with configured executable and POTCAR.
- MACE, external ASE factories, and registered calculator plugins.
- FC2/FC3 through phono3py finite displacements or HiPhive fitting.
- Thirdorder FC3 and Fourthorder FC4 for ShengBTE/FourPhonon workflows.

NEP, VASP, and suitable ASE calculators support relaxation. An ASE calculator
used for cell relaxation must provide stress as well as energy and forces.

Thermal transport
-------------------

- phono3py three-phonon RTA and iterative LBTE.
- phono3py SMM19 Wigner transport.
- FourPhonon three-plus-four-phonon transport with selectable solver.
- Optional isotope scattering and configured boundary scattering.
- Film/wire geometry normalization during plotting; raw HDF5 is unchanged.

Temperature-dependent phonons
-------------------------------

QHA uses isotropic volume scans to obtain thermal expansion and thermodynamic
properties. Fixed-volume Phonopy SSCHA generates temperature-renormalized FC2.
The ``qha-sscha`` route runs SSCHA at the QHA equilibrium volume for each
temperature, without reoptimizing that volume. The resulting FC2 is an auxiliary
harmonic matrix, not a free-energy Hessian. Optional bubble postprocessing
calculates input-FC3 diagonal on-shell frequency shifts in separate output files;
FC2 and transport remain unchanged. See :doc:`input_files` for approximation
details.

The command ``scph`` names the Phonopy stochastic SSCHA route implemented here;
it is not an ALAMODE perturbative SCPH calculation. Transport may use the
renormalized FC2 with the configured FC3 (and FC4 for coupled four-phonon
transport). These approximations have their own convergence requirements.

Execution and analysis
------------------------

``run`` executes a preset or custom stage list. Individual stage commands reuse
existing results. Slurm submission is available for force arrays, LBTE, and
FourPhonon.
Cached force jobs are reused only when their recorded inputs match.

``plot`` and ``compare`` read existing results. ``converge`` prepares parameter
sweeps. ``report`` collects results and recorded calculation settings.

Experimental dynamics
-----------------------

The ``tdbte`` stage evolves homogeneous three-phonon populations at fixed
frequencies from a prebuilt energy-shell kernel. It records energy,
equilibrium-control and entropy diagnostics, but physical rate normalization
is not independently validated. It does not construct kernels from FC2/FC3
or model laser absorption. See :doc:`tdbte` before using this research feature.
