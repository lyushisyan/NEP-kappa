Capabilities
==============

NEP-kappa organizes force calculations, force constants, transport, and result
analysis through one YAML workflow. Force providers and transport solvers are
configured separately.

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
QHA+SSCHA evaluates this renormalization at each QHA equilibrium volume.

The command ``scph`` names the Phonopy stochastic SSCHA route implemented here;
it is not an ALAMODE perturbative SCPH calculation. Transport may use the
renormalized FC2 with the configured FC3 (and FC4 for coupled four-phonon
transport). These approximations have their own convergence requirements.

Execution and analysis
------------------------

Workflow presets reduce the normal interface to ``init`` and ``run``.
Stage commands support reuse and debugging. Slurm supports force arrays, LBTE,
and FourPhonon submission; available scheduling differs by stage.
Cached force jobs are reused only when their recorded inputs match.

Plotting and comparison consume completed phonon/transport data; convergence
studies prepare isolated parameter sweeps. Provenance and reports retain
settings, input identities, outputs, and execution state.

Choose a template in :doc:`examples` and consult :doc:`input_files`
for supported combinations and limits.
