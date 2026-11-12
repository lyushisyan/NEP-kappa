Spatial nonlinear TD-BTE
==========================

Scope and input
-----------------

An optional ``tdbte.spatial`` mapping enables a one-dimensional periodic real
space domain with the existing nonlinear, mode-coupled three-phonon collision
integral at every cell. The phonon q mesh remains three-dimensional. Omitting
this mapping retains the original homogeneous workflow and integration settings.
No RTA, Callaway, linearized evolution, four-phonon or electron coupling is used.

From the repository root, validate the complete example and then compute:

.. code-block:: bash

   nepkappa info examples/tdbte-spatial.yaml --for run
   nepkappa run examples/tdbte-spatial.yaml
   nepkappa plot examples/tdbte-spatial.yaml
   nepkappa report calculations/3C-SiC/tdbte-spatial

``run`` generates FC2/FC3 before dynamics. For existing force constants, use
the compact input below with ``nepkappa tdbte``. No force jobs are started by
that command. Paths resolve from the launch directory. Replace the source
directory with matching ``fc2.hdf5``, ``fc3.hdf5`` and ``phono3py_disp.yaml``:

.. code-block:: yaml

   tdbte:
     force-constants: calculations/3C-SiC/fc
     mesh: [3, 3, 3]
     temperature: 300
     duration-ps: 0.1
     max-step-ps: 0.01
     samples: 6
     spatial:
       length-nm: 64
       cells: 16
       direction: [1, 0, 0]
       boundary: periodic
       cfl: 0.8
       batch-size: 8
       initial:
         kind: ttg
         delta-temperature-K: 10
   output:
     result-dir: calculations/3C-SiC/spatial-ttg

A completed packaged kernel can instead be selected with ``tdbte.kernel``;
omit ``force-constants`` and ``mesh``. Its recorded matching force constants
remain necessary to calculate velocities. ``spatial.phonons`` may point to
relocated matching force constants; all three source hashes are checked.
Legacy kernels without native grid/source identity are insufficient for spatial
velocity preparation. No NAC data, shifted/generalized grids, or 1D/2D q meshes
are accepted by the spatial route.

The numbers above demonstrate the workflow, not convergence or physical rates.
Use distinct result directories: existing ``result-dir/tdbte`` is refused.

Initial conditions
----------------------

``temperature`` supplies the reference Bose state. The supported ``initial``
models are:

* ``kind: ttg``: :math:`T=T_0+\Delta T\cos(2\pi x/L)`, followed by a Bose
  distribution at each cell. ``delta-temperature-K`` is :math:`\Delta T`.
  Require at least four cells and :math:`|\Delta T|<T_0`.
* ``kind: gaussian``: thermal Gaussian pulse; specify ``delta-temperature-K``,
  ``sigma-nm`` and optionally ``center-nm`` (default :math:`L/2`). The code
  sums periodic images and normalizes the center to one. This matches a
  Gaussian when its width is small compared with the domain and remains smooth
  at the periodic seam for wider pulses.
* ``kind: selective``: multiply the reference Bose occupations of explicitly
  selected ``branches`` by :math:`1+\delta g(x)`. ``amplitude`` is the relative
  occupation increase :math:`\delta`, not kelvin or laser fluence. Specify
  ``sigma-nm`` and optionally ``center-nm``; the periodic Gaussian profile is
  used. ``profile: uniform`` excites these branches equally throughout space.
  Branches may alternatively be supplied by top-level ``tdbte.branches``.
* ``kind: uniform``: local Bose state at ``temperature + delta-temperature-K``.
  It is useful as a stationary numerical control.

For example, replace only ``spatial.initial`` to obtain these two pulses:

.. code-block:: yaml

   # Thermal pulse
   initial:
     kind: gaussian
     delta-temperature-K: 30
     sigma-nm: 8
     center-nm: 32

.. code-block:: yaml

   # Selective optical excitation for two-atom 3C-SiC
   initial:
     kind: selective
     branches: [3, 4, 5]
     amplitude: 0.2
     sigma-nm: 8
     center-nm: 32

Indices are frequency-ordered native branch IDs; they do not track polarization
through crossings and are not a universal optical-mode classifier. Excitation
of an isolated member of a degenerate group can depend on the eigenvector basis.

Streaming, collisions and units
----------------------------------

``direction`` is a Cartesian three-vector, normalized internally, independently
of primitive lattice axes. Velocities come from the analytic FC2 dynamical-matrix
derivative, projected into the same native eigenvectors as the collision
backend, in native BZ representatives. Reciprocal translations change
eigenvector phases; using only modulo-grid coordinates with those eigenvectors
would be inconsistent. The output stores both indexing addresses and native
fractional q coordinates. Mode ordering is q-major, native branch-minor.

The conversion from the native THz-Angstrom derivative to nm/ps is 0.1.
Fixed Gamma translations have zero occupation and zero streaming velocity.
The primitive-cell volume is in nm cubed. No vacuum, cross-sectional-area or
conductivity-plot normalization enters this bulk 3D-q occupation equation.

Each step performs half a nonlinear DOP853 collision interval, one full
first-order upwind streaming interval, and another half collision interval.
``max-step-ps`` caps the split step and the actual step also obeys
:math:`\max|v_x|\Delta t/\Delta x\leq\mathrm{cfl}\leq1`.
The upwind stage is a conservative convex combination, preserving positivity
under that condition. ``spatial.collisions: false`` disables only collisions
for a ballistic control.

Collision tolerances are ``spatial.rtol: 1e-11`` and ``atol: 1e-14`` by default.
If a DOP853 stage produces nonpositive active occupations, the collision
interval is retried from its original state with a smaller internal step,
at most eight halvings. Retries are counted; exhaustion raises an error.
No occupation clipping or scattering-rate rescaling is performed. This explicit
method can be expensive for stiff collisions; it does not implement an implicit
stiff solver.

Strang ordering does not make the full method second order: the chosen
forward-Euler upwind streaming step gives first-order time and space errors.
Numerical spatial diffusion must be checked against cell size and CFL.
``samples`` sets saved times, or ``spatial.snapshot-times-ps`` supplies an
increasing explicit list including zero and ``duration-ps``. Steps land on
saved times, so dense output requests can shorten split steps. Compare runs
with identical saved times when testing time-step convergence.

All cells share one immutable collision kernel. Each RHS evaluation opens
each event chunk once and processes cell batches inside that mapping.
``batch-size`` bounds temporary entropy arrays; it does not change physics.
The DOP853 state/stages still require memory proportional to cells times modes.
No dense mode collision matrix or cell-by-event array is built.

Outputs and continuity
--------------------------

``result-dir/tdbte/spatial.hdf5`` stores only requested snapshots. Its shape
conventions and names expose the units:

.. list-table:: Spatial datasets
   :header-rows: 1

   * - Dataset
     - Shape / interpretation
   * - ``occupation``
     - (saved times, cells, modes), dimensionless
   * - ``temperature_K``
     - (times, cells), temperature matching local thermal energy
   * - ``heat_flux_eV_nm2_ps``
     - (times, cells), physical cell-centered phonon heat flux
   * - ``branch_energy_eV_nm3``
     - (times, cells, branches), thermal energy density; no zero-point energy
   * - ``branch_heat_flux_eV_nm2_ps``
     - (times, cells, branches), physical branch heat flux
   * - ``branch_collision_power_eV_nm3_ps``
     - (times, cells, branches), instantaneous nonlinear collision energy exchange
   * - ``branch_face_flux_integral_eV_nm2``
     - (times, cells, branches), numerical RIGHT-face flux integrated since preceding snapshot
   * - ``branch_collision_energy_integral_eV_nm3``
     - (times, cells, branches), actual collision energy change over both split half steps
   * - ``branch_continuity_residual_eV_nm3``
     - (times, cells, branches), saved-interval discrete energy-budget residual
   * - ``departure_from_local_BE``
     - Energy-weighted occupation L1 difference from local energy-equivalent Bose state, divided by local energy
   * - ``temperature_harmonics_real_K`` / ``imag_K``
     - Four complex Fourier coefficients; validity mask excludes Nyquist/aliased harmonics
   * - ``ttg_cosine_amplitude_K``
     - Twice the real fundamental coefficient; initialized to the TTG temperature amplitude

The file also stores coordinates, frequencies, zero-mode mask, branch IDs,
Cartesian velocities, native q representatives, primitive cell/volume,
native eigenvectors, physical constants, source/kernel identity and the audit.
Multiply heat flux by :math:`1.602176634\times10^{11}` for W/m squared and
energy density by :math:`1.602176634\times10^8` for J/m cubed; both factors are
stored as attributes.

The branch budget uses exactly the face flux used by streaming:

.. math::

   e_s(t_j)-e_s(t_{j-1})
   +\nabla_x\!\cdot\!\int_{t_{j-1}}^{t_j}F_s^{\rm upwind}dt
   =\int_{t_{j-1}}^{t_j}P_s^{3\rm ph}dt.

Thus spatial transport and local scattering are separated without subtracting
two finite-difference derivatives of sparsely saved data. Physical
cell-centered heat flux and numerical upwind face flux are distinct datasets;
the former approaches the continuum quantity as the grid is refined.
``audit.json`` records energy, branch continuity, collision energy, entropy,
occupation positivity, CFL, RHS count and retries. Failure preserves diagnostic
files with an unsuccessful audit. ``nepkappa report`` includes these checks;
``nepkappa plot`` draws spatial temperature, heat flux, branch energies and
temperature harmonics without loading the full occupation dataset.

Validation and physical limits
---------------------------------

Run the opt-in numerical study from the repository root with a fresh directory:
the command requires a local copy of the script in ``tools/``, which is ignored
by Git and is not included in a fresh checkout.

.. code-block:: bash

   OMP_NUM_THREADS=1 RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
     python tools/validate_spatial_tdbte.py \
     --output calculations/spatial-tdbte-validation

The nonlinear six-mode algebraic model contains both :math:`2\leftrightarrow1+1`
and :math:`3\leftrightarrow1+2` events in all mirrored propagation directions.
It has arbitrary coefficients and is not SiC. Its Fourier diffusivity is an
independent closed-form reference derived from the Bose polynomials and heat
capacity, not fitted to time evolution. No auxiliary linearized evolution or
RTA solver is introduced. This study tests analytic ballistic advection,
spatial/time convergence, a weak-amplitude long-wave Fourier limit, Gaussian
pulses, selective excitation and three TTG amplitudes with paired ballistic
controls. Nonlinear occupation initialization and energy-to-temperature
conversion can generate signals even without collisions; those controls are
necessary before attributing harmonics to collision-induced nonlinearity.

For optional real SiC numerical demonstrations, add ``--kernel`` pointing to
a completed 3D packaged SiC kernel. ``--material-only`` skips the algebraic
studies. Real runs use 300 K, a 3D native q grid, eight spatial cells in 64 nm
and 0.02 ps duration, with both collisions and ballistic controls. They exercise
actual FC-based velocities and nonlinear collisions but do not establish
material relaxation rates or a converged material Fourier limit.

All outputs retain ``validated_for_physical_dynamics: false``. Complete absolute
rate, q-mesh, force-supercell, spatial and time convergence remain necessary
for quantitative material claims. The native tetrahedron steady-state operator
and this weak energy-shell operator differ at finite mesh; an unrelated
conductivity value must not be assumed to be a Fourier benchmark. Degenerate
off-diagonal velocity matrix elements are measured and recorded; scalar
populations retain diagonal native-basis velocities and omit coherences.
This basis approximation requires separate checks for degenerate-mode dynamics.
The off-shell small-q audit is optional (``--acoustic-check`` in the collision
audit tool), not a prerequisite for software development.

Recorded numerical results (2026-11-07)
-----------------------------------------

The complete opt-in report is available as
:download:`spatial validation JSON <_static/tdbte-spatial-validation-20261107.json>`.
These runs used single-threaded Rayon, OpenMP and OpenBLAS on Python 3.12,
phono3py/Phonopy 4.1.0 and the existing Numba event reductions.
The full repository regression passed 882 cases with 24 warnings and no failures
or skips. The strict Sphinx HTML build passed. Public checks include real
non-Bravais FC2 group velocities compared with frequency finite differences,
batch/scalar parity, single-cell homogeneous parity, HDF5/plot/report routing,
all initial profiles, branch continuity and the complete algebraic benchmark.

* Analytic collision-free advection errors decrease at approximately first
  order over 32, 64 and 128 cells, for positive and negative velocities.
* With the nonlinear toy at fixed spatial resolution, split steps 0.04, 0.02
  and 0.01 ps give occupation L2 errors 0.01089, 0.005514 and 0.002871 against
  the 0.00125 ps reference. The overall first-order temporal error is visible.
* The toy's independently derived Fourier diffusivity is 0.507074 nm squared
  per ps at 300 K. For a 64 nm TTG, 0.01 K perturbation and 20 ps duration,
  128, 256 and 512 cells give decay-rate errors of 41.86%, 17.26% and 4.95%.
  This is a convergence demonstration and a finite-grid 10% acceptance check,
  not a claim that the continuum limit is exactly reached.
* Thermal Gaussian, selective Gaussian and TTGs at 1, 30 and 90 K amplitudes
  pass energy, branch continuity, collision energy, entropy and positivity
  checks, each paired with a collision-free control in the toy model.
* TTG and Gaussian separately converge under spatial refinement and time-step
  refinement. Their occupation RMS errors against 128-cell nonlinear references
  are respectively (0.01451, 0.007185, 0.002555) and
  (0.02268, 0.01249, 0.004687) on 16, 32 and 64 cells. Both also approach their
  analytic collision-free mode-occupation profiles under refinement.
* Eight actual SiC runs use the same FC2/FC3 hashes as :doc:`tdbte_validation`,
  a 3x3x3 q mesh (162 modes and 833,976 shell events), eight spatial cells and
  0.02 ps duration. Two TTG amplitudes (1 and 90 K), a 30 K Gaussian pulse and
  20% selective optical excitation were each run with and without nonlinear
  collisions. Energy drift is at most 1.36e-16 of total energy; branch-budget
  residual is at most 1.70e-16 of the maximum local energy density. Maximum
  calculated Cartesian speed component is 7.69 nm/ps. These short coarse-grid
  checks do not establish physical SiC transport convergence.

The native SiC degenerate velocity off-diagonal diagnostic is 0.527 nm/ps in
the transport direction. This is recorded, not discarded or repaired, and
precludes treating the scalar native-basis calculation as a validation of
degenerate-mode coherence dynamics.

The material HDF5 files and YAML inputs are retained under
``calculations/3C-SiC/spatial-tdbte-material-final-20261107/`` and the full
synthetic convergence report under
``calculations/3C-SiC/spatial-tdbte-numerical-complete-final-20261107/``. The opt-in tool
can reproduce both with a fresh output directory and the recorded kernel.

Implementation files
------------------------

* New solver modules: ``tdbte_spatial.py``, ``tdbte_streaming.py`` and
  ``tdbte_observables.py`` in ``src/nepkappa/``.
* Existing collision integration: ``tdbte.py``, ``tdbte_builder.py``,
  ``tdbte_storage.py`` and ``tdbte_accumulate.py``. Homogeneous integration
  settings and the event flux formula are retained.
* Existing application integration: ``config.py``, ``config_schema.py``,
  ``config_models.py``, ``application.py``, ``provenance.py``, ``report.py``
  and ``tdbte_plot.py``.
* New public example: ``examples/tdbte-spatial.yaml``. The example catalog,
  input reference, TD-BTE guide, development guide and repository README are updated.
* New verification: ``tests/public/test_tdbte_spatial.py`` and
  ``tools/validate_spatial_tdbte.py``. The public configuration catalog test and
  the local example-catalog test include the new example. The earlier collision
  audit tool now makes its off-shell acoustic check optional.
