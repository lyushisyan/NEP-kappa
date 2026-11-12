Input files and parameters
==========================

``structure``: crystal structure
--------------------------------

``structure`` specifies the crystal structure, relaxation and geometric
normalization for low-dimensional systems.

.. list-table:: Structure parameters
   :header-rows: 1
   :widths: 38 62

   * - Parameter
     - Description
   * - ``poscar``
     - Crystal structure in VASP POSCAR format.
   * - ``relaxation``
     - Enable structure relaxation using a boolean or detailed settings.
   * - ``dimensionality``
     - System dimensionality: ``1``, ``2`` or ``3``; default ``3``.
   * - ``effective-thickness``
     - Effective thickness of a two-dimensional system, in Å.
   * - ``effective-area``
     - Effective cross-sectional area of a one-dimensional system, in Å².
   * - ``vacuum-axis``
     - Optional vacuum direction for a two-dimensional system.
   * - ``periodic-axis``
     - Optional periodic direction for a one-dimensional system.

Structure relaxation
~~~~~~~~~~~~~~~~~~~~~~~~

Set ``relaxation: true`` to enable relaxation. VASP also supports multiple stages:

.. code-block:: yaml

   structure:
     poscar: POSCAR
     relaxation:
       enabled: true
       workdir: vasp-relax
       stages:
         coarse: {ediffg: -0.05, prec: Normal}
         fine: {ediffg: -0.01, prec: Accurate}

.. note::

   ``nepkappa run`` performs relaxation and subsequent calculations in order.
   For direct ``fc2``, ``fc2fc3`` or ``fc2fc3fc4`` calculations with relaxation
   enabled, first run ``nepkappa relax`` to generate ``POSCAR_relaxed``.
   QHA reads the specified input structure and requires ``relaxation: false``.

Geometric normalization for low-dimensional systems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For films and nanowires containing vacuum, normalize the physical volume
using an effective thickness or cross-sectional area.

.. list-table:: Low-dimensional settings
   :header-rows: 1
   :widths: 20 72 38

   * - System
     - Required settings
     - Optional settings
   * - Two-dimensional film
     - ``dimensionality: 2`` and ``effective-thickness``
     - ``vacuum-axis``
   * - One-dimensional nanowire
     - ``dimensionality: 1`` and ``effective-area``
     - ``periodic-axis``

These settings normalize volumetric heat capacity and conductivity during
plotting. They do not modify the original phono3py HDF5 data.

.. note::

   QHA, SSCHA and FourPhonon currently support only three-dimensional systems
   (``dimensionality: 3``). Geometric normalization does not extend their scope.

``calculator``: interatomic interactions
----------------------------------------

``calculator`` selects the energy and atomic-force method. NEP-kappa supports
NEP, VASP, MACE and other ASE-compatible calculators.

.. list-table:: Calculator types
   :header-rows: 1
   :widths: 16 16 64

   * - Calculator
     - ``name``
     - Main settings
   * - NEP
     - ``nep``
     - Set the potential with ``model``; CPU evaluation is the default.
   * - MACE
     - ``mace``
     - Select ``model`` or ``foundation``; ``device`` and ``dtype`` control the device and numerical precision.
   * - VASP
     - ``vasp``
     - Configure ``vasp-command``, ``vasp-path`` and ``potcar-path``.
   * - ASE
     - ``ase``
     - Configure other calculators with ``factory``, ``kwargs`` and ``model-files``.

.. _nep-gpu-force-batches:

NEP GPU acceleration
~~~~~~~~~~~~~~~~~~~~~~~~

Use GPU batches to evaluate forces for many displaced structures:

.. code-block:: yaml

   calculator:
     name: nep
     model: /path/to/nep.txt
     nep-force-backend: gpu-batch
     nep-batch-size: 128
     nep-command: /path/to/nep

Main settings:

- ``nep-force-backend``: ``cpu`` or ``gpu-batch``.
- ``nep-batch-size``: maximum number of structures per batch.
- ``nep-command``: path to GPUMD's ``nep`` executable.

GPU batches support native phono3py FC2/FC3 and HiPhive displacement forces.
QHA, SSCHA, Thirdorder and Fourthorder are currently unsupported.

.. note::

   GPU batches require Calorine 3.5 or later, GPUMD's ``nep`` executable
   and a compatible GPU environment.

``force-constant``: force constants
-----------------------------------

``force-constant`` controls second-, third- and fourth-order force-constant
generation. It supports finite displacements, HiPhive fitting and different backends.

.. list-table:: Force-constant parameters
   :header-rows: 1
   :widths: 35 65

   * - Parameter
     - Description
   * - ``dim-fc2`` / ``dim-fc3`` / ``dim-fc4``
     - Supercell dimensions for FC2, FC3 and FC4, respectively.
   * - ``fc3-backend``
     - Defaults to phono3py; ``thirdorder`` is also supported.
   * - ``cutoff-fc3``
     - Third-order displacement cutoff; its meaning depends on the backend.
   * - ``cutoff-fc4``
     - Fourthorder interaction cutoff.
   * - ``use-hiphive``
     - Enable HiPhive force-constant fitting.
   * - ``format``
     - Output format: ``phono3py``, ``shengbte`` or ``both``.
   * - ``thirdorder-command`` / ``fourthorder-command``
     - Executable paths or launch commands for Thirdorder and Fourthorder.
   * - ``qha`` / ``sscha`` / ``four-phonon``
     - Enable and configure the corresponding calculation modules.

Supercells and cutoffs
~~~~~~~~~~~~~~~~~~~~~~~~

``dim-fc2``, ``dim-fc3`` and ``dim-fc4`` independently control the supercells
for each order. Supercell size affects cost and convergence; test it for your material.

``cutoff-fc3`` depends on the FC3 backend:

- **phono3py**: distance cutoff for displaced-atom pairs, in Å.
- **Thirdorder**: a negative integer selects a neighbor shell; a positive value is a distance cutoff in nm.

``cutoff-fc4`` is passed to Fourthorder; negative integers select neighbor shells.
HiPhive uses ``n-structures``, ``rattle-std``, ``cutoffs`` and ``min-dist``
to control training structures and fitting. ``cutoffs`` is in Å.

QHA and SSCHA
~~~~~~~~~~~~~~~~~~~~~~~~

NEP-kappa supports QHA, fixed-volume SSCHA and SSCHA at QHA equilibrium volumes.

.. list-table:: Temperature-dependent phonons
   :header-rows: 1
   :widths: 25 75

   * - Method
     - Calculations
   * - QHA
     - Volume scans and equation-of-state fits give temperature-dependent equilibrium volumes and thermal-expansion effects.
   * - SSCHA
     - Finite-temperature sampling at fixed volume self-consistently updates auxiliary harmonic force constants.
   * - QHA + SSCHA
     - QHA determines equilibrium volumes before fixed-volume SSCHA is run at each volume.

Main settings:

- ``qha.volume-ratios``: volume scan range.
- ``qha.temps``: QHA temperature range.
- ``sscha.temps``: SSCHA sampling temperatures.
- ``sscha.snapshots``: sampled structures per iteration.
- ``sscha.iterations``: self-consistency iterations.
- ``run-transport``: calculate conductivity at the corresponding temperatures.

See :doc:`temperature` for methods and complete examples.

``kappa``: phonon scattering and conductivity
---------------------------------------------

``kappa`` configures phonon scattering, BTE solvers and lattice conductivity.
NEP-kappa supports phono3py and FourPhonon, with different scattering mechanisms
and solver combinations.

Common parameters
~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: Thermal-transport parameters
   :header-rows: 1
   :widths: 36 64

   * - Parameter
     - Description
   * - ``engine``
     - Calculation engine: ``phono3py`` or ``fourphonon``.
   * - ``mesh``
     - Brillouin-zone q mesh, for example ``[31, 31, 31]``.
   * - ``temps``
     - Calculation temperatures, in K.
   * - ``method``
     - phono3py solver: ``rta`` or ``lbte``.
   * - ``method-3ph`` / ``method-4ph``
     - Three- and four-phonon solvers in FourPhonon.
   * - ``isotope``
     - Include isotope scattering.
   * - ``bfmp``
     - Boundary-scattering mean free path, in μm.
   * - ``wigner``
     - Enable Wigner coherence transport.
   * - ``force-constants``
     - Existing force-constant source directory; ordinary direct phono3py transport only.
   * - ``command``
     - Custom phono3py launch command.

Engines and solvers
~~~~~~~~~~~~~~~~~~~~~~~~

phono3py primarily handles three-phonon transport; FourPhonon includes
three- and four-phonon scattering. Both support RTA and iterative LBTE,
with different allowed combinations.

.. list-table:: Supported transport methods
   :header-rows: 1
   :widths: 25 30 45

   * - Engine
     - Solver
     - Calculations
   * - phono3py
     - ``rta``
     - Three-phonon relaxation-time approximation.
   * - phono3py
     - ``lbte``
     - Iterative solution of the three-phonon linearized BTE.
   * - FourPhonon
     - ``rta / rta``
     - RTA for both three- and four-phonon scattering.
   * - FourPhonon
     - ``lbte / rta``
     - Iterative three-phonon solution with four-phonon RTA.
   * - FourPhonon
     - ``lbte / lbte``
     - Iterative solution for both three- and four-phonon scattering.

RTA approximates collisions using a relaxation time for each mode; iterative
LBTE includes linear coupling between modes. FourPhonon selects the two
solvers using ``method-3ph`` and ``method-4ph``. ``rta / lbte`` is unsupported.

Wigner coherence transport
~~~~~~~~~~~~~~~~~~~~~~~~~~

``wigner: true`` enables Wigner transport, separating conductivity into
particle and coherence contributions:

.. math::

   \kappa_{\mathrm{total}}
   = \kappa_{\mathrm{P}} + \kappa_{\mathrm{C}}.

The particle term describes transport from phonon populations; the coherence
term describes off-diagonal contributions between vibrational modes.
Wigner is not a separate BTE solver; its supported scope depends on the engine:

- **phono3py**: three-phonon RTA or LBTE, with SMM19 coherence contributions.
- **FourPhonon**: requires Wigner_Park; currently RTA for both three- and four-phonon scattering.

See :doc:`thermal_transport` for inputs and examples.

Calculation temperatures
~~~~~~~~~~~~~~~~~~~~~~~~

``temps`` accepts a single temperature or an equally spaced range:

.. list-table:: Temperature settings
   :header-rows: 1
   :widths: 35 65

   * - Format
     - Example
   * - Single temperature ``[T]``
     - ``temps: [300]``
   * - Range ``[Tmin, Tmax, step]``
     - ``temps: [100, 1000, 50]``

Temperatures are in K and must be positive. Ranges include both endpoints;
the positive step must land exactly on the upper endpoint. Arbitrary lists
of unequally spaced temperatures are not supported.

For coupled QHA or SSCHA transport, conductivity temperatures must also meet
the sampling-range requirements described in :doc:`temperature`.

Reusing force constants
~~~~~~~~~~~~~~~~~~~~~~~~

With FC2/FC3 already generated, skip atomic-force calculations and use phono3py
to calculate conductivity at different temperatures or q meshes. For example:

.. code-block:: yaml

   kappa:
     force-constants: calculations/3C-SiC/fc
     engine: phono3py
     mesh: [41, 41, 41]
     temps: [300]
     method: rta

   output:
     result-dir: calculations/3C-SiC/kappa-m41

Check the input and run:

.. code-block:: bash

   nepkappa info kappa.yaml --for kappa
   nepkappa kappa kappa.yaml

``force-constants`` selects the source directory and ``output.result-dir``
the new result directory. The source must contain matching files:

- ``fc2.hdf5``: second-order force constants.
- ``fc3.hdf5``: third-order force constants.
- ``phono3py_disp.yaml``: structure, displacements and associated metadata.

When source and output directories differ, required inputs are copied to the
output directory and the original files are preserved.

Custom phono3py commands
~~~~~~~~~~~~~~~~~~~~~~~~

Advanced users can set ``kappa.command`` to a complete phono3py launch command
instead of the automatically generated command line.
Its mesh and temperatures must agree with the YAML settings. When a custom
command is used, NEP-kappa does not append other calculation arguments.

The command runs directly in the result directory without a shell, so pipes
and redirection are unsupported. See :doc:`thermal_transport` for further examples.

``tdbte``: nonequilibrium phonon dynamics
-----------------------------------------

``tdbte`` solves the time-dependent BTE for the evolution of nonequilibrium
phonon populations and energy exchange between modes. It does not require a
prior steady-state transport solution, but needs matching FC2, FC3 and phonon
structure metadata.

Calculation routes
~~~~~~~~~~~~~~~~~~~~~~~~

Two routes are supported:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Route
     - Purpose
   * - Homogeneous TD-BTE
     - Time evolution of nonequilibrium populations, energy exchange and thermalization.
   * - One-dimensional spatial TD-BTE
     - Phonon propagation and nonlinear three-phonon collisions in periodic one-dimensional real space.

Both routes share the three-phonon collision kernel. The optional
``tdbte.spatial`` settings enable spatial transport; the phonon q mesh remains
three-dimensional.

Common parameters
~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: TD-BTE parameters
   :header-rows: 1
   :widths: 35 65

   * - Parameter
     - Description
   * - ``force-constants``
     - Directory containing existing FC2/FC3 and matching metadata.
   * - ``kernel``
     - Existing collision kernel to reuse without rebuilding.
   * - ``mesh``
     - Phonon q mesh used for the collision integral.
   * - ``temperature``
     - Reference temperature for the initial Bose distribution.
   * - ``branches``
     - Excited phonon branches, indexed from zero.
   * - ``excitation``
     - Relative increase in the selected modes' populations.
   * - ``duration-ps``
     - Total evolution time, in ps.
   * - ``max-step-ps``
     - Maximum integration step, in ps.
   * - ``samples``
     - Number of saved time points.
   * - ``spatial``
     - Length, mesh, boundary and initial conditions for one-dimensional spatial transport.

In homogeneous TD-BTE, ``excitation`` multiplies the selected branches' initial
populations by ``1 + excitation``. It is not a temperature increase or laser fluence.

.. note::

   Current TD-BTE excludes four-phonon scattering, Wigner coherence and electron coupling.

See :doc:`tdbte` for homogeneous dynamics and :doc:`tdbte_spatial` for spatial transport.

``plot``: plotting results
--------------------------

``plot`` controls visualization of phonon properties, scattering, conductivity
and nonequilibrium dynamics. It reads saved data without repeating force or
thermal-transport calculations.

.. list-table:: Plot parameters
   :header-rows: 1
   :widths: 35 65

   * - Parameter
     - Description
   * - ``layout``
     - ``separate``, ``combined`` or ``both`` for individual and combined plots.
   * - ``path``
     - High-symmetry dispersion path; ``seekpath`` generates it automatically.
   * - ``path-points`` / ``path-segments``
     - Custom dispersion path.
   * - ``temperature``
     - Scattering-rate plot temperature; matching saved data are selected by default.
   * - ``kappa``
     - Conductivity component: ``x``, ``y``, ``z`` or ``all``.
   * - ``tau``
     - Scattering mechanisms to plot.
   * - ``scattering-rate-min``
     - Lower scattering-rate display limit, in ps⁻¹.
   * - ``dpi``
     - PNG image resolution.
   * - ``tdbte-times``
     - Three snapshot times for homogeneous TD-BTE dispersion plots.

Plot contents and required data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Different plots require different calculation outputs:

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Plot
     - Required data
   * - Dispersion and DOS
     - FC2 and matching structure metadata.
   * - Volumetric heat capacity and group velocities
     - FC2 and corresponding phonon data.
   * - Phonon scattering rates
     - Phonon linewidths or scattering-rate data.
   * - Lattice thermal conductivity
     - Completed transport results.
   * - Nonequilibrium dynamics
     - Saved TD-BTE trajectories.

``layout: both`` generates individual and combined plots.

Scattering-rate plots
~~~~~~~~~~~~~~~~~~~~~~~~

``tau`` supports:

- ``phph``: phonon–phonon scattering.
- ``isotope``: isotope scattering.
- ``total``: phonon–phonon scattering, isotope scattering and their sum.
- ``normal`` / ``umklapp``: N / U scattering.
- ``nu``: show both N and U scattering.
- ``all``: all available scattering information.

phono3py linewidths are converted to scattering rates using ``4*pi*gamma``.
``total`` excludes boundary scattering. ``nu`` requires both N and U data.
For FourPhonon, three- and four-phonon rates are plotted separately rather than
automatically summed.

Conductivity plots
~~~~~~~~~~~~~~~~~~~~~~~~

Ordinary LBTE results containing RTA data allow both solutions to be compared
in one plot. Wigner results show particle, coherence and total conductivity.
``kappa: all`` selects the mean of the three diagonal components, not all
components of the conductivity tensor.

``output``: result output
-------------------------

``output`` configures the result directory and progress display.

.. list-table:: Output parameters
   :header-rows: 1
   :widths: 35 65

   * - Parameter
     - Description
   * - ``result-dir``
     - Directory for calculation outputs.
   * - ``progress``
     - Enable dynamic progress output.

Use separate result directories for different materials, temperatures or q
meshes to organize, compare and reuse data. The log is saved as ``run.log``
inside the result directory. After completion, ``nepkappa plot`` and
``nepkappa report`` generate plots and reports from the saved results.

.. _parallel-execution:

``parallel``: parallel execution
--------------------------------

The optional ``parallel`` section configures resources for force-constant and
transport stages. NEP-kappa supports Slurm scheduling with separate resources
for each stage. For example, distribute force jobs and phono3py transport:

.. code-block:: yaml

   parallel:
     force-constant:
       backend: slurm
       jobs: 16
       submit: false

     kappa:
       backend: slurm
       jobs: 4
       max-concurrent: 4
       cpus-per-task: 8
       collect-cpus-per-task: 8
       submit: false

.. list-table:: Common Slurm settings
   :header-rows: 1
   :widths: 35 65

   * - Parameter
     - Description
   * - ``backend``
     - Execution backend, for example ``slurm``.
   * - ``jobs``
     - Number of worker tasks.
   * - ``max-concurrent``
     - Maximum number of simultaneously running tasks.
   * - ``cpus-per-task``
     - CPU cores per task.
   * - ``collect-cpus-per-task``
     - CPU cores for the collection stage.
   * - ``submit``
     - Automatically submit Slurm jobs.
   * - ``extra-sbatch``
     - Additional Slurm resource options.

``submit: false`` generates scripts without calling ``sbatch``, but preparation
may still perform local work. Both phono3py RTA and LBTE support distributing
irreducible q points across multiple Slurm jobs.
