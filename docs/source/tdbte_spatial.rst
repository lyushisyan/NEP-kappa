TD-BTE: nonequilibrium phonon transport
=======================================

The previous chapter studied energy exchange with homogeneous TD-BTE.
This chapter adds one real-space coordinate, coupling phonon propagation
to nonlinear three-phonon collisions.

NEP-kappa uses a three-dimensional q mesh and a periodic one-dimensional
real-space mesh to simulate temperature, heat flux and phonon populations
in space and time. A transient thermal grating (TTG) in 3C-SiC illustrates
the input, method and analysis.

Calculation method
------------------------

Spatial TD-BTE contains both drift and collisions:

.. math::

   \frac{\partial n_\lambda(x,t)}{\partial t}
   +
   v_{\lambda,x}
   \frac{\partial n_\lambda(x,t)}{\partial x}
   =
   \mathcal C_\lambda^{3\mathrm{ph}}[n(x,t)].

Here:

- :math:`n_\lambda(x,t)` is the phonon population in space and time.
- :math:`v_{\lambda,x}` is the group-velocity component along transport.
- :math:`\mathcal C_\lambda^{3\mathrm{ph}}` is the nonlinear three-phonon collision integral.

The added drift term describes energy carried by phonons between spatial positions.

.. note::

   One-dimensional refers to the real-space transport coordinate. Phonons
   still use the full three-dimensional q mesh. This differs from the
   ``mesh: [N, 1, 1]`` collision model for one-dimensional nanowires.

Transient thermal gratings
--------------------------

A TTG starts with a periodic temperature distribution:

.. math::

   T(x,0)
   =
   T_0+\Delta T
   \cos\left(\frac{2\pi x}{L}\right).

:math:`T_0` is the background temperature, :math:`\Delta T` the initial
modulation amplitude and :math:`L` the spatial period. Initially, each cell
has a Bose–Einstein distribution at its local temperature. The system then
evolves without continuing external heating as propagation and collisions
change the initial temperature profile.

Prepare the calculation
------------------------

Continue in the previous chapter's working directory:

.. code-block:: bash

   conda activate nepkappa
   cd ~/work/3C-SiC

Check FC2, FC3 and phonon metadata:

.. code-block:: bash

   ls results/fc2.hdf5
   ls results/fc3.hdf5
   ls results/phono3py_disp.yaml

Spatial TD-BTE uses these files for frequencies, group velocities and the
collision kernel. As in the homogeneous model, no prior conductivity
calculation is required.

Configure the input
------------------------

Create ``tdbte-spatial.yaml``:

.. code-block:: yaml

   tdbte:
     force-constants: results
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
         profile: ttg
         modes: all
         excitation: temperature
         delta: 10

   plot:
     layout: both
     dpi: 150

   output:
     result-dir: results/tdbte-spatial

The example uses a three-dimensional 3×3×3 q mesh and a periodic 64 nm
real-space domain along x. Sixteen cells give a grid spacing of:

.. math::

   \Delta x=\frac{L}{N_x}=4\ \mathrm{nm}.

The initial temperature is a cosine modulation of amplitude 10 K around 300 K.
The small q mesh and short evolution time demonstrate the spatial workflow.

Main spatial parameters
~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: Spatial TD-BTE parameters
   :header-rows: 1
   :widths: 36 64

   * - Parameter
     - Description
   * - ``length-nm``
     - Length of the one-dimensional periodic domain, in nm.
   * - ``cells``
     - Number of real-space cells.
   * - ``direction``
     - Cartesian vector defining the transport direction.
   * - ``boundary``
     - Spatial boundary condition; currently ``periodic``.
   * - ``cfl``
     - CFL number controlling drift-integration stability.
   * - ``threads``
     - CPU threads for cell-parallel collisions (default 1; requires Numba above 1).
   * - ``batch-size``
     - Number of spatial cells per collision-integral batch.
   * - ``initial``
     - Spatial distribution of the initial temperature or populations.

``length-nm`` and ``cells`` determine real-space resolution, whereas ``mesh``
determines the q mesh and phonon-mode resolution for collisions.

Run the calculation
------------------------

Check the input and run spatial TD-BTE:

.. code-block:: bash

   nepkappa info tdbte-spatial.yaml --for tdbte
   nepkappa tdbte tdbte-spatial.yaml

The program reads force constants and group velocities, builds the collision
kernel and initializes the temperature profile. It then propagates phonons
between cells while evaluating nonlinear collisions inside each cell.
Results are saved in ``results/tdbte-spatial``.

Multicore collision evaluation
-------------------------------

Spatial collisions can run on multiple CPU threads. Set these keys under
``tdbte.spatial``:

.. code-block:: yaml

   threads: 8
   batch-size: 16

Each worker processes independent spatial cells with the same shared event
chunk. Event accumulation within a cell retains its serial order, avoiding
concurrent deposition into the same mode. The default ``threads: 1`` keeps
the serial backend. More than one thread requires Numba; unavailable backends
or requests above the Numba thread limit raise an error. The requested thread
count is recorded in ``audit.json`` under ``spatial_options``.

For a single-node Slurm job, request ``--ntasks=1 --cpus-per-task=8`` and set
``NUMBA_NUM_THREADS=8`` before Python starts. Keep BLAS threads at one to avoid
nested threading; ``OMP_NUM_THREADS`` alone does not enable this feature.
``batch-size`` limits the cells processed concurrently: with 16 cells,
``threads: 8`` and ``batch-size: 16`` is a useful configuration to benchmark.
More threads than cells per batch cannot provide additional cell parallelism.

Time integration, drift and event-chunk loading remain serial. Scaling depends
on collision cost, memory bandwidth and storage; an eight-thread job is not
guaranteed to run eight times faster. Homogeneous TD-BTE still uses serial
single-state collision evaluation. This option does not implement MPI or
multi-node execution.

Time integration
------------------------

Spatial TD-BTE splits collision and drift processes. Each time step contains:

1. A half step of nonlinear three-phonon collisions.
2. A full step of spatial drift.
3. A half step of nonlinear three-phonon collisions.

Collisions use DOP853; drift uses a conservative first-order upwind scheme.
The step is limited by both ``max-step-ps`` and the CFL condition:

.. math::

   \max_\lambda |v_{\lambda,x}|
   \frac{\Delta t}{\Delta x}
   \leq \mathrm{cfl}\leq1.

Because drift is first-order upwind, the overall method is first-order in
time and space. Check numerical diffusion caused by the mesh and CFL settings.

Results and analysis
------------------------

After completion, generate plots and a report:

.. code-block:: bash

   nepkappa plot tdbte-spatial.yaml
   nepkappa report results/tdbte-spatial

Plots are saved in ``results/tdbte-spatial/plots``.

Temperature evolution in space and time
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Spatial TD-BTE saves populations at each position and time. The program
inverts the total phonon thermal energy in each cell to obtain an
energy-equivalent temperature:

.. math::

   E(x,t)=E_{\mathrm{BE}}[T(x,t)].

:math:`T(x,t)` is a local energy-equivalent temperature; it does not require
local equilibrium within a cell. For TTG initial conditions, follow the
periodic temperature perturbation over time. Its decay depends on propagation,
collisions and spatial length scales.

Heat flux and energy transport
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Spatial TD-BTE also calculates the heat-flux distribution. Modes contribute
according to their frequencies, group velocities and populations.
Temperature and heat-flux profiles reveal differences between nonequilibrium
transport and ordinary diffusion. Branch-resolved energy densities and heat
fluxes distinguish the contributions of individual branches.

TTG temperature harmonics
~~~~~~~~~~~~~~~~~~~~~~~~~

The initial TTG is a single cosine perturbation. Fourier analysis tracks
the evolution of spatial harmonics. For the fundamental:

.. math::

   A_1(t)
   =
   \frac{2}{L}
   \int_0^L
   [T(x,t)-\overline T(t)]
   \cos\left(\frac{2\pi x}{L}\right)\,dx.

Analyze :math:`A_1(t)` to compare modulation decay at different spatial periods.
In the diffusive limit with constant thermal diffusivity:

.. math::

   A_1(t)
   =
   A_1(0)
   \exp\left[
   -\alpha\left(\frac{2\pi}{L}\right)^2t
   \right].

Here, :math:`\alpha` is the thermal diffusivity. This provides a diffusive
reference, while spatial TD-BTE can examine departures from simple diffusion.

Independent initial-condition choices
--------------------------------------

``spatial.initial`` separates three independent choices:

.. list-table:: Initial-condition options
   :header-rows: 1
   :widths: 25 35 40

   * - Key
     - Values
     - Meaning
   * - ``profile``
     - ``ttg``, ``gaussian``, ``uniform``
     - Spatial envelope; default ``ttg``.
   * - ``modes``
     - ``all``, ``branches``, ``frequency-window``
     - Mode selection; default ``all``.
   * - ``excitation``
     - ``temperature``, ``relative``
     - Occupation prescription; default ``temperature``.

All 18 combinations are supported. ``uniform`` means spatially uniform;
``all`` means all active modes. Branch and frequency selections apply at
every q point represented in the kernel; zero-frequency rigid modes stay fixed.

Let :math:`p(x)` be the spatial envelope and :math:`S` the selected modes.
For ``ttg``, :math:`p(x)=\cos(2\pi x/L)`; for ``uniform``, :math:`p(x)=1`.
For ``gaussian``, the periodic Gaussian has width ``sigma-nm`` (required)
and center ``center-nm`` (default :math:`L/2`), normalized to one at its center.
The two prescriptions, for :math:`\lambda\in S`, are:

.. math::

   n_\lambda(x,0)=
   \begin{cases}
   n_{\mathrm{BE}}(\omega_\lambda,T_0+\Delta T p(x)),
      &\text{temperature},\\
   n_{\mathrm{BE}}(\omega_\lambda,T_0)[1+A p(x)],
      &\text{relative}.
   \end{cases}

Unselected modes remain at the background Bose occupation at :math:`T_0`.
With selected modes, ``temperature`` specifies the temperature of their
initial Bose occupations, not the energy-equivalent temperature of the whole
cell. Even ``relative`` applied to all modes generally creates a nonthermal
frequency distribution.

Both excitation prescriptions use ``delta``. With ``excitation: temperature``,
``delta`` is in K (default 1 K); with ``excitation: relative``, it is a
dimensionless fractional occupation change (default ``tdbte.excitation``, 0.01).
For example, ``delta: 0.2`` means 0.2 K or 20%, respectively.
TTG requires :math:`|\Delta T|<T_0` or :math:`|A|<1`; the negative cosine
half-period decreases occupations relative to the background. Gaussian and
uniform profiles require :math:`T_0+\min(\Delta T,0)>0` or :math:`A>-1`.

Selective TTG with temperature excitation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Replace ``spatial.initial`` with:

.. code-block:: yaml

   initial:
     profile: ttg
     modes: branches
     branches: [3, 4, 5]  # Zero-based IDs; verify for the material/kernel.
     excitation: temperature
     delta: 10

Set ``modes: all`` and omit ``branches`` for the full-mode thermal TTG.

Gaussian excitation of a frequency window
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   initial:
     profile: gaussian
     sigma-nm: 8
     center-nm: 32
     modes: frequency-window
     frequency-range-THz: [10, 15]
     excitation: relative
     delta: 0.2

The frequency window includes both endpoints, with :math:`0\le f_{min}<f_{max}`.
Only active modes in that interval are changed, by up to 20% at the pulse center.
An empty selection or a branch absent from the kernel is rejected.
Use ``profile: uniform`` and omit the Gaussian parameters for spatially
uniform excitation. Different windows with the same relative amplitude do
not generally receive the same injected energy.

Legacy inputs
~~~~~~~~~~~~~~

Existing ``kind: ttg``, ``gaussian`` and ``uniform`` inputs retain their thermal,
full-mode behavior. ``kind: selective`` retains relative excitation of the
specified branches, with a Gaussian or uniform spatial ``profile``.
Legacy ``delta-temperature-K`` (thermal kind) and ``amplitude`` (selective kind)
are translated to ``delta``. Supplying both ``delta`` and its legacy equivalent
is rejected. New combinations use ``delta`` and the independent keys above,
and omit ``kind``.
Combining legacy ``kind`` with ``modes`` or ``excitation`` is rejected.
Previously ignored ``branches``/``profile`` under a legacy thermal ``kind``
now raise an error rather than silently exciting all modes.

Main outputs
------------------------

Space- and time-resolved results are saved to:

.. code-block:: text

   results/tdbte-spatial/tdbte/spatial.hdf5

.. list-table:: Main output datasets
   :header-rows: 1
   :widths: 43 57

   * - Dataset
     - Contents
   * - ``occupation``
     - Populations at each position, time and phonon mode.
   * - ``temperature_K``
     - Local energy-equivalent temperature.
   * - ``heat_flux_eV_nm2_ps``
     - Spatially resolved phonon heat flux.
   * - ``branch_energy_eV_nm3``
     - Branch energy densities.
   * - ``branch_heat_flux_eV_nm2_ps``
     - Branch heat fluxes.
   * - ``branch_collision_power_eV_nm3_ps``
     - Branch energy-exchange rates from three-phonon collisions.
   * - ``ttg_cosine_amplitude_K``
     - Fundamental TTG cosine temperature amplitude.
   * - ``departure_from_local_BE``
     - Departure from the local equal-energy Bose distribution.

``audit.json`` records checks of energy conservation, branch energy balances,
entropy, population positivity and CFL. The program plots temperature, heat
flux, branch energies and temperature harmonics to analyze spatial
nonequilibrium transport.
