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
         kind: ttg
         delta-temperature-K: 10

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

Other initial conditions
------------------------

Spatial TD-BTE also supports thermal Gaussian pulses, selective phonon
excitation and uniform-temperature initial states.

.. list-table:: Initial conditions
   :header-rows: 1
   :widths: 25 75

   * - Type
     - Physical meaning
   * - ``ttg``
     - Spatially periodic temperature perturbation.
   * - ``gaussian``
     - Localized Gaussian temperature pulse.
   * - ``selective``
     - Spatially selective population excitation of specified branches.
   * - ``uniform``
     - Spatially uniform Bose distribution for numerical comparisons.

Thermal Gaussian pulses
~~~~~~~~~~~~~~~~~~~~~~~~

Replace ``spatial.initial`` with:

.. code-block:: yaml

   initial:
     kind: gaussian
     delta-temperature-K: 30
     sigma-nm: 8
     center-nm: 32

This sets a pulse centered at 32 nm in a periodic 64 nm domain, allowing
temperature and heat-flux propagation to be studied.

Selective phonon excitation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Excite only selected branches, such as the three optical branches of 3C-SiC:

.. code-block:: yaml

   initial:
     kind: selective
     branches: [3, 4, 5]
     amplitude: 0.2
     sigma-nm: 8
     center-nm: 32

``amplitude: 0.2`` increases the selected reference populations by up to 20%
according to the spatial envelope; it is not a 20 K temperature increase.
The calculation includes both nonequilibrium energy exchange and real-space
transport. Compared with homogeneous selective excitation, it also shows how
energy spreads from a localized region.

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
