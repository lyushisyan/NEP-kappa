Harmonic phonons and temperature renormalization
====================================================

This chapter covers harmonic phonons, quasiharmonic thermal expansion (QHA),
fixed-volume stochastic self-consistent harmonic approximation (SSCHA), and
SSCHA at QHA equilibrium volumes. Transport calculations follow in
:doc:`thermal_transport`.

Workflow and preparation
----------------------------

Continue from the :doc:`starting` working directory:

.. code-block:: bash

   mkdir -p ~/work/3C-SiC/phonons
   cd ~/work/3C-SiC/phonons
   cp ../POSCAR .
   cp ../nep_3C-SiC.txt .

Download the input for each section into this directory. All four routes use
separate result directories and the same NEP model.

.. list-table::
   :header-rows: 1
   :widths: 22 50 28

   * - Route
     - Calculation
     - Command
   * - Harmonic
     - Relaxation followed by FC2 generation.
     - ``relax``, then ``fc2``
   * - QHA
     - Volume scan, free-energy fit and equilibrium-volume phonons.
     - ``fc2``
   * - Fixed-volume SSCHA
     - Initial FC2, thermal sampling and self-consistent FC2 updates.
     - ``run``
   * - QHA+SSCHA
     - QHA equilibrium volumes followed by SSCHA at each temperature.
     - ``run``

Harmonic phonon properties
------------------------------

Input and execution
~~~~~~~~~~~~~~~~~~~~~~~

Save :download:`harmonic.yaml <inputs/harmonic.yaml>`:

.. literalinclude:: inputs/harmonic.yaml
   :language: yaml

.. code-block:: bash

   nepkappa info harmonic.yaml --for relax
   nepkappa relax harmonic.yaml
   nepkappa info harmonic.yaml --for fc2
   nepkappa fc2 harmonic.yaml

The separate ``fc2`` command reads the optimized structure because relaxation
is enabled. The outputs include ``results/harmonic/POSCAR_relaxed``, FC2 and
matching phonon metadata. No third-order force constants are needed here.

Plot and interpret
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   nepkappa info harmonic.yaml --for plot
   nepkappa plot harmonic.yaml
   ls results/harmonic/plots

FC2-only plotting produces four panels: dispersion, DOS, volumetric heat
capacity and group velocity. ``layout: both`` writes four separate PNGs and a
2×2 ``combined.png``. The mesh controls Brillouin-zone integration; the band
path controls the displayed dispersion.

QHA: thermal expansion and phonons
--------------------------------------

Input and execution
~~~~~~~~~~~~~~~~~~~~~~~

Save :download:`qha.yaml <inputs/qha.yaml>`. It starts from the optimized
harmonic structure and scans five volume ratios:

.. literalinclude:: inputs/qha.yaml
   :language: yaml

.. code-block:: bash

   nepkappa info qha.yaml --for fc2
   nepkappa fc2 qha.yaml

``equilibrium-phonons: true`` additionally calculates FC2 at each fitted
thermal-equilibrium volume. ``run-transport: false`` keeps this example focused
on phonon and thermodynamic properties. The scan must bracket the equilibrium
volumes over the requested temperature range.

Plot and interpret
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   nepkappa info qha.yaml --for plot
   nepkappa plot qha.yaml
   ls results/qha/plots

The four-panel overview contains free energy versus volume, equilibrium volume,
volumetric thermal expansion, and temperature-dependent dispersion with DOS.
PNG and PDF outputs include ``qha_combined`` and the four individual figures.
A scan point with volume ratio 1 supplies the static-volume phonon reference.
QHA frequency changes arise from volume changes.

SSCHA: fixed-volume phonon renormalization
----------------------------------------------

Input and execution
~~~~~~~~~~~~~~~~~~~~~~~

Save :download:`sscha.yaml <inputs/sscha.yaml>`:

.. literalinclude:: inputs/sscha.yaml
   :language: yaml

.. code-block:: bash

   nepkappa info sscha.yaml --for run
   nepkappa run sscha.yaml

This route uses Phonopy's stochastic SSCHA implementation with NEP energies
and forces at fixed volume. It samples 300, 600 and 900 K. Converge the snapshot
count, iteration count, supercell and sampling mesh for the material.

Plot and interpret
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   nepkappa info sscha.yaml --for plot
   nepkappa plot sscha.yaml
   ls results/sscha/plots

With transport disabled, the four-panel overview overlays all requested SSCHA
temperatures in dispersion, DOS and group velocity. The heat-capacity panel
shows a point evaluated from each temperature's FC2, alongside the fixed-volume
harmonic reference. ``plot.temperature`` does not select a single temperature
for this overview.

The exported FC2 is an auxiliary harmonic matrix, not the free-energy Hessian.
Its heat capacity is a harmonic-model estimate. Optional FC3 bubble shifts are
stored separately and do not modify this FC2 or the conductivity.

QHA+SSCHA: expansion and renormalization
--------------------------------------------

Input and execution
~~~~~~~~~~~~~~~~~~~~~~~

Save :download:`qha-sscha.yaml <inputs/qha-sscha.yaml>`:

.. literalinclude:: inputs/qha-sscha.yaml
   :language: yaml

.. code-block:: bash

   nepkappa info qha-sscha.yaml --for run
   nepkappa run qha-sscha.yaml

QHA first determines the equilibrium volumes. SSCHA then updates FC2 at the
corresponding temperature and volume, keeping that volume fixed during its
iterations. The SSCHA temperatures must lie within the QHA range.

Plot and interpret
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   nepkappa info qha-sscha.yaml --for plot
   nepkappa plot qha-sscha.yaml
   ls results/qha-sscha/plots

The four-panel overview compares SSCHA spectra with harmonic references at
matching QHA volumes. The heat-capacity reference instead uses the fixed-volume
QHA scan point with volume ratio 1; SSCHA values are separate colored points.
This sequence combines QHA expansion with SSCHA renormalization without
minimizing the volume using SSCHA free energy.

For temperature-dependent conductivity, see the coupled transport settings in
:doc:`input_files` and :doc:`tutorial`.
