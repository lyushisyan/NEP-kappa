Quick start: three-phonon RTA and Wigner transport
======================================================

This example uses the bundled 3C-SiC NEP model to relax the structure, generate
FC2/FC3, and calculate three-phonon RTA and SMM19 Wigner thermal transport.

Prepare the example
-----------------------

From the installed repository root, copy the structure and potential into a
working directory:

.. code-block:: bash

   conda activate nepkappa
   mkdir -p ~/work/3C-SiC
   cp examples/structures/3C-SiC/POSCAR_primitive ~/work/3C-SiC/POSCAR
   cp potentials/3C-SiC/nep_3C-SiC.txt ~/work/3C-SiC/
   cd ~/work/3C-SiC

Paths in YAML are relative to this launch directory. Later chapters use
``phonons/``, ``transport/`` and ``tdbte/`` subdirectories to keep results separate.

Prepare the input
---------------------

Save :download:`quickstart.yaml <inputs/quickstart.yaml>` as ``input.yaml`` in
this directory:

.. literalinclude:: inputs/quickstart.yaml
   :language: yaml

The example uses 3×3×3 force-constant supercells, a 31×31×31 phonon mesh and
100–1000 K in 50 K steps. ``wigner: true`` enables the particle and coherence
contributions; ``isotope: true`` includes isotope scattering.

Run the calculation
-----------------------

.. code-block:: bash

   nepkappa info input.yaml --for run
   nepkappa run input.yaml

The screen shows relaxation, force-constant generation and transport progress,
followed by completion and result paths. To show detailed external-program
output, use ``nepkappa run input.yaml -v`` instead. Both modes save detailed
output in ``results/run.log``.

For a long workstation run, use the :ref:`tmux workflow <transport-tmux>`.
On a cluster, use an allocation or a :ref:`Slurm batch script <transport-slurm>`;
replace its final command with ``nepkappa run input.yaml`` and configure
resources for this phono3py calculation.

Inspect results
-------------------

.. code-block:: bash

   tail -n 20 results/run.log
   ls -lh results/fc2.hdf5 results/fc3.hdf5 results/phono3py_disp.yaml
   ls -lh results/kappa-m313131.hdf5

``POSCAR_relaxed`` contains the optimized structure. The HDF5 force constants
and matching ``phono3py_disp.yaml`` can be reused in later chapters.
``kappa-m313131.hdf5`` stores conductivity and mode data for the specified mesh;
its filename changes when the mesh changes.

Plot and summarize
----------------------

.. code-block:: bash

   nepkappa info input.yaml --for plot
   nepkappa plot input.yaml
   nepkappa report input.yaml
   ls results/plots

With ``layout: both`` and ``tau: total``, this example produces six separate
plots and ``combined.png``. The combined figure contains dispersion, DOS,
volumetric constant-volume heat capacity, group velocity, scattering rate and
thermal conductivity. The conductivity panel distinguishes particle, coherence
and total Wigner contributions. ``plot.temperature`` selects the scattering-rate
data; the conductivity panel retains all calculated temperatures.

The report is saved as ``results/report.md`` and ``results/report.yaml``.
Plotting and reporting read saved results without recalculating forces or
transport. Further examples are listed in :doc:`examples`.
