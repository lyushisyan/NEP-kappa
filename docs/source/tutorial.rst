Advanced usage
==============

The main workflows are in :doc:`harmonic_temperature`,
:doc:`thermal_transport` and :doc:`tdbte`. This page collects stage reuse,
calculator choices and execution variants. See :doc:`input_files` for parameters.

Run stages or reuse force constants
---------------------------------------

A normal three-phonon workflow is equivalent to:

.. code-block:: bash

   nepkappa relax input.yaml
   nepkappa fc2fc3 input.yaml
   nepkappa kappa input.yaml

Plot and summarize completed results separately:

.. code-block:: bash

   nepkappa plot input.yaml
   nepkappa report input.yaml

Use ``fc2`` instead of ``fc2fc3`` for harmonic force constants only.
With relaxation enabled, separately invoked FC stages expect ``POSCAR_relaxed``
in the result directory. Run ``relax`` first, or set relaxation to false when
the supplied POSCAR is already relaxed.

For a new mesh or temperature range, keep the matching ``fc2.hdf5``,
``fc3.hdf5``, and ``phono3py_disp.yaml`` together and rerun only ``kappa``.
Inspect that stage with ``nepkappa info input.yaml --for kappa``.
Choose ``force-constant.format: both`` when both phono3py HDF5 and ShengBTE
text outputs are needed; ``shengbte`` alone is not sufficient for phono3py transport.

Cached displacement forces are reused only when their fingerprints match the
structure, calculator settings, model/POTCAR, software, and workflow inputs.
Inspect ``job-input.json`` and ``job-input.sha256`` to audit reuse. Reusing
forces is distinct from restarting an interrupted external solver.

Cutoffs and force-constant fitting
--------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 25 35 40

   * - Route
     - Key
     - Meaning
   * - Native phono3py FC3
     - ``cutoff-fc3``
     - Displaced-pair distance, angstrom.
   * - Thirdorder FC3
     - ``cutoff-fc3``
     - Negative integer: neighbor shell; positive value: nm.
   * - HiPhive FC2/FC3
     - ``cutoffs``
     - Fitting cutoffs, angstrom, by order.
   * - Fourthorder FC4
     - ``cutoff-fc4``
     - Passed to Fourthorder; negative values select neighbor shells.
   * - Native FC2
     - ``dim-fc2``
     - No independent real-space cutoff; increase the represented supercell.

The phono3py pair cutoff limits displacement pairs, rather than truncating all
tensor elements at that radius. Check cutoff and supercell convergence together.
For HiPhive fitting options, see :doc:`input_files`.

Use another calculator
--------------------------

Start from ``examples/vasp-rta-3ph.yaml`` for VASP and configure the executable, MPI command,
and POTCAR source. ``potcar-path`` can be a combined file or a library directory;
entries are assembled in POSCAR species order. Check stage-specific INCAR
settings and converge the electronic parameters before generating forces.

For MACE, install the optional dependency and supply a checkpoint in a copied
six-section input. Use ``calculator.device: cuda`` only
when the selected environment supports it. External ASE calculators can be
configured with ``factory`` and ``kwargs``; see :doc:`input_files`.

For films, set a physical effective thickness in your input. For wires,
set ``dimensionality: 1`` and the physical ``effective-area``. Meshes and
supercells must respect non-periodic directions. Geometry normalization affects
plotted transport values, not the original phono3py HDF5.

Slurm execution
-------------------

The ten public inputs keep site-specific scheduler resources out of their
scientific settings. Use an external batch script for a whole workflow, or
add an optional top-level ``parallel.force-constant`` for displaced-force
arrays and ``parallel.kappa`` for distributed phono3py LBTE or FourPhonon.
Configure allocation sizes, MPI/OpenMP settings, environment setup, and a
shared result directory for the actual cluster. Start with ``submit: false``
when generating NEP-kappa-managed scripts for inspection.

``info`` only inspects configuration. Executing a stage with
``submit: false`` may perform preparation before writing scripts. A complete
workflow may also run earlier local stages. After reviewing generated scripts,
set ``submit: true`` to submit and inspect the jobs with the site's Slurm tools.
Force-array collection jobs continue the remaining workflow when supported.

Plot existing results
-------------------------

For FC2-only work, place ``fc2.hdf5`` with its matching ``phono3py_disp.yaml``,
``phonopy.yaml``, or ``phonopy_disp.yaml`` in the result directory. A minimal
plot input is:

.. code-block:: yaml

   output:
     result-dir: calculations/my-material/results
   kappa:
     mesh: [21, 21, 21]
     temps: [0, 1000, 100]
   plot:
     layout: both
     path: seekpath

From the directory relative to which those paths are defined:

.. code-block:: bash

   nepkappa info plot.yaml --for plot
   nepkappa plot plot.yaml

This produces dispersion, DOS, volume heat capacity, and group velocity,
as separate PNGs and a 2-by-2 combination. No potential or FC3 is required.
The mesh and temperatures control harmonic-property sampling and require
convergence checks. Existing same-named plot files may be replaced.

With transport HDF5 present, select its matching mesh. The scattering-rate
panel needs linewidth data; cumulative conductivity needs ``mode_kappa``.
Unsupported or incomplete transport files raise an error.
Adding a ``plot`` section alone does not generate figures; run the explicit
``nepkappa plot`` command after the necessary outputs exist.

``report`` summarizes existing results. Custom manuscript figures use separate
scripts. See :doc:`input_files` for built-in panels
and :doc:`troubleshooting` for file-selection problems.

.. _change-the-transport-approximation:
.. _three-and-four-phonon-conductivity:
.. _qha-and-qha-volume-sscha:
.. _time-dependent-bte:

Method-specific workflows
-------------------------

The numbered chapters contain the complete inputs and analysis instructions:

* :doc:`harmonic_temperature`: harmonic phonons, QHA, SSCHA and QHA+SSCHA.
* :doc:`thermal_transport`: LBTE, same-backend 3ph/4ph comparison, tmux and Wigner Slurm submission.
* :doc:`tdbte`: optical excitation, numerical audits and six-panel relaxation plots.

For conductivity at QHA volumes or with SSCHA FC2, use the transport templates
in :doc:`examples`. Set the corresponding ``run-transport`` switch and ensure
conductivity temperatures are covered by the QHA/SSCHA grids, as specified in
:doc:`input_files`.
