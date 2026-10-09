Calculation guide
==================

The commands below assume a configured YAML input. See :doc:`starting` for
the 3C-SiC examples and :doc:`input_files` for parameter definitions.

Run stages or reuse force constants
-------------------------------------

A normal three-phonon workflow is equivalent to:

.. code-block:: bash

   nepkappa stage relax input.yaml
   nepkappa stage fc2fc3 input.yaml
   nepkappa stage kappa input.yaml

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
text outputs are needed; ``shengbte`` alone is not sufficient for ``kappa``.

Cached displacement forces are reused only when their fingerprints match the
structure, calculator settings, model/POTCAR, software, and workflow inputs.
Inspect ``job-input.json`` and ``job-input.sha256`` to audit reuse. Reusing
forces is distinct from restarting an interrupted external solver.

Change the transport approximation
------------------------------------

In a copy of ``examples/nep-rta-wigner-3ph.yaml``, edit the existing ``kappa``
section. Set ``wigner: false`` if only the diagonal LBTE result is needed:

.. code-block:: yaml

   kappa:
     method: lbte
     mesh: [21, 21, 21]
     temps: [300]

Use ``method: rta`` for RTA. The two retained three-phonon Wigner examples
show both methods. ``temps`` is either ``[T]`` or
``[minimum, maximum, step]``; three entries are not an arbitrary temperature list.
Do not append a second section with the same YAML key.

Cutoffs and force-constant fitting
------------------------------------

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
------------------------

Start from ``examples/vasp-rta-3ph.yaml`` for VASP and configure the executable, MPI command,
and POTCAR source. ``potcar_path`` can be a combined file or a library directory;
entries are assembled in POSCAR species order. Check stage-specific INCAR
settings and converge the electronic parameters before generating forces.

For MACE, install the optional dependency and supply a checkpoint in a copied
six-section input. Use ``calculator.device: cuda`` only
when the selected environment supports it. External ASE calculators can be
configured with ``factory`` and ``kwargs``; see :doc:`input_files`.

For films, set a physical effective thickness in your input. For wires,
set ``dimensionality: 1`` and the physical ``effective_area``. Meshes and
supercells must respect non-periodic directions. Geometry normalization affects
plotted transport values, not the original phono3py HDF5.

Three- and four-phonon conductivity
-------------------------------------

``examples/nep-rta-3ph-4ph.yaml`` computes combined three-plus-four-phonon
conductivity. If a separate pure-three-phonon reference is also needed, use
matching FC2/FC3 in a three-phonon ``kappa`` stage. A custom stage plan can
request both channels:

.. code-block:: yaml

   workflow:
     preset: custom
     steps: [fc2fc3, kappa, fc4, kappa4]

The ``kappa`` stage computes pure three-phonon transport. ``kappa4`` computes
combined three-plus-four-phonon transport. The ``four-phonon`` preset runs the
combined route but does not independently add a pure-three-phonon reference.

Thirdorder/Fourthorder and FourPhonon executables are additional prerequisites.
For existing ShengBTE-format IFCs, use ``nepkappa stage kappa4 input.yaml``;
for FC4 alone, use ``nepkappa stage fc4 input.yaml`` after configuring the input.
The three retained 3ph+4ph templates demonstrate RTA/RTA, LBTE/RTA, and
LBTE/LBTE choices.
Auto-generated CONTROL files disable non-analytic corrections; a calculation
requiring dielectric and Born-charge data needs a matching custom CONTROL.

QHA and QHA-volume SSCHA
-----------------------------------

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - Route
     - Template
     - Physical change
   * - QHA-volume 3ph RTA
     - ``nep-qha-rta-3ph.yaml``
     - New FC2/FC3 and conductivity at each equilibrium volume.
   * - SSCHA at QHA volumes
     - ``nep-qha-sscha-rta-3ph.yaml``
     - SSCHA FC2 at each temperature's QHA equilibrium volume.

Run a configured template with ``nepkappa run examples/<name>.yaml``.
The coupled route regenerates force constants at each QHA volume, which stays
fixed during SSCHA. Exported FC2 is the auxiliary harmonic matrix.
``scph.bubble: true`` writes separate input-FC3 diagonal on-shell shifts;
``nepkappa stage bubble input.yaml`` applies this step to existing SSCHA results.
Neither operation updates the conductivity with bubble corrections.
See :doc:`input_files` for the method's limits.

QHA needs at least five increasing volume ratios, an equilibrium-volume range
that brackets the fit, and dynamically sensible phonons. It writes
``qha-summary.yaml``, energy-volume and temperature-property tables, and plots.
SSCHA requires converged supercells, snapshot counts, iterations, transient
iterations, and sampling mesh. Temperature-specific results are stored in the
configured ``phonopy-sscha`` work directory.

The public ``scph`` command remains available for custom fixed-volume SSCHA
inputs, but this catalog does not include a standalone SSCHA example. In the
coupled six-section input, enable the nested QHA and SSCHA settings:

.. code-block:: yaml

   force-constant:
     qha:
       enabled: true
     sscha:
       enabled: true
       run-transport: true
     four-phonon: false

Enabling ``four-phonon`` additionally generates FC4 and runs FourPhonon at
each volume. SSCHA temperatures
must lie within the available QHA temperature range. Nested force-job arrays
are unsupported for this coupled route; use an outer batch allocation.
Use ``nepkappa report input.yaml`` to summarize completed stages.

Slurm execution
-----------------

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
set ``submit: true`` to submit and inspect with ``nepkappa status input.yaml``.
Force-array collection jobs continue the remaining workflow when supported.

Plot existing results
-----------------------------------------------------

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

With transport HDF5 present, select its matching mesh. Scattering/lifetime
panels need linewidth data; cumulative conductivity needs ``mode_kappa``.
Unsupported or incomplete transport files raise an error.
Adding a ``plot`` section alone does not add a plotting stage
to a workflow: use the explicit command above or a custom stage plan.

``report`` summarizes existing results. Custom manuscript figures use separate
scripts. See :doc:`input_files` for built-in panels
and :doc:`troubleshooting` for file-selection problems.

Time-dependent BTE
--------------------------------

Use ``examples/tdbte.yaml`` with matching FC2/FC3 and phono3py metadata. Set
``tdbte.force-constants``, ``mesh`` and the excitation, then validate with
``--for tdbte``. The stage builds a kernel and solves it in disk-backed chunks;
physical relaxation-rate validation remains separate from numerical audits.
See :doc:`tdbte` for excitation definitions, numerical audits, and outputs.
