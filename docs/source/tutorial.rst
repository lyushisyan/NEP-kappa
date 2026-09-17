Workflow recipes
==================

Start with :doc:`starting` for the first calculation. This page explains how
to adapt that workflow and reuse results. Exact keys and units are in
:doc:`input_files`; complete templates are listed in :doc:`examples`.

Run stages or reuse force constants
-------------------------------------

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
text outputs are needed; ``shengbte`` alone is not sufficient for ``kappa``.

Cached displacement forces are reused only when their fingerprints match the
structure, calculator settings, model/POTCAR, software, and workflow inputs.
Inspect ``job-input.json`` and ``job-input.sha256`` to audit reuse. Reusing
forces is distinct from restarting an interrupted external solver.

Change the transport approximation
------------------------------------

In a copy of ``examples/nep-rta.yaml``, edit the existing ``kappa`` section:

.. code-block:: yaml

   kappa:
     method: lbte
     mesh: [21, 21, 21]
     temps: [300]

Use ``method: rta`` for RTA. For Wigner transport, start from
``examples/wigner.yaml``. ``temps`` is either ``[T]`` or
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
     - ``pair-cutoff-fc3``
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

These options have different meanings. A displaced-pair cutoff is not a
universal radius applied to every force-constant tensor. Converge the chosen
cutoff together with the supercell. HiPhive requires enough independent
rattled structures; see ``examples/nep-hiphive.yaml``.

Use another calculator
------------------------

Start from ``vasp-rta.yaml`` for VASP and configure the executable, MPI command,
and POTCAR source. ``potcar_path`` can be a combined file or a library directory;
entries are assembled in POSCAR species order. Check stage-specific INCAR
settings and converge the electronic parameters before generating forces.

``mace-rta.yaml`` uses CPU inference with a user-supplied checkpoint.
Install the optional MACE dependency and use ``calculator.device: cuda`` only
when the selected environment supports it. External ASE calculators can be
configured with ``factory`` and ``kwargs``; see :doc:`input_files`.

For films, use ``film.yaml`` with a physical effective thickness. For wires,
set ``dimensionality: 1`` and the physical ``effective_area``. Meshes and
supercells must respect non-periodic directions. Geometry normalization affects
plotted transport values, not the original phono3py HDF5.

Three- and four-phonon conductivity
-------------------------------------

``examples/bas-three-four-phonon.yaml`` requests both channels:

.. code-block:: yaml

   workflow:
     preset: custom
     steps: [fc2fc3, kappa, fc4, kappa4]

The ``kappa`` stage computes pure three-phonon transport. ``kappa4`` computes
combined three-plus-four-phonon transport. The ``four-phonon`` preset runs the
combined route but does not independently add a pure-three-phonon reference.

Thirdorder/Fourthorder and FourPhonon executables are additional prerequisites.
For existing ShengBTE-format IFCs, use ``nepkappa kappa4 input.yaml``;
for FC4 alone, use ``nepkappa fc4 examples/vasp-fc4.yaml`` after configuring VASP.
``fourphonon.yaml`` demonstrates Slurm transport and sampling options.
Auto-generated CONTROL files disable non-analytic corrections; a calculation
requiring dielectric and Born-charge data needs a matching custom CONTROL.

QHA, SSCHA, and their combination
-----------------------------------

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - Route
     - Template
     - Physical change
   * - QHA
     - ``qha.yaml``
     - Thermal expansion from a volume-dependent harmonic free energy.
   * - SSCHA
     - ``sscha.yaml``
     - Temperature-renormalized FC2 at a fixed cell.
   * - QHA+SSCHA
     - ``qha-sscha.yaml``
     - SSCHA at each temperature's QHA equilibrium volume.

Run any of these templates with ``nepkappa run examples/<name>.yaml``.
The coupled route regenerates matching force constants at each volume; it does
not add two independently calculated frequency shifts.

QHA needs at least five increasing volume ratios, an equilibrium-volume range
that brackets the fit, and dynamically sensible phonons. It writes
``qha-summary.yaml``, energy-volume and temperature-property tables, and plots.
SSCHA requires converged supercells, snapshot counts, iterations, transient
iterations, and sampling mesh. Temperature-specific results are stored in the
configured ``phonopy-sscha`` work directory.

The public ``scph`` command is the Phonopy stochastic SSCHA route. Fixed-volume
transport uses temperature-dependent FC2 with the configured FC3. In the
coupled workflow, set the independent transport switches:

.. code-block:: yaml

   qha-sscha:
     three-phonon: true
     four-phonon: false

Enable the second switch to generate FC4 and run FourPhonon at each volume.
Both can be enabled, or both disabled for phonon-only work. SSCHA temperatures
must lie within the available QHA temperature range. Nested force-job arrays
are unsupported for this coupled route; use an outer batch allocation.
Use ``nepkappa report input.yaml`` to summarize completed stages.

Slurm execution
-----------------

Choose ``slurm-vasp.yaml``, ``slurm-lbte.yaml``, or ``fourphonon.yaml`` for
force arrays, distributed LBTE, or FourPhonon transport respectively.
Configure allocation sizes, MPI/OpenMP settings, environment setup, and a
shared result directory. All templates set ``submit: false``.

``validate`` and ``info`` only inspect configuration. Executing a stage with
``submit: false`` may perform preparation before writing scripts. A complete
workflow may also run earlier local stages. After reviewing generated scripts,
set ``submit: true`` to submit and inspect with ``nepkappa status input.yaml``.
Force-array collection jobs continue the remaining workflow when supported.

Compare models and test convergence
-------------------------------------

In ``examples/compare.yaml``, keep two or more labeled result directories.
Each must provide matching phonon metadata, FC2, and transport outputs. Then:

.. code-block:: bash

   nepkappa compare examples/compare.yaml

For a mesh sweep, copy ``examples/converge-qmesh.yaml`` and set its ``base``.
``base`` and ``study.directory`` resolve relative to the study YAML; paths
inside the base workflow still resolve from the launch directory.
The swept dotted key must already exist in the base input.

.. code-block:: bash

   nepkappa converge examples/converge-qmesh.yaml

With ``study.execute: false`` this writes case inputs and analysis without
running simulations. After inspection, enable execution to run the cases.
Once complete, set it back to false and repeat for the CSV, figure, and summary.
The last configured case is the numerical reference, not the exact physical
answer; incomplete cases prevent a convergence declaration.
