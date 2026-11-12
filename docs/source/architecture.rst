Calculation stages and architecture
===================================

Scientific routes
-----------------

NEP-kappa has two routes sharing force-constant preparation:

.. code-block:: text

   structure + calculator -> displacements -> forces -> FC2 / FC3 / FC4
                                                         |
                   +-------------------------------------+----------------+
                   |                                                      |
      Phonon properties and thermal transport                   Nonlinear TD-BTE
      FC2 -> dispersion, DOS, heat capacity                      FC2 + FC3 -> kernel
      FC2 + FC3 -> RTA / LBTE / Wigner                            -> n(q, branch, t)
      optional FC4 -> four-phonon transport                       optional spatial grid
      optional QHA / SSCHA -> temperature-dependent cases         -> n(x, q, branch, t)

The routes are independent: a calculated conductivity is not a TD-BTE input.
QHA and SSCHA own their physical volume, temperature and iteration loops.
Their options are not automatically transferred to TD-BTE. The spatial route
uses a one-dimensional periodic real-space grid and a three-dimensional phonon
q grid. See :doc:`tdbte_spatial` for numerical order and physical limitations.

Execution and dependency ownership
----------------------------------

.. code-block:: text

   cli.py: commands, terminal/log streams, command provenance
      -> config.py: YAML -> expanded settings -> argparse -> WorkflowConfig
           config_schema.py: accepted keys, strict YAML, scalar conversion
           config_plan.py: scientific switches and ordered steps
           config_validation.py: operation-specific prerequisites
           config_models.py: mutable configuration and typed section snapshots
      -> application.py: stage selection, timing, deferred execution
           stages/: structure relaxation and force constants
           transport.py / fourphonon.py / tdbte.py: selected solver
           qha.py / sscha.py / qha_transport.py / qha_sscha.py: scientific loops
      -> plots and reports

``command_registry.py`` defines public commands and input-selected operations.
``application.py`` supplies a runner factory to coupled temperature workflows;
those scientific modules do not import the application layer. The factory
creates a stage runner for each case while sharing the execution record.
Programmatic users executing those cases must supply ``runner_factory`` too.

``workflow.py`` is the shared context for structures, calculator instances,
force jobs and force-constant exports. It does not run plans. Stages and Slurm
workers use that context; ``calculators.py`` owns ASE-compatible backend
construction and ``adapters/vasp.py`` owns VASP I/O and execution.
``stages/base.py`` provides only prerequisite validation and execution.

Low-level identities live in ``identity.py``; provenance and calculator code
both depend on it. Temperature grids and directory names live in
``temperatures.py``. Phono3py command sampling and conventional result names
live in ``command_sampling.py``. These helpers do not import scientific
workflows. ``transport_data.py`` owns HDF5 interpretation, unit/geometry
normalization and scattering-channel selection; ``plot.py`` consumes those
results for presentation. Some old import locations remain available through
imports for internal Python callers.

Configuration and execution contracts
-------------------------------------

Argparse remains the sole source of scalar coercions, defaults and choices.
Its argument declarations are grouped by scientific section in ``config.py``.
YAML duplicate keys and unknown options are rejected before execution.
``WorkflowConfig`` is the canonical mutable state; ``WorkflowSections`` is an
uncached typed snapshot. Case workflows copy the configuration before changing
volume, temperature or output paths.

``run`` requires the six sections for its selected scientific route. Direct
commands require only their own sections. TD-BTE artifact reuse uses the
``tdbte`` section and does not prepare new forces. Keep backend-specific
operations out of CLI handlers and preserve command-specific preflight checks.

``CommandExecution.status`` distinguishes scientific completion from successful
preparation or submission: ``complete``, ``prepared`` and ``submitted``.
A deferred stage stops the plan. ``scheduler.py`` submits dependent jobs;
``slurm.py`` prepares scripts and manifests; ``run_state.py`` records their
lifecycle. Nested QHA/SSCHA transport requires synchronous execution and rejects
unexpected deferred results before continuing to later temperatures.

Force caches and reproducibility
--------------------------------

``force_identity.py`` defines schema-3 force inputs and an explicit dependency
list. The key includes the displaced structure, calculator parameters and
model/POTCAR contents, selected force-evaluation implementations, and installed
force-provider Python/native-file identities. Normalized Python syntax ignores
comments and docstrings. Workflow force methods are selected explicitly;
plotting, fitting and transport implementations are excluded from that selection.
Review the dependency list when adding a force-affecting method or backend.

``ForceArtifactStore`` also includes exact ASE arrays, cell and periodicity.
Slurm task preparation writes a checksummed ``atoms.json`` so workers restore
the exact staged state instead of losing arrays or precision through POSCAR.
Arbitrary atom info/constraints disable persistent reuse because their hidden
inputs cannot be certified by this contract.

Persistent reuse is currently certified for NEP, the ASE LJ/EMT/Morse factories,
and directly identified native VASP commands. Opaque plugins, MACE factories
and VASP launcher/wrapper commands reuse only within one calculation or Slurm
pipeline. A fresh scope token prevents reuse across independent runs; workers
and the collector share that token, so collection does not repeat worker forces. Extending certification requires an explicit provider,
model and implementation dependency contract. Installed provider files must
remain immutable for a process lifetime: their hashes are memoized within the
process. Model input contents are checked on subsequent fingerprint calls.
This is conservative dependency tracking, not proof of numerical equivalence.

``artifacts.py`` validates force shapes/finiteness and verifies the input JSON,
input digest and output SHA-256 before reuse. Both ``forces.npy`` and
Thirdorder/Fourthorder ``vasprun.xml`` are covered. Publication withdraws the
validity marker before replacing output, installs the record, and publishes the
marker last. Old records are misses; they are never migrated into valid new
records. Array outputs are replaced only after validation. Order-job backends
may overwrite XML during evaluation, so their validity marker is withdrawn
before invoking the backend. This supports single-writer interruption recovery,
not concurrent writers or
power-loss durability. Use a fresh result directory for scientific comparisons.

``provenance.py`` independently records broad software versions, including
SciPy, input identities, command history and artifact hashes. These records
remain broader than force-cache dependencies. Original command provenance is
not retroactively changed when a scheduled job finishes.

TD-BTE modules
--------------

``tdbte_builder.py`` prepares native modes/vertices and group velocities;
``tdbte_quadrature.py`` integrates energy shells. ``tdbte_storage.py`` validates
and streams versioned, checksummed event blocks. ``tdbte_accumulate.py`` supplies
NumPy reference and optional Numba reductions, without a dense collision matrix.
``tdbte.py`` handles homogeneous evolution and shared artifact loading.

``tdbte_streaming.py`` owns conservative periodic advection;
``tdbte_observables.py`` owns temperature, energy, heat flux and entropy;
``tdbte_spatial.py`` combines nonlinear collisions and streaming, audits the
result and saves HDF5 snapshots. Numerical conservation checks and independent
physical-rate/convergence validation are separate acceptance criteria.

Testing and maintenance
-----------------------

Public tests include synthetic numerical baselines, energy-shell moments,
transport normalization, cache corruption/interruption, scientific stage
selection, worker/collector agreement and deferred execution contracts.
CI explicitly exercises NumPy and Numba collision paths. Local research
snapshots and material calculations remain separate; see :doc:`development`
and :doc:`release_baseline`.

Prefer small modules with one concrete responsibility. Add shared helpers only
when multiple real consumers need them; avoid duplicate dispatchers, generic
plugin frameworks for stages, or new configuration conversion systems.
