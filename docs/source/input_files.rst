Input Files
=============

NEP-kappa has two calculation input shapes. **Static** calculations use six
top-level sections: ``structure``, ``calculator``, ``force-constant``,
``kappa``, ``plot``, and ``output``. A standard NEP three-phonon input is:

.. code-block:: yaml

   structure:
     poscar: examples/structures/Si/POSCAR_bulk
     relaxation: true
   calculator:
     name: nep
     nep_model: potentials/Si/Si_Bulk_Fan.txt
   force-constant:
     dim-fc2: [3, 3, 3]
     dim-fc3: [3, 3, 3]
     qha: false
     sscha: false
     four-phonon: false
   kappa:
     engine: phono3py
     mesh: [21, 21, 21]
     temps: [100, 1000, 50]
     method: rta
   plot:
     layout: both
     path: seekpath
   output:
     result_dir: calculations/example-runs/nep-rta

The switches in ``force-constant`` select the static route automatically:

.. list-table::
   :header-rows: 1
   :widths: 39 61

   * - Enabled switches
     - ``nepkappa run`` stages
   * - None
     - Structure preparation, FC2/FC3, three-phonon conductivity
   * - ``qha``
     - QHA volume scan
   * - ``qha`` with ``kappa.qha-volumes: true``
     - QHA volume scan, then FC2/FC3 and phono3py RTA conductivity at each target temperature's equilibrium volume
   * - ``sscha``
     - FC2/FC3 and fixed-volume SSCHA; transport only with ``run-transport: true``
   * - ``four-phonon``
     - Structure preparation, FC2/FC3, FC4, and four-phonon conductivity
   * - ``qha`` and ``sscha``
     - QHA followed by SSCHA at QHA volumes
   * - All three
     - QHA-volume SSCHA with FC4 and four-phonon transport at each temperature

Each switch may be ``false``, ``true``, or a mapping with ``enabled`` and
method-specific settings. For example:

.. code-block:: yaml

   force-constant:
     dim-fc2: [3, 3, 3]
     dim-fc3: [3, 3, 3]
     qha:
       enabled: true
       volume-ratios: [0.94, 0.97, 1.00, 1.03, 1.06]
       temps: [0, 900, 10]
     sscha:
       enabled: true
       temps: [300, 900, 300]
       snapshots: 200
       iterations: 4
       run-transport: true
     four-phonon:
       enabled: false

Set ``kappa.engine: phono3py`` for three-phonon transport and
``kappa.engine: fourphonon`` when ``force-constant.four-phonon`` is enabled.
Older inputs without ``engine`` infer the same choice from the four-phonon
switch. When ``four-phonon`` is enabled, its mapping also accepts FC4 options
such as ``dim-fc4`` and FourPhonon executable options such as ``command``.
The default command is ``ShengBTE_cpu``; point it to the actual
FourPhonon executable for your installation. A four-phonon switch executes
both FC4 generation and conductivity. Standalone QHA plus four-phonon and
standalone SSCHA plus four-phonon are rejected because those combinations do
not have a supported corrected-FC2 transport path. The QHA+SSCHA route is a
sequential QHA-volume approximation, not free-energy volume optimization.
With ``parallel.kappa.backend: slurm`` and ``submit: false``, the FourPhonon stage
prepares a job script for inspection instead of submitting the calculation.
``qha`` and ``sscha`` are grouped under ``force-constant`` for input
organization, but each launches its full temperature-dependent workflow.
For QHA-coupled RTA transport, set ``kappa.qha-volumes: true``. NEP-kappa then
recomputes FC2 and FC3 at every requested QHA equilibrium volume and runs
phono3py at the corresponding temperature. This accounts for isotropic thermal
expansion but does not itself add explicit anharmonic frequency renormalization.
The requested ``kappa.temps`` must lie inside the QHA temperature range.
Per-temperature results are saved under ``qha-kappa/TxxxxK/``;
``nepkappa plot`` uses ``plot.temperature`` to select one case.

Parallel settings are optional and separate from the six scientific sections:

.. code-block:: yaml

   parallel:
     force-constant:
       backend: slurm
       jobs: 16
       submit: false
     kappa:
       backend: slurm
       submit: false

``parallel.force-constant`` controls displaced-structure force jobs.
``parallel.kappa`` controls phono3py LBTE or FourPhonon jobs according to
``kappa.engine``. For FourPhonon it can also hold ``mpi-processes``,
``mpi-launcher``, ``omp-threads``, and ``omp-stacksize``. Older nested
parallel locations remain accepted; do not specify the same setting in both
places. QHA-coupled calculations run their per-temperature force stages inside
one allocation and do not support nested Slurm arrays.

The ``structure.relaxation`` spelling controls optimization within the structure
section; the older ``relaxation.enabled`` section remains accepted. If both
are present, they must agree. QHA reads ``structure.poscar`` directly, so the
static form rejects ``structure.relaxation: true`` when QHA is enabled.
The ``plot`` section only configures plotting:
run ``nepkappa plot input.yaml`` separately to create figures. Optional
parameters use parser defaults; supercell, q mesh, and temperatures are
research choices that require convergence checks.

**Dynamic** TD-BTE calculations use a separate, shorter input:

.. code-block:: yaml

   tdbte:
     force-constants: calculations/SiC/fc
     mesh: [5, 5, 5]
     temperature: 300
     branches: [3, 4, 5]
   output:
     result-dir: calculations/SiC/tdbte

With ``tdbte`` and ``output`` alone, ``nepkappa run`` selects the dynamic
route. It reuses existing matching FC2/FC3 and metadata; it does not need a
calculator or regenerate force constants. Do not mix static and dynamic
sections in one new input. See :doc:`tdbte` for the physical model and its
limits.

The earlier ``workflow.stages`` interface remains available for advanced and
existing inputs. It
describes a calculation in six phases: structure, forces, force constants,
temperature, transport, and analysis. Each phase selects a method and can
enable supported features. Detailed numerical and external-program settings
remain in the corresponding sections:

- ``workflow``: a stage-first plan, an older preset, or an advanced step list
- ``structure``: input POSCAR
- ``calculator``: NEP, VASP, a generic ASE calculator, or an installed plugin
- ``relaxation``: structure relaxation settings
- ``force-constant``: force-constant generation settings
- ``qha``: isotropic volume scan and quasi-harmonic settings
- ``scph``: Phonopy stochastic SSCHA settings
- ``qha-sscha``: independent three/four-phonon transport switches for the coupled route
- ``kappa``: thermal-conductivity settings passed to ``phono3py``
- ``fourphonon``: combined three-plus-four-phonon transport settings
- ``plot``: plot layout, band path, relaxation-time channel, and kappa component
- ``tdbte``: dynamics from FC2/FC3 or an existing energy-shell kernel; see :doc:`tdbte`
- ``output``: progress display and result directory

YAML validation is strict. Unknown sections and keys are rejected before a
calculation starts, and likely misspellings include a suggested supported key.
Free-form VASP ``vasp_kwargs`` and relaxation-stage mappings remain available
for native VASP options.

Command behavior
------------------

Run a stage-first input with:

.. code-block:: bash

   nepkappa run input.yaml

``nepkappa init input.yaml`` generates the six-section static form for all
five supported initializer goals. Existing preset and stage-first inputs remain
valid for compatibility.

``workflow.stages`` is a mapping in fixed calculation order. An omitted or
disabled phase does not run. ``structure.method: input`` uses the input
structure; ``relax`` executes optimization. ``forces.method`` selects the
calculator. ``force-constants`` selects a generation method and ordered force
constant orders. ``temperature`` can be omitted or choose QHA, SSCHA, or
QHA-volume SSCHA. ``transport`` chooses one solver route, and ``analysis``
enables optional result operations.

.. code-block:: yaml

   workflow:
     stages:
       structure: {method: input}
       forces: {method: nep}
       force-constants:
         method: finite-displacement
         orders: [2, 3]
         features: {compact: true}
       temperature: {enabled: false}
       transport:
         method: three-phonon-rta
         features: {isotope: false}
       analysis:
         method: plot
         features: {report: true}

The nine public static examples use the six-section form rather than this
compatibility form. See ``examples/nep-rta-wigner-3ph-4ph.yaml`` for the
four-phonon Wigner route. Supported stage-first choices are:

.. list-table::
   :header-rows: 1
   :widths: 22 49 29

   * - Stage
     - Methods
     - Features
   * - ``structure``
     - ``input``, ``relax``
     - Relaxation details in ``relaxation``
   * - ``forces``
     - NEP, VASP, MACE, or an installed ASE calculator/plugin name
     - Calculator details in ``calculator``
   * - ``force-constants``
     - ``finite-displacement``, ``hiphive``, ``thirdorder``; ``orders`` is
       ``[2]``, ``[2, 3]``, or ``[2, 3, 4]``
     - ``compact``; numerical settings in ``force-constant``
   * - ``temperature``
     - ``qha``, ``sscha``, ``qha-sscha``
     - ``bubble``, ``three-phonon``, ``four-phonon`` where supported
   * - ``transport``
     - ``three-phonon-rta``, ``three-phonon-lbte``,
       ``three-phonon-wigner``, ``four-phonon-rta``,
       ``four-phonon-wigner``
     - ``isotope``; ``boundary-mfp`` for three-phonon or ``nonanalytic``
       for four-phonon
   * - ``analysis``
     - ``tdbte``, ``plot``, or ``report``
     - Enable additional operations with ``features: {report: true}``, etc.

The stage compiler rejects conflicting choices between the stage map and
detailed sections. Temperature workflows and a separate transport stage cannot
be chained in this interface: QHA does not automatically update conductivity,
and standalone SSCHA transport runs within its own stage. For QHA-volume
SSCHA, enable its three/four-phonon transport under
``temperature.features``. The older ``workflow.preset`` choices
(``three-phonon``, ``four-phonon``, ``qha``, ``scph``, ``qha-sscha``) and
``workflow.preset: custom`` with ``workflow.steps`` remain accepted.

When ``force-constants.orders`` includes 4, the default export format becomes
``both`` so the FourPhonon inputs are available. After a Slurm FourPhonon job
finishes, any following analysis stages resume from the completed transport
stage. Inspect the generated job script before submission.

The following commands expose individual stages for advanced use:

.. code-block:: bash

   nepkappa stage relax input.yaml
   nepkappa stage fc2 input.yaml
   nepkappa stage fc2fc3 input.yaml
   nepkappa stage fc4 input.yaml
   nepkappa stage qha input.yaml
   nepkappa stage scph input.yaml
   nepkappa stage bubble input.yaml
   nepkappa stage qha-sscha input.yaml
   nepkappa stage tdbte input.yaml
   nepkappa stage kappa input.yaml
   nepkappa stage kappa4 input.yaml
   nepkappa plot input.yaml
   nepkappa info input.yaml
   nepkappa report calculations/runs/calculation

- ``nepkappa stage relax`` relaxes the structure and writes ``POSCAR_relaxed`` to ``output.result_dir``.
- ``nepkappa stage fc2`` generates ``phono3py_disp.yaml`` and ``fc2.hdf5`` only.
- ``nepkappa stage fc2fc3`` generates ``phono3py_disp.yaml``, ``fc2.hdf5``, and ``fc3.hdf5``.
- ``nepkappa stage fc4`` generates FourPhonon ``FORCE_CONSTANTS_4TH`` using ``Fourthorder_vasp.py``.
- ``nepkappa stage qha`` computes isotropic quasi-harmonic thermal properties over a volume scan.
- ``nepkappa stage scph`` runs Phonopy stochastic SSCHA.
- ``nepkappa stage qha-sscha`` runs SSCHA at volumes interpolated from a completed QHA fit.
- ``nepkappa stage kappa`` computes thermal conductivity using existing ``phono3py_disp.yaml``, ``fc2.hdf5``, and ``fc3.hdf5``.
- ``nepkappa stage kappa4`` runs FourPhonon with existing ShengBTE-format force constants.
- ``nepkappa plot`` creates harmonic plots from FC2 and matching metadata, adding transport panels when compatible conductivity data exist.
- ``nepkappa stage bubble`` postprocesses completed SSCHA results with diagonal on-shell frequency shifts; it does not update transport.
- ``nepkappa stage tdbte`` builds a kernel from matching force constants and a q mesh,
  or reuses a kernel, then evolves populations in chunks. Numerical audits do
  not independently validate physical relaxation rates.
- ``nepkappa report`` writes ``report.yaml`` and ``report.md`` from an existing result tree.
- ``nepkappa run`` expands and executes the selected workflow preset. Without a
  ``workflow`` section it preserves the legacy ``relax`` + ``fc2fc3`` +
  ``kappa`` behavior.
- ``nepkappa info`` prints the parsed configuration without running a calculation.

``nepkappa stage fc2fc3`` computes and writes FC2 first, then starts the FC3
displacement, force, and export stage.

When ``relaxation.enabled`` is ``true``, ``nepkappa stage fc2`` and
``nepkappa stage fc2fc3`` read
``POSCAR_relaxed`` from ``output.result_dir``. Run ``nepkappa stage relax`` first, or
use ``nepkappa run``.

Configuration validation is command-aware. For example, ``nepkappa stage scph``
validates its calculator, structure, and ``scph`` sections, while
``nepkappa stage kappa`` validates transport settings without requiring QHA or
force-generation settings to be complete.

Finite-displacement FC2/FC3 and HiPhive force evaluations are resumable. NEP
and external ASE force jobs are stored in
``output.result_dir/force-jobs/<stage>/<index>``; VASP jobs remain in
``output.result_dir/vasp-runs/<stage>/<index>`` by default. Each completed job
records ``POSCAR``, ``forces.npy``, ``job-input.sha256``, and
``job-input.json``. A rerun reuses only a force array whose audited input
fingerprint still matches, and recalculates missing, invalid, or stale entries.

Repository-provided examples
------------------------------

See :doc:`examples` for the maintained workflow catalog. Public structures
are in ``examples/structures/`` and reusable models are in ``potentials/``.
Run examples from the repository root. Private cluster settings and research
results belong in the ignored ``calculations/`` directory.

``fourphonon``
----------------

``nepkappa stage kappa4`` stages ShengBTE-format ``FORCE_CONSTANTS_2ND``,
``FORCE_CONSTANTS_3RD``, and ``FORCE_CONSTANTS_4TH`` and runs FourPhonon.
Set ``harmonic-format: espresso`` when the harmonic input is an official
``espresso.ifc2`` file; this mode requires a matching custom ``CONTROL``.
Solver choices are ``rta``, ``3ph-iterative``, and ``full-iterative``.
In six-section static inputs, set both ``kappa.method-3ph`` and
``kappa.method-4ph`` to ``rta`` or ``lbte``. NEP-kappa maps ``rta/rta`` to
``rta``, ``lbte/rta`` to ``3ph-iterative``, and ``lbte/lbte`` to
``full-iterative``. FourPhonon has no ``rta/lbte`` mode. Do not also set
``kappa.method`` or a conflicting ``force-constant.four-phonon.solver``.
This mapping follows the `FourPhonon manual
<https://github.com/FourPhonon/FourPhonon/blob/main/Manual.md>`_.
Sampling settings expose FourPhonon's 3ph/4ph scattering and phase-space
process counts. MPI+OpenMP and Slurm resources are configured in the same
section. Results are normalized to ``kappa4-rta.dat`` and
``kappa4-iterative.dat`` and summarized in ``fourphonon-summary.yaml``.

Set ``fourphonon.wigner: true`` to use an executable built from FourPhonon's
`Wigner_Park branch <https://github.com/FourPhonon/FourPhonon/tree/Wigner_Park>`_.
This branch adds the four-phonon scattering rate to the three-phonon and isotope
rates before evaluating the RTA population conductivity and the off-diagonal
Wigner coherence term. It writes separate population, coherence, and total
tensors; NEP-kappa checks that the total equals the sum of the first two and
normalizes them to ``kappa4-rta.dat``, ``kappa4-rta-coherence.dat``, and
``kappa4-rta-total.dat``. ``nepkappa report`` lists all three components.
In the branch source, the RTA scattering rate and the coherence resonance
factor have the forms

.. math::

   \Gamma_s = \Gamma_{s,3\mathrm{ph}} + \Gamma_{s,\mathrm{iso}}
   + \Gamma_{s,4\mathrm{ph}}, \qquad
   L_{ss'} = \frac{\Gamma_s+\Gamma_{s'}}
   {4(\omega_s-\omega_{s'})^2+(\Gamma_s+\Gamma_{s'})^2}.

The branch reconstructs ``rate`` from the RTA response vector when computing
coherence; a mode with exactly zero diagonal group velocity receives zero
reconstructed rate. This implementation detail should be considered when
comparing very flat modes or other Wigner implementations.

This mode requires ``solver: rta``, ``harmonic-format: shengbte``, and all four
``sample-*`` settings set to ``-1``. The Wigner_Park branch does not support
the newer FourPhonon sampling flags. Point ``fourphonon.command`` to the
Wigner_Park executable, not a regular FourPhonon build. NEP-kappa invokes that
external executable; it does not implement the Wigner formula itself. The
FourPhonon Wigner_Park calculation is distinct from phono3py's experimental
``kappa.wigner`` (SMM19) path. See ``examples/nep-rta-wigner-3ph-4ph.yaml`` for a
complete input. Converge the force constants, q mesh, and scattering settings
before interpreting material results.

Automatically generated CONTROL files disable non-analytic corrections. Polar
crystals such as SiC should set ``fourphonon.control`` to a validated CONTROL
containing dielectric and Born effective-charge data.

Example outputs are written to ``calculations/example-runs/...`` and are not
tracked by Git.

Minimal NEP example
---------------------

.. code-block:: yaml

   structure:
     poscar: examples/structures/Si/POSCAR_bulk
     dimensionality: 3

   calculator:
     name: nep
     nep_model: potentials/Si/Si_Bulk_Fan.txt

   relaxation:
     enabled: true

   force-constant:
     dim-fc2: [3, 3, 3]
     dim-fc3: [3, 3, 3]
     fc3-backend: phono3py
     format: phono3py
     use_hiphive: false
     compact-fc: true

   kappa:
     mesh: [21, 21, 21]
     temps: [100, 1000, 50]
     method: rta
     isotope: false
     bfmp: 1.0e6
     wigner: false

   plot:
     layout: both
     path: seekpath
     tau: total
     kappa: all
     temperature: 300
     dpi: 300

   output:
     progress: true
     result_dir: calculations/example-runs/nep-rta

HiPhive example
-----------------

.. code-block:: yaml

   structure:
     poscar: examples/structures/Si/POSCAR_film
     dimensionality: 2
     effective_thickness: 10.0
     vacuum_axis: z

   calculator:
     name: nep
     nep_model: potentials/Si/Si_NWs_XuKe.txt

   relaxation:
     enabled: false

   force-constant:
     dim-fc2: [4, 4, 1]
     dim-fc3: [4, 4, 1]
     use_hiphive: true
     n_structures: 500
     rattle_std: 0.03
     min_dist: 2.2
     cutoffs: [5.0, 4.0]

   kappa:
     mesh: [21, 21, 1]
     temps: [100, 1000, 50]
     method: rta
     isotope: false
     bfmp: 1.0e6
     wigner: false

   output:
     progress: true
     result_dir: calculations/example-runs/film-hiphive

VASP example
--------------

.. code-block:: yaml

   structure:
     poscar: examples/structures/Si/POSCAR_bulk

   calculator:
     name: vasp
     vasp_command: env -u DISPLAY -u XAUTHORITY mpirun -np 24 /path/to/vasp_std
     vasp_path: /path/to/vasp_std
     potcar_path: /path/to/potpaw_PBE.64

   relaxation:
     enabled: true
     workdir: vasp-relax
     stages:
       coarse:
         nsw: 80
         ibrion: 2
         isif: 3
         ediff: 1.0e-5
         ediffg: -0.05
         prec: Normal
       fine:
         nsw: 150
         ibrion: 2
         isif: 3
         ediff: 1.0e-6
         ediffg: -0.01
         prec: Accurate

   force-constant:
     workdir: vasp-runs
     vasp_kwargs:
       encut: 650
       kspacing: 0.2
       kgamma: true
       kpar: 2
     ncore: 3
     lreal: false
     dim-fc2: [3, 3, 3]
     dim-fc3: [3, 3, 3]
     use_hiphive: false
     compact-fc: true

   kappa:
     mesh: [21, 21, 21]
     temps: [100, 1000, 50]
     method: rta
     isotope: false
     bfmp: 1.0e6
     wigner: false

   output:
     progress: true
     result_dir: calculations/example-runs/vasp-rta

``structure``
---------------

``poscar`` is the input POSCAR path. Relative paths are interpreted from the
directory where ``nepkappa`` is launched, usually the repository root.

.. code-block:: yaml

   structure:
     poscar: examples/structures/Si/POSCAR_bulk
     dimensionality: 3

``dimensionality`` controls effective-volume corrections used by
``nepkappa plot``:

- ``3``: bulk material. No effective-geometry correction is applied.
- ``2``: film or slab. Set ``effective_thickness`` in Angstrom and optionally
  ``vacuum_axis`` (``x``, ``y``, or ``z``; default ``z``).
- ``1``: nanowire. Set ``effective_area`` in Angstrom^2 and optionally
  ``periodic_axis`` (``x``, ``y``, or ``z``; default ``z``).

For films, NEP-kappa scales cell-normalized volume heat capacity and thermal
conductivity by ``cell_thickness / effective_thickness``. For nanowires, it
scales them by ``cell_cross_section / effective_area``. These corrections are
applied only during plotting; the original ``kappa-m*.hdf5`` file is not
modified.

``calculator``
----------------

For NEP calculations:

.. code-block:: yaml

   calculator:
     name: nep
     nep_model: potentials/Si/Si_Bulk_Fan.txt

For VASP calculations:

.. code-block:: yaml

   calculator:
     name: vasp
     vasp_command: mpirun -np 24 /path/to/vasp_std
     vasp_path: /path/to/vasp_std
     potcar_path: /path/to/potpaw_PBE.64

``vasp_command`` is used when a full launcher command is needed. If it is not
set, NEP-kappa uses ``vasp_path``. ``potcar_path`` may point to a ready POTCAR
file or a potential-library directory. For multi-element POSCAR files, NEP-kappa
concatenates POTCAR chunks in POSCAR element order.

Changing the VASP installation does not require editing Python source. Update
``calculator.vasp_command`` in the input. If both ``vasp_command`` and
``vasp_path`` are present, the command wins: changing only ``vasp_path`` has
no effect. Using only ``vasp_path`` does not automatically launch MPI ranks.
Prefer absolute executable and POTCAR paths on the actual execution host.

The command is split into arguments, not evaluated by a shell. Put environment
setup such as ``module load`` in the batch script or force-job
``force-constant.parallel.preamble``, not in ``vasp_command``. Shell variables,
``~``, pipes and redirections are not expanded there. Force-job preambles do
not configure separately launched relaxation or QHA stages.

On the target host, ``command -v vasp_std`` checks the current PATH; no result
may mean a module has not been loaded. VASP and licensed POTCAR files are not
bundled with NEP-kappa. Parser validation does not establish executable,
MPI-library or compute-node availability. The :doc:`input_assistant` can help
inspect these separately without launching a calculation.

MACE is available as a native backend. A local checkpoint uses:

.. code-block:: yaml

   calculator:
     name: mace
     model: models/Si_MACE.model
     device: cuda
     dtype: float64

   relaxation:
     enabled: false

MACE-MP foundation models use the installed ``mace_mp`` convenience factory:

.. code-block:: yaml

   calculator:
     name: mace
     foundation: mp
     model: medium
     device: cuda
     dtype: float64

Omit ``model`` to accept the installed factory's current default. Local
checkpoints are SHA256-hashed. Automatically downloaded foundation models are
identified by factory, model name, MACE version, and configuration; use a local
checkpoint for byte-for-byte provenance. Extra constructor options can be
provided under ``calculator.kwargs``.

For another machine-learning potential, install its Python package and point to
an ASE-compatible calculator class or factory:

.. code-block:: yaml

   calculator:
     name: ase
     factory: some_package.calculators:MLCalculator
     kwargs:
       model: models/potential.model
       device: cuda
       dtype: float64
     model-files:
       - models/potential.model

   relaxation:
     enabled: false

The configured callable receives ``kwargs`` and must return an ASE-compatible
calculator. Existing file paths nested in ``kwargs`` are hashed automatically;
``model-files`` explicitly identifies additional model artifacts for cache and
provenance hashing. External calculators can evaluate finite-displacement,
HiPhive, thirdorder, and Fourthorder structures. They also support structure
relaxation through ASE when the calculator provides energy, forces, and stress
for cell optimization. Use ``relaxation.enabled: false`` for a pre-relaxed
structure.

Plugin packages may expose a short name with a Python entry point such as:

.. code-block:: toml

   [project.entry-points."nepkappa.calculators"]
   my_potential = "my_package.calculators:make_calculator"

The workflow can then set ``calculator.name: my_potential`` and omit
``factory``. The entry-point callable follows the same ``factory(**kwargs)``
contract.

``relaxation``
----------------

.. code-block:: yaml

   relaxation:
     enabled: true

For VASP relaxation, optional staged settings can be provided:

.. code-block:: yaml

   relaxation:
     enabled: true
     workdir: vasp-relax
     stages:
       coarse:
         ediffg: -0.05
         prec: Normal
       fine:
         ediffg: -0.01
         prec: Accurate

``enabled: false`` means the input POSCAR is copied to ``POSCAR_relaxed`` during
the relax stage.

``force-constant``
--------------------

Finite-displacement route:

.. code-block:: yaml

   force-constant:
     dim-fc2: [3, 3, 3]
     dim-fc3: [3, 3, 3]
     fc3-backend: phono3py
     format: phono3py
     use_hiphive: false
     compact-fc: true

``dim-fc2`` sets the phonon supercell used for FC2. ``dim-fc3`` sets the
supercell used for FC3. The deprecated ``dim`` key is still accepted for
compatibility and is interpreted as ``dim-fc2``, ``dim-fc3``, and ``dim-fc4``
when the explicit keys are not set.
``compact-fc`` defaults to ``true`` and writes phono3py v4 compact FC2/FC3
arrays. Set it to ``false`` to write full supercell force-constant arrays.
``format`` controls FC2/FC3 export files:

- ``phono3py``: write ``fc2.hdf5`` and, in ``fc2fc3`` mode, ``fc3.hdf5``
- ``shengbte``: write ``FORCE_CONSTANTS_2ND`` and, in ``fc2fc3`` mode,
  ``FORCE_CONSTANTS_3RD``
- ``both``: write both phono3py HDF5 files and ShengBTE text files

The current ``nepkappa stage kappa`` command reads phono3py HDF5 files. Use
``format: phono3py`` or ``format: both`` when the next stage is
``nepkappa stage kappa``. Use ``format: shengbte`` when the next stage is a
ShengBTE/FourPhonon workflow.

For native phono3py FC3, ``cutoff-fc3`` is an optional displaced-pair
distance cutoff in Angstrom. For ``fc3-backend: thirdorder``, the same
``cutoff-fc3`` key follows a different convention: a negative integer selects a neighbor shell and a
positive value is a distance in nm. HiPhive uses its own real-space
``cutoffs`` list in Angstrom. Native finite-displacement FC2 has no independent
real-space cutoff: ``dim-fc2`` controls the represented interaction range.
Blindly zeroing fitted FC2 elements after the calculation is not provided
because it can break translational/rotational sum rules.

The same ``cutoff-fc3`` key is used for both FC3 backends, but its units and
meaning follow the selected backend. Omitting it means no pair cutoff for
phono3py and the third neighbor shell (``-3``) for Thirdorder. The old
``pair-cutoff-fc3`` name remains a phono3py compatibility alias; do not use
both names with different values. A phono3py pair cutoff limits displacement
pairs, not all FC3 tensor elements beyond a radius.

Thirdorder FC3 route:

.. code-block:: yaml

   force-constant:
     dim-fc3: [4, 4, 4]
     fc3-backend: thirdorder
     cutoff-fc3: -3
     thirdorder-command: thirdorder_vasp.py
     fc3-workdir: fc3-thirdorder-runs
     format: shengbte

With ``fc3-backend: thirdorder``, ``nepkappa stage fc2fc3`` generates FC2 first,
runs ``thirdorder_vasp.py sow``, evaluates every ``3RD.POSCAR.*`` structure,
and passes the ordered XML list to ``thirdorder_vasp.py reap``. VASP jobs keep
their native ``vasprun.xml`` files. NEP jobs write a minimal force-only XML
with forces in eV/Angstrom, matching the parser in ``thirdorder_vasp.py``.
``format: both`` additionally converts ``FORCE_CONSTANTS_3RD`` to a full
phono3py ``fc3.hdf5``.

FourPhonon FC4 route:

.. code-block:: yaml

   force-constant:
     dim-fc4: [2, 2, 2]
     cutoff-fc4: -2
     fourthorder-command: Fourthorder_vasp.py
     fc4-workdir: fc4-runs

``nepkappa stage fc4`` runs ``Fourthorder_vasp.py sow`` in
``output.result_dir/fc4-workdir``, calculates forces for every generated
``4TH.POSCAR.*`` structure with VASP, NEP, or an external ASE/plugin calculator, and pipes the resulting XML list
into ``Fourthorder_vasp.py reap``. VASP supplies native ``vasprun.xml`` files;
NEP uses the same minimal force-only XML adapter as thirdorder. The final
``FORCE_CONSTANTS_4TH`` file is copied to ``output.result_dir``.
``dim-fc4`` defaults to ``dim-fc3``. ``cutoff-fc4`` is passed directly to
Fourthorder; negative values follow Fourthorder's neighbor-shell convention.

HiPhive route:

.. code-block:: yaml

   force-constant:
     dim-fc2: [3, 3, 3]
     dim-fc3: [3, 3, 3]
     use_hiphive: true
     n_structures: 100
     rattle_std: 0.03
     min_dist: 2.0
     cutoffs: [5.0, 4.0]

Slurm force arrays
~~~~~~~~~~~~~~~~~~

Native phono3py FC2/FC3 and HiPhive force evaluations can be split across a
Slurm array:

.. code-block:: yaml

   force-constant:
     dim-fc2: [3, 3, 3]
     dim-fc3: [3, 3, 3]
     fc3-backend: phono3py
     parallel:
       backend: slurm
       jobs: 16
       max-concurrent: 4
       partition: normal4
       nodes: 1
       ntasks: 16
       cpus-per-task: 1
       time: "7-00:00:00"
       memory: 96G

The initial command writes displaced structures and task lists under
``output.result-dir/force-slurm``, submits the array, and submits FC fitting
with an ``afterok`` dependency. ``jobs`` is the maximum number of array
elements; structures are distributed round-robin within them.
``max-concurrent`` throttles simultaneous elements. Set ``submit: false`` to
generate scripts without calling ``sbatch``. A rerun uses the audited force
cache and only recalculates missing or stale structures. The force-array route
currently supports native phono3py and HiPhive, not thirdorder FC3.

For VASP force calculations, ``workdir`` controls the subdirectory under
``output.result_dir`` and ``vasp_kwargs`` supplies important INCAR/KPOINTS
settings:

.. code-block:: yaml

   force-constant:
     workdir: vasp-runs
     vasp_kwargs:
       encut: 650
       kspacing: 0.2
       kgamma: true
       lreal: false

``qha``
---------

.. code-block:: yaml

   qha:
     volume-ratios: [0.94, 0.96, 0.98, 1.00, 1.02, 1.04, 1.06]
     dim-fc2: [3, 3, 3]
     mesh: [31, 31, 31]
     temps: [0, 1000, 10]
     eos: vinet
     pressure: 0.0
     relax-internal: true
     relax-fmax: 0.001
     relax-steps: 200
     displacement-distance: 0.01
     imaginary-frequency-tolerance: -0.1
     cutoff-frequency: 0.0001

``volume-ratios`` are ratios relative to the input cell volume and must contain
at least five positive, unique values in ascending order. ``dim-fc2`` defaults
to ``force-constant.dim-fc2``. ``mesh`` is the harmonic q-point mesh and
``temps`` is ``[minimum, maximum, step]`` in kelvin. Supported equations of
state are ``vinet``, ``birch_murnaghan``, and ``murnaghan``; ``pressure`` is in
GPa. ``relax-internal`` relaxes atomic positions while keeping each scaled cell
fixed. The imaginary-frequency tolerance is in THz.
``cutoff-frequency`` excludes numerical near-zero acoustic modes from the
thermal sums and is also expressed in THz.

For VASP, ``vasp-relax-kwargs`` and ``vasp-static-kwargs`` may contain extra
INCAR mappings for the fixed-cell relaxation and static-energy stages. Generic
ASE/plugin calculators and NEP use the same QHA workflow through ASE. QHA is
currently limited to three-dimensional periodic structures.

``scph``
----------

The workflow samples the canonical ensemble of the current
harmonic FC2, evaluates every supercell with the configured NEP/MACE/ASE
calculator, and refits FC2 using Phonopy and symfc:

.. code-block:: yaml

   scph:
     initial-fc2: calculations/example-runs/harmonic/fc2.hdf5
     born: BORN
     workdir: phonopy-sscha
     temps: [300, 900, 300]
     snapshots: 1000
     iterations: 10
     transient: 2
     sscha-mesh: [20, 20, 20]
     random-seed: 271828
     cutoff-frequency: 0.01
     fc-calculator: symfc
     save-datasets: false
     run-transport: true
     transport-fc3: calculations/example-runs/harmonic/fc3.hdf5
     transport-metadata: calculations/example-runs/harmonic/phono3py_disp.yaml

``initial-fc2`` defaults to ``output.result-dir/fc2.hdf5``. ``born``,
``transport-fc3``, and ``transport-metadata`` are optional. Each temperature is
checkpointed below ``workdir/T####K`` and writes ``mlpsscha.hdf5``, the
post-transient average as ``force_constants.hdf5`` and ``fc2.hdf5``,
``phonopy_sscha.yaml``, and ``summary.yaml``. Re-running unchanged settings
resumes completed iterations.

With ``run-transport: true``, the ``kappa`` section supplies the phono3py mesh
and method. phono3py is run separately with each FC2(T) and the fixed FC3. This
fixed-FC3 treatment is an approximation and should be stated when reporting
temperature-renormalized transport.

.. important::

   The exported FC2 is an **auxiliary harmonic** matrix. This workflow does not
   calculate the free-energy Hessian. Bubble postprocessing is off by default;
   the optional approximation below does not change the exported FC2.
   Three-phonon scattering in the subsequent transport step does not automatically
   correct the phonon frequencies by the real part of that self-energy.
   The historical command name ``scph`` remains for input compatibility; it does
   not select an FC4-based SCPH-plus-bubble solver.

Optional on-shell bubble frequency correction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Add these keys to the existing ``scph`` section (do not create a second section):

.. code-block:: yaml

   scph:
     bubble: true
     bubble-mesh: [9, 9, 9]
     bubble-epsilons: [0.05, 0.1]  # THz, principal-value regularization
     # bubble-grid-points: [0]    # Optional phono3py BZ grid indices for a pilot

``nepkappa run input.yaml`` / ``nepkappa stage scph input.yaml`` applies the correction
after each SSCHA temperature when ``bubble`` is true. To reuse already completed
SSCHA results without evaluating an ASE calculator or running transport:

.. code-block:: bash

   nepkappa validate input.yaml --for bubble
   nepkappa stage bubble input.yaml
   nepkappa report input.yaml

The explicit ``bubble`` command runs postprocessing regardless of the automatic
``scph.bubble`` switch. It uses ``scph.temps`` and ``scph.workdir`` to locate
``T####K/fc2.hdf5``, ``phonopy_sscha.yaml``, and a completed ``summary.yaml``.
The usual ``scph.transport-metadata`` and ``scph.transport-fc3`` overrides select
matching phono3py metadata and FC3, otherwise those files are taken from
``output.result-dir``. For a ``qha-sscha`` preset/section, the command instead
selects the matching FC3 and metadata inside each ``qha-sscha/T####K`` case;
fixed-volume overrides are not used. No sampling or force-constant generation
is triggered by this postprocessing command.

The calculation evaluates the real part of the cubic bubble on the auxiliary
frequencies, using the supplied FC3. All modes are retained together so that
phono3py's degenerate-mode averaging is preserved. This is an **input-FC3,
diagonal, one-shot on-shell approximation**, not an SSCHA ensemble-averaged
cubic vertex or a free-energy Hessian. It does not include off-diagonal mode
mixing, solve a frequency-dependent Dyson equation, or calculate a spectral
function. In ordinary-frequency THz units the saved estimates are
``frequency_linear = nu + Delta`` and
``frequency_on_shell = sqrt(nu**2 + 2*nu*Delta)``. No extra ``2*pi`` is applied
to phono3py's THz output. The square-root estimate is not a self-consistent root.

Outputs are isolated under ``T####K/bubble/<input-hash>/``:

- ``inputs.yaml``: input hashes, settings and backend versions;
- ``gp-*.hdf5``: restartable per-q checkpoints;
- ``bubble.hdf5``: auxiliary frequencies, Delta, both frequency estimates,
  q coordinates, epsilons and validity masks;
- ``bubble-summary.yaml``: method limits and diagnostic warnings, included in
  ``nepkappa report``.

Delta and corrected arrays have axes ``(epsilon, grid_point, band)``; auxiliary
frequencies have axes ``(grid_point, band)``. The squared-frequency array has
units THz squared. Modes below ``scph.cutoff-frequency`` are excluded from
frequency estimates, and nonpositive corrected squared frequencies are saved
as invalid/NaN, not clipped. Meshes with imaginary auxiliary frequencies below
-1e-4 THz are rejected. NAC parameters are taken only from the saved SSCHA YAML;
there is no directional LO limit at Gamma in this first implementation.

The default mesh and epsilon values are pilot settings, not convergence results.
Converge both together. Explicit grid-point selection is not a full-BZ integral;
otherwise all irreducible mesh points and their weights are saved. A conservative
2 GiB bound on the dense interaction array rejects oversized jobs before FC3
loading; this does not bound total RSS. Large-cell calculations need a tiled
backend instead of silently reducing the physical model.

**No bubble linewidth is added to an existing three-phonon scattering rate.**
Neither FC2 nor existing thermal conductivity is modified. Consistent
bubble-corrected transport requires further development, not simply feeding
these shifted numbers into an unchanged group velocity and scattering table.
See the `phono3py self-energy conventions
<https://phonopy.github.io/phono3py/command-options.html#imaginary-and-real-parts-of-self-energy>`_
and the `SSCHA spectral tutorial <https://sscha.eu/Tutorials/tutorial_spectral/>`_.

For SSCHA at QHA equilibrium volumes, use explicit transport switches:

.. code-block:: yaml

   qha-sscha:
     three-phonon: true
     four-phonon: false

``three-phonon`` runs phono3py with FC2(T,V(T)) and FC3(V(T)).
``four-phonon`` runs FourPhonon with FC2(T,V(T)), FC3(V(T)), and FC4(V(T)); it
therefore requires a ``fourphonon`` section and a working Fourthorder/FourPhonon
installation. The legacy ``scph.run-transport`` option is still accepted as an
alias for ``qha-sscha.three-phonon`` when the new switch is omitted.
Here FC2 is explicitly SSCHA-renormalized, whereas FC3 and FC4 are regenerated
at the QHA equilibrium volume but are not themselves temperature-renormalized.

The full ``scph`` preset expands to ``fc2fc3 -> scph``. Use
``nepkappa stage scph input.yaml`` with explicit existing file paths to reuse a
completed harmonic calculation. ``snapshots``, ``iterations``, ``transient``,
supercell size, and the fitting mesh all require convergence testing.

``qha-sscha`` first runs QHA, interpolates the equilibrium primitive-cell
volume at every ``scph.temps`` point, isotropically rescales the input cell,
optionally relaxes internal coordinates, and regenerates matching FC2. When
``qha-sscha.three-phonon: true`` it also regenerates FC3 at every QHA volume
before running SSCHA and phono3py. With ``qha-sscha.four-phonon: true`` it also
regenerates FC4, exports the temperature-renormalized FC2 in ShengBTE format,
and runs FourPhonon. Outputs are stored under
``qha-sscha/T####K/``. Nested force-constant Slurm arrays are not supported;
neither are nested FourPhonon submissions. Submit the entire ``nepkappa run``
as one batch job instead.

Interpretation of the volume coupling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The legacy name ``QHA+SSCHA`` denotes **SSCHA at QHA equilibrium volumes**.
No QHA and SSCHA frequencies, force constants, or conductivities are added.
The potential is sampled again at each rescaled structure. Nevertheless, the
volume is prescribed by QHA: the SSCHA free energy does not feed back into the
equilibrium volume, lattice shape, or finite-temperature centroid optimization.
Optional ``qha.relax-internal`` is a static fixed-cell relaxation, not a thermal
SSCHA centroid optimization.

This is a sequential approximation, not fully self-consistent anharmonic thermal
expansion. Validate the QHA volumes against SSCHA free-energy/stress calculations
or suitable experimental data before claiming quantitative accuracy. Adding a
bubble spectral correction would not, by itself, correct these QHA volumes.
Once a cell is optimized using SSCHA free energy, an additional QHA expansion
correction must not be applied to it.

New summaries record these limitations in ``approximation``. Reports preserve
that metadata; older results without it are marked as unrecorded rather than
being silently reclassified as a more complete calculation.

``kappa``
-----------

.. code-block:: yaml

   kappa:
     mesh: [21, 21, 21]
     temps: [100, 1000, 50]
     method: rta
     isotope: false
     bfmp: 1.0e6
     wigner: false

``temps`` accepts either one temperature or ``[tmin, tmax, tstep]``.
``method`` can be ``rta`` or ``lbte``. ``isotope: true`` adds phono3py's
``--isotope`` option. ``bfmp`` maps to phono3py's ``--boundary-mfp`` option;
its unit is micrometer, and the phono3py CLI default is ``1.0e6``. ``wigner:
true`` enables phono3py's experimental SMM19 Wigner implementation through
``--tt smm19``.

Custom phono3py command
~~~~~~~~~~~~~~~~~~~~~~~

Advanced users can bypass NEP-kappa's automatic kappa command builder and write
the phono3py command directly:

.. code-block:: yaml

   kappa:
     command: phono3py phono3py_disp.yaml --fc2 --fc3 --br --nu --mesh 21 21 21 --tmin 100 --tmax 1000 --tstep 50

The command runs inside ``output.result_dir``. Relative paths such as
``phono3py_disp.yaml``, ``fc2.hdf5``, and ``fc3.hdf5`` therefore refer to files
in the result directory. When ``command`` is set, NEP-kappa skips the automatic
``mesh``/``method``/scattering option builder and runs this command instead.
The command is split like a normal command line; shell features such as pipes
and redirects are not interpreted.

Distributed LBTE with Slurm
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The normal ``nepkappa stage kappa`` command can submit an LBTE calculation split over
irreducible grid points:

.. code-block:: yaml

   kappa:
     mesh: [21, 21, 21]
     temps: [300]
     method: lbte
     parallel:
       backend: slurm
       jobs: 32

NEP-kappa first runs phono3py with ``--wgp`` to obtain irreducible grid points.
It then submits a ``--write-phonon`` preparation job, a ``--write-pp`` Slurm
array, and a final ``--read-pp`` collection job. Slurm ``afterok`` dependencies
ensure that each stage starts only after the previous stage succeeds.

``jobs`` controls the number of grid-point chunks and defaults to 32.
``max-concurrent`` limits
the number of array tasks running simultaneously. The ``collect-*`` options can
reserve more memory and CPUs for construction and diagonalization of the full
LBTE collision matrix. ``partition``, ``account``, ``extra-sbatch``, and
``preamble`` are optional cluster-specific settings. Each array element is one
single-node, single-task phono3py process; use ``jobs`` rather than ``ntasks``
to distribute grid points. When resource overrides are omitted, Slurm's default
time, memory, CPU, partition, and account settings are preserved.

Generated scripts, grid-point lists, logs, and ``submission.yaml`` are written
under ``output.result_dir/lbte-slurm``. Set ``submit: false`` to generate these
files without invoking ``sbatch``. The result directory has to be on a shared
filesystem visible to all Slurm nodes.

``plot``
----------

.. code-block:: yaml

   plot:
     layout: separate
     path: seekpath
     tau: total
     kappa: all
     temperature: 300
     dpi: 300

``layout`` controls how figures are written:

- ``separate``: write one PNG per requested figure
- ``combined``: write one automatically arranged multi-panel ``combined.png`` figure
- ``both``: write separate PNG files and ``combined.png``

NEP-kappa chooses plots from the available data. With ``fc2.hdf5`` and matching
structure metadata (``phono3py_disp.yaml``, ``phonopy.yaml``, or
``phonopy_disp.yaml``), it generates ``dispersion``, ``dos``, ``heat_capacity``
and ``group_velocity`` without requiring FC3 or a kappa HDF5 file. A bare FC2
array is insufficient: matching cell, primitive and supercell information is
required. Four harmonic figures use a 2-by-2 combined layout.

For harmonic-only plotting, ``kappa.mesh`` controls DOS, group-velocity and
thermal-property sampling, and ``kappa.temps`` specifies the heat-capacity
temperatures as ``[T]`` or ``[start, stop, step]``. These settings do not launch
a conductivity calculation. Heat capacity is converted from Phonopy's molar
units to J m^-3 K^-1 using the primitive-cell volume (and the configured film
or wire effective geometry). Imaginary modes are retained in the dispersion;
if present, a warning explains that the heat capacity excludes nonpositive
modes and must not be interpreted as stable-phase thermodynamics.

When transport data are present, ``relaxation_time``, ``scattering_rate`` and
``kappa`` are included. ``cumulative_kappa`` additionally requires
``mode_kappa``. Comparisons use the figures supported by every dataset, so
harmonic-only datasets can also be compared. ``layout: separate``,
``combined`` and ``both`` apply to both harmonic-only and transport plots.
For harmonic-only comparisons with no configured mesh, the fallback is
21-by-21-by-21 and heat-capacity temperatures are 100--1000 K in 100 K steps;
these are plotting defaults, not convergence claims.

Transport-file selection is explicit: a configured mesh must match an existing
``kappa-m*.hdf5`` filename when conductivity files are present. Without a mesh,
multiple candidate files raise an ambiguity error instead of selecting one
silently. If linewidths are absent, scattering-rate and lifetime plots are
omitted while other supported plots remain available. A requested plot
temperature not in the file uses the nearest stored temperature with a warning.
Unsupported tensor dimensions, split-grid files and incomplete transport
datasets produce explicit errors rather than silently mixing incompatible data.
For phono3py Wigner results, ``kappa_intra`` and ``kappa_inter`` are shown as
particle and coherence contributions alongside the stored total ``kappa``;
the reader checks that the components add to the total. Wigner component plots
require a complete standard kappa HDF5 file, not a component-only file.

When ``fourphonon/fourphonon-summary.yaml`` is present, ``nepkappa plot`` also
reads completed FourPhonon conductivity tables. It includes each available
RTA or iterative solution in ``kappa.png`` and, for Wigner_Park, particle,
coherence and total. A matching phono3py kappa HDF5 adds a separate 3ph-only
reference curve. Missing reference calculations are never inferred. The
``BTE.w_3ph`` and ``BTE.w_4ph`` files in the nearest temperature directory
produce two separate scattering-rate figures. Builds that additionally save
``BTE.w_3ph4ph_NU`` or ``BTE.w_3ph_NU`` produce an N/U overlay; ordinary
FourPhonon outputs do not supply that decomposition.

``path`` controls the high-symmetry path used for phonon dispersion:

- ``seekpath``: determine the path automatically with seekpath
- ``custom``: use the user-defined points and segments below

.. code-block:: yaml

   plot:
     path: custom
     path_points:
       G: [0.0, 0.0, 0.0]
       X: [0.5, 0.0, 0.5]
       U: [0.625, 0.25, 0.625]
       K: [0.375, 0.375, 0.75]
       L: [0.5, 0.5, 0.5]
       W: [0.5, 0.25, 0.75]
     path_segments:
       - [G, X]
       - [X, U]
       - [K, G]
       - [G, L]
       - [L, W]
       - [W, X]

When two adjacent segments are disconnected, NEP-kappa combines the labels at
the break point, e.g. ``[X, U]`` followed by ``[K, G]`` is shown as ``U|K``.

``tau`` controls the relaxation-time channel:

- ``total``: total scattering rate from ``gamma``
- ``normal``: normal-process scattering rate from ``gamma_N``
- ``umklapp``: Umklapp-process scattering rate from ``gamma_U``
- ``nu``: plot N and U channels together, without the total channel; both
  ``gamma_N`` and ``gamma_U`` must be present
- ``all``: plot total, N, and U channels together when available

``kappa`` controls the thermal-conductivity components shown in the kappa
figure:

- ``x``: plot only ``kappa_xx``
- ``y``: plot only ``kappa_yy``
- ``z``: plot only ``kappa_zz``
- ``all``: plot ``kappa_xx``, ``kappa_yy``, ``kappa_zz``, and their average for
  a single solution; when Wigner components or multiple transport solutions
  are present, compare their spatial averages on one axis

``temperature`` selects the target temperature for the relaxation-time and
scattering-rate plots. The scattering rate is ``4 pi gamma`` in ps^-1, the
inverse of the plotted lifetime convention ``tau = 1/(4 pi gamma)``.
The same temperature is used for ``cumulative_kappa.png``; its per-mode values
are divided by the full q-point mesh size before frequency accumulation, as in
phono3py's total-kappa reduction.
NEP-kappa uses the closest temperature available in ``kappa-m*.hdf5``.

``output``
------------

.. code-block:: yaml

   output:
     progress: true
     result_dir: calculations/example-runs/nep-rta

All generated files are written inside ``result_dir``. The terminal output is
also saved to ``run.log`` in that directory.

``plot`` output
-----------------

``nepkappa plot`` reads matching phonon metadata and ``fc2.hdf5`` from
``output.result_dir``. A compatible ``kappa-m*.hdf5`` adds transport panels;
it is not required for harmonic-only figures. The command writes figures under
``output.result_dir/plots``. The figures follow a publication-oriented style
with larger axis labels, tick labels, line widths, and marker sizes. Subplot
titles are intentionally omitted so the figures are easier to compose in papers.

- ``dispersion.png``: phonon dispersion along seekpath high-symmetry lines
- ``dos.png``: phonon density of states
- ``heat_capacity.png``: volume heat capacity
- ``group_velocity.png``: group velocity magnitude in km/s
- ``relaxation_time.png``: relaxation time at the nearest available ``plot.temperature`` (default 300 K), when linewidths exist
- ``scattering_rate.png``: available total, Normal, and/or Umklapp scattering rates at that temperature
- ``scattering_rate_3ph.png`` and ``scattering_rate_4ph.png``: separate
  FourPhonon rates, when their ``BTE.w_*`` files exist
- ``scattering_rate_nu.png``: optional FourPhonon N and U overlay when an
  ``*_NU`` file was written by the solver
- ``cumulative_kappa.png``: cumulative conductivity against phonon frequency when ``mode_kappa`` is available
- ``kappa.png``: selected conductivity component, or an average comparison
  of available Wigner contributions and/or 3ph+4ph solver schemes
- ``combined.png``: automatically arranged multi-panel figure when ``layout`` is ``combined`` or ``both``
