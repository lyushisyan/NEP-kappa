Input Files
=============

NEP-kappa 2.0.0 uses YAML input files. Each YAML file describes one workflow by
grouping settings into the same stages used by the command line:

- ``workflow``: high-level preset or an advanced custom stage plan
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
- ``output``: progress display and result directory

The ``compare`` command uses a smaller YAML file with ``datasets``, ``compare``,
and ``plot`` sections. The original two-way ``reference``/``candidate`` schema
remains supported for compatibility.

YAML validation is strict. Unknown sections and keys are rejected before a
calculation starts, and likely misspellings include a suggested supported key.
Free-form VASP ``vasp_kwargs`` and relaxation-stage mappings remain available
for native VASP options.

Command behavior
------------------

Normal users select ``workflow.preset`` and use only:

.. code-block:: bash

   nepkappa run input.yaml

Presets are ``three-phonon``, ``four-phonon``, ``qha``, ``scph``, and
``qha-sscha``. The
following commands expose individual stages for advanced use:

.. code-block:: bash

   nepkappa relax input.yaml
   nepkappa fc2 input.yaml
   nepkappa fc2fc3 input.yaml
   nepkappa fc4 input.yaml
   nepkappa qha input.yaml
   nepkappa scph input.yaml
   nepkappa qha-sscha input.yaml
   nepkappa kappa input.yaml
   nepkappa kappa4 input.yaml
   nepkappa plot input.yaml
   nepkappa compare compare.yaml
   nepkappa info input.yaml
   nepkappa report calculations/runs/calculation

- ``nepkappa relax`` relaxes the structure and writes ``POSCAR_relaxed`` to ``output.result_dir``.
- ``nepkappa fc2`` generates ``phono3py_disp.yaml`` and ``fc2.hdf5`` only.
- ``nepkappa fc2fc3`` generates ``phono3py_disp.yaml``, ``fc2.hdf5``, and ``fc3.hdf5``.
- ``nepkappa fc4`` generates FourPhonon ``FORCE_CONSTANTS_4TH`` using ``Fourthorder_vasp.py``.
- ``nepkappa qha`` computes isotropic quasi-harmonic thermal properties over a volume scan.
- ``nepkappa scph`` runs Phonopy stochastic SSCHA.
- ``nepkappa qha-sscha`` runs SSCHA at volumes interpolated from a completed QHA fit.
- ``nepkappa kappa`` computes thermal conductivity using existing ``phono3py_disp.yaml``, ``fc2.hdf5``, and ``fc3.hdf5``.
- ``nepkappa kappa4`` runs FourPhonon with existing ShengBTE-format force constants.
- ``nepkappa plot`` creates standard plots from ``fc2.hdf5`` and ``kappa-m*.hdf5``.
- ``nepkappa compare`` overlays DFT and multiple potential-model result directories in the same standard figures.
- ``nepkappa converge`` generates and analyzes a parameter sweep from one base workflow YAML.
- ``nepkappa report`` writes ``report.yaml`` and ``report.md`` from an existing result tree.
- ``nepkappa run`` expands and executes the selected workflow preset. Without a
  ``workflow`` section it preserves the legacy ``relax`` + ``fc2fc3`` +
  ``kappa`` behavior.
- ``nepkappa info`` prints the parsed configuration without running a calculation.

``nepkappa fc2fc3`` computes and writes FC2 first, then starts the FC3
displacement, force, and export stage.

When ``relaxation.enabled`` is ``true``, ``nepkappa fc2`` and
``nepkappa fc2fc3`` read
``POSCAR_relaxed`` from ``output.result_dir``. Run ``nepkappa relax`` first, or
use ``nepkappa run``.

Configuration validation is command-aware. For example, ``nepkappa scph``
validates its calculator, structure, and ``scph`` sections, while
``nepkappa kappa`` validates transport settings without requiring QHA or
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

``nepkappa kappa4`` stages ShengBTE-format ``FORCE_CONSTANTS_2ND``,
``FORCE_CONSTANTS_3RD``, and ``FORCE_CONSTANTS_4TH`` and runs FourPhonon.
Set ``harmonic-format: espresso`` when the harmonic input is an official
``espresso.ifc2`` file; this mode requires a matching custom ``CONTROL``.
Solver choices are ``rta``, ``3ph-iterative``, and ``full-iterative``.
Sampling settings expose FourPhonon's 3ph/4ph scattering and phase-space
process counts. MPI+OpenMP and Slurm resources are configured in the same
section. Results are normalized to ``kappa4-rta.dat`` and
``kappa4-iterative.dat`` and summarized in ``fourphonon-summary.yaml``.

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

The current ``nepkappa kappa`` command reads phono3py HDF5 files. Use
``format: phono3py`` or ``format: both`` when the next stage is
``nepkappa kappa``. Use ``format: shengbte`` when the next stage is a
ShengBTE/FourPhonon workflow.

For native phono3py FC3, ``pair-cutoff-fc3`` is an optional displaced-pair
distance cutoff in Angstrom. For ``fc3-backend: thirdorder``, use
``cutoff-fc3`` instead: a negative integer selects a neighbor shell and a
positive value is a distance in nm. HiPhive uses its own real-space
``cutoffs`` list in Angstrom. Native finite-displacement FC2 has no independent
real-space cutoff: ``dim-fc2`` controls the represented interaction range.
Blindly zeroing fitted FC2 elements after the calculation is not provided
because it can break translational/rotational sum rules.

Thirdorder FC3 route:

.. code-block:: yaml

   force-constant:
     dim-fc3: [4, 4, 4]
     fc3-backend: thirdorder
     cutoff-fc3: -3
     thirdorder-command: thirdorder_vasp.py
     fc3-workdir: fc3-thirdorder-runs
     format: shengbte

With ``fc3-backend: thirdorder``, ``nepkappa fc2fc3`` generates FC2 first,
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

``nepkappa fc4`` runs ``Fourthorder_vasp.py sow`` in
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

For the coupled QHA+SSCHA workflow, use explicit transport switches:

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
``nepkappa scph input.yaml`` with explicit existing file paths to reuse a
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

The normal ``nepkappa kappa`` command can submit an LBTE calculation split over
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

NEP-kappa always generates the seven standard plots: ``dispersion``, ``dos``,
``heat_capacity``, ``group_velocity``, ``relaxation_time``,
``scattering_rate``, and ``kappa``. When every input HDF5 contains
``mode_kappa``, ``cumulative_kappa`` is added automatically.

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
- ``all``: plot ``kappa_xx``, ``kappa_yy``, ``kappa_zz``, and their average

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

``nepkappa plot`` reads ``phono3py_disp.yaml``, ``fc2.hdf5``, and
``kappa-m*.hdf5`` from ``output.result_dir`` and writes figures under
``output.result_dir/plots``. The figures follow a publication-oriented style
with larger axis labels, tick labels, line widths, and marker sizes. Subplot
titles are intentionally omitted so the figures are easier to compose in papers.

- ``dispersion.png``: phonon dispersion along seekpath high-symmetry lines
- ``dos.png``: phonon density of states
- ``heat_capacity.png``: volume heat capacity
- ``group_velocity.png``: group velocity magnitude in km/s
- ``relaxation_time.png``: relaxation time near 300 K
- ``scattering_rate.png``: total, Normal, and/or Umklapp scattering rate near 300 K
- ``cumulative_kappa.png``: cumulative conductivity against phonon frequency when ``mode_kappa`` is available
- ``kappa.png``: selected thermal conductivity component or all diagonal components plus average
- ``combined.png``: automatically arranged multi-panel figure when ``layout`` is ``combined`` or ``both``

Multi-model comparison
------------------------

``nepkappa compare compare.yaml`` reads completed DFT and machine-learning
potential result directories, then overlays any number of them in the same
seven standard figures. Each result directory must contain ``phono3py_disp.yaml``,
``fc2.hdf5``, and ``kappa-m*.hdf5``.

.. code-block:: yaml

   datasets:
     - label: DFT
       directory: calculations/example-runs/dft
     - label: NEP-1
       directory: calculations/example-runs/nep-1
     - label: NEP-2
       directory: calculations/example-runs/nep-2
     - label: MACE
       directory: calculations/example-runs/mace

   compare:
     compare_dir: comparison

   plot:
     layout: both
     path: seekpath
     tau: total
     temperature: 300
     kappa: all
     dpi: 300

``compare.compare_dir`` receives ``compare.log`` and a ``plots/`` directory.
The ``layout``, ``path``, ``tau``, ``temperature``, ``kappa``, and ``dpi``
settings follow the same rules as ``nepkappa plot``.

Convergence-study input
-------------------------

``nepkappa converge`` reads a small study YAML rather than a normal workflow
YAML. It copies a base workflow, changes one dotted YAML parameter for every
configured value, and gives every case an isolated result directory.

.. code-block:: yaml

   base: input.yaml
   parameter: kappa.mesh
   values:
     - [15, 15, 15]
     - [21, 21, 21]
     - [27, 27, 27]
     - [31, 31, 31]

   study:
     directory: studies/qmesh
     execute: false
     temperature: 300
     component: average
     tolerance: 0.02

``parameter`` must already exist in the base YAML. Typical paths include
``kappa.mesh``, ``force-constant.dim-fc2``,
``force-constant.dim-fc3``, ``force-constant.dim-fc4``, and
``scph.snapshots``. Unknown paths and duplicate values are rejected.

``execute: false`` prepares inputs and analyzes whichever results already
exist. ``execute: true`` also invokes the base workflow for every case. This
works with local workflows and Slurm-enabled base inputs. The target component
may be ``x``, ``y``, ``z``, or ``average``. The HDF5 temperature closest to
``study.temperature`` is used.

The final configured value is used as the reference. ``convergence.csv`` and
``convergence.png`` show the available results, while
``convergence-summary.yaml`` records relative errors and completion state. To
avoid a premature scientific conclusion, the earliest converged case is
reported only when all configured cases are complete and that case plus every
later case is within ``tolerance`` of the final reference.
