Input Files
===========

NEP-kappa uses six sections for calculations starting from a structure.
Phonon-property and thermal-transport calculations use:
``structure``, ``calculator``, ``force-constant``, ``kappa``, ``plot``, and
``output``. TD-BTE replaces ``kappa`` with ``tdbte``. All six section names
must appear; a section without options may be written as ``{}``. Optional
force-job execution resources may be placed in a seventh ``parallel`` section.
For existing TD-BTE force constants or kernels, a compact two-section
``tdbte``/``output`` reuse input is also supported. Do not combine ``kappa``
and ``tdbte`` in one input.

A compact 3C-SiC phonon/transport input is:

.. code-block:: yaml

   structure:
     poscar: examples/structures/3C-SiC/POSCAR_primitive
     relaxation: true
   calculator:
     name: nep
     nep-model: potentials/3C-SiC/nep_3C-SiC.txt
   force-constant:
     dim-fc2: [3, 3, 3]
     dim-fc3: [3, 3, 3]
   kappa:
     engine: phono3py
     mesh: [31, 31, 31]
     temps: [100, 1000, 50]
     method: rta
   plot:
     layout: both
   output:
     result-dir: calculations/example-runs/3c-sic-nep-rta

These numerical settings demonstrate the workflow; they are not converged
parameters. ``nepkappa info input.yaml`` parses and checks the configuration
without starting a calculation. It does not validate the physical quality of
a potential, external executables, or convergence.

Command behavior
----------------

.. code-block:: text

   nepkappa info input.yaml     check and show the input
   nepkappa run input.yaml      run the input-selected calculation
   nepkappa relax input.yaml    relax the structure
   nepkappa fc2 input.yaml      generate FC2, or scan QHA volumes when enabled
   nepkappa fc2fc3 input.yaml   generate FC2 and FC3
   nepkappa kappa input.yaml    run the selected transport route
   nepkappa tdbte input.yaml    run nonequilibrium phonon dynamics
   nepkappa plot input.yaml     plot completed phonon/transport results
   nepkappa report input.yaml   write a result summary

The six sections describe a complete phonon/transport calculation, but a direct command
runs only its named part. If ``structure.relaxation`` is enabled, run
``nepkappa relax`` before a separate ``fc2`` or ``fc2fc3`` command; these
commands then read ``POSCAR_relaxed``. ``nepkappa run`` carries out the whole
input-selected sequence. QHA, SSCHA, and four-phonon settings are switches
within that one run; they select additional stages and their supported
dependencies. ``plot`` is run separately after the necessary data exist.
Internal FC4, SSCHA, bubble, and coupled-transport steps are selected by input
switches; they have no separate public command.
For a QHA-only calculation, set ``force-constant.qha.enabled: true`` and
``force-constant.qha.run-transport: false``, then execute
``nepkappa fc2 input.yaml``. This scans the configured volumes and generates
their FC2 and static energies; it does not calculate conductivity. For
QHA+SSCHA, enable both sections and use ``nepkappa run input.yaml`` so QHA
finishes before SSCHA reads its temperature-dependent equilibrium volumes.

Phonon and thermal-transport sections
-------------------------------------

``structure`` sets ``poscar`` and ``relaxation``. For a VASP relaxation with
custom stages, ``relaxation`` may be a mapping:

.. code-block:: yaml

   structure:
     poscar: examples/structures/3C-SiC/POSCAR_primitive
     relaxation:
       enabled: true
       workdir: vasp-relax
       stages:
         coarse: {ediffg: -0.05, prec: Normal}
         fine: {ediffg: -0.01, prec: Accurate}

``dimensionality`` is 1, 2, or 3. For films set ``effective-thickness`` in
Angstrom and optionally ``vacuum-axis``; for wires set ``effective-area`` in
Angstrom squared and optionally ``periodic-axis``. These geometry factors
correct plotted heat capacity and conductivity. They do not modify the original
phono3py HDF5 file. QHA reads the specified POSCAR directly and therefore
requires ``relaxation: false`` in this input.

The present QHA, SSCHA, and FourPhonon implementations accept only
``dimensionality: 3``. This is a capability limit of those stages, independent
of whether ``run`` orchestrates them: the QHA scan scales the full cell volume,
and the current SSCHA/FourPhonon routes have not been validated with effective
film thickness or wire area. Setting a low-dimensional geometry factor alone
does not establish a physically consistent low-dimensional calculation.

``calculator`` selects the force and energy backend. NEP needs ``nep-model``.
VASP uses ``name: vasp``, a ``vasp-command`` or ``vasp-path``, and
``potcar-path``. The executable and licensed POTCAR files must exist on the
execution host. Native MACE accepts ``model``, ``foundation``, ``device``, and
``dtype``. Another ASE-compatible calculator can use ``name: ase`` with
``factory``, ``kwargs``, and optional ``model-files``. An installed calculator
plugin may supply a short name. An input check does not launch these programs.

``force-constant`` sets ``dim-fc2`` and ``dim-fc3`` independently. By default,
phono3py finite displacements create compact HDF5 FC2/FC3; ``use-hiphive``
selects fitting instead. ``fc3-backend: thirdorder`` selects Thirdorder for
FC3. ``cutoff-fc3`` is a positive displaced-pair cutoff in Angstrom for native
phono3py, or a negative integer neighbor shell / positive nm cutoff for
Thirdorder. HiPhive accepts ``n-structures``, ``rattle-std``, ``cutoffs``, and
``min-dist``. ``format`` may be ``phono3py``, ``shengbte``, or ``both``.
Phono3py transport needs HDF5 FC2/FC3; FourPhonon needs ShengBTE-format FC2/FC3
and FC4. The FC4 supercell and cutoff are ``dim-fc4`` and ``cutoff-fc4``.
VASP force calculations can use ``workdir`` and ``vasp-kwargs`` inside this
section.

QHA, SSCHA, and four-phonon work are switches inside ``force-constant``. Each
is ``false``, ``true``, or a mapping with ``enabled`` and options. For example:

.. code-block:: yaml

   force-constant:
     dim-fc2: [3, 3, 3]
     dim-fc3: [3, 3, 3]
     qha:
       enabled: true
       volume-ratios: [0.94, 0.97, 1.00, 1.03, 1.06]
       temps: [0, 1000, 10]
     sscha:
       enabled: true
       run-transport: true
       temps: [300, 900, 300]
       snapshots: 200
       iterations: 4
     four-phonon: false
   kappa:
     engine: phono3py
     temps: [300, 900, 300]
     method: rta

QHA scans isotropically scaled volumes, fits an equation of state, and reports
thermal expansion. Its options include ``dim-fc2``, ``mesh``, ``eos``,
``pressure``, ``relax-internal``, and optional VASP relaxation/static kwargs.
For QHA-derived conductivity, set ``force-constant.qha.run-transport: true``:
NEP-kappa
regenerates FC2 and FC3 at each requested temperature's QHA equilibrium volume
and runs phono3py RTA there. With QHA enabled and SSCHA disabled, this route
supports phono3py RTA only. Leave ``run-transport: false`` to finish after QHA.
QHA alone does not add explicit anharmonic frequency renormalization.

SSCHA uses the configured NEP/MACE/ASE calculator to refit an auxiliary
harmonic FC2 at each requested temperature. Common options are ``snapshots``,
``iterations``, ``transient``, ``sscha-mesh``, and ``initial-fc2``.
When QHA and SSCHA are combined, set
``force-constant.sscha.run-transport: true`` to use both the QHA volume and
the SSCHA FC2 for thermal transport. An explicit opposing QHA switch is rejected.
Enabling ``bubble`` saves separate diagonal on-shell
frequency shifts from the input FC3; those shifts do not update the exported
FC2 or conductivity. Standalone SSCHA has no public example; the maintained
QHA+SSCHA example runs SSCHA at each QHA equilibrium volume. This sequential
route does not minimize volume using SSCHA free energy. FC3/FC4 are regenerated
at that volume, but are not themselves temperature-renormalized.

With ``four-phonon`` enabled, ``nepkappa run`` generates FC4 and computes
three-plus-four-phonon transport through FourPhonon. Its mapping accepts
``command``, ``control``, ``harmonic-format``, ``sample-*`` settings, and solver
resources. Automatically generated CONTROL files disable non-analytic
corrections; polar crystals such as SiC need validated dielectric and Born-charge
inputs when those corrections matter. A FourPhonon executable and Fourthorder
are external dependencies. The ``Wigner_Park`` executable is required for the
four-phonon Wigner mode. Its coherence term uses combined three-phonon,
isotope, and four-phonon linewidths; ``nepkappa report`` records population,
coherence, and total conductivity separately. See
``examples/nep-rta-wigner-3ph-4ph.yaml`` for a complete configuration.

``kappa`` sets the transport engine, q mesh, temperatures, and solver. It does
not choose whether QHA or SSCHA corrections are used: their ``run-transport``
switches do that. Ordinary three-phonon and enabled four-phonon inputs run
transport by default. An enabled QHA or SSCHA section defaults to
``run-transport: false``. When QHA and SSCHA are both enabled, either true
switch selects their combined transport route; explicitly conflicting values
are rejected.
``temps`` accepts ``[T]`` or
``[minimum, maximum, step]``. For phono3py, use
``engine: phono3py`` and ``method: rta`` or ``lbte``. Optional ``isotope``,
``bfmp`` (micrometers), and ``wigner`` select scattering and SMM19 coherence.
For FourPhonon, use ``engine: fourphonon`` and set ``method-3ph`` and
``method-4ph`` separately. Supported pairs are ``rta/rta``, ``lbte/rta``, and
``lbte/lbte``; ``rta/lbte`` is unavailable. Four-phonon Wigner currently needs
``rta/rta`` and an appropriate Wigner_Park executable.

Temperature settings have distinct purposes and are validated together:

* ``force-constant.qha.temps`` is the QHA volume-scan grid, written as
  ``[Tmin, Tmax, step]`` with ``0 <= Tmin < Tmax``.
* ``force-constant.sscha.temps`` is the SSCHA sampling grid, written as
  ``[Tmin, Tmax, step]`` with ``0 <= Tmin <= Tmax``. A single SSCHA point can
  be written ``[T, T, step]`` with a positive step.
* ``kappa.temps`` is the set of conductivity temperatures: ``[T]`` or
  ``[Tmin, Tmax, step]``, with temperatures strictly above 0 K. An arbitrary
  list of unevenly spaced temperatures is not currently accepted.

Every step must be positive and divide its interval exactly, including the
end point. With QHA and SSCHA enabled, the SSCHA range must fit inside the
QHA range; SSCHA points may lie between QHA grid points because equilibrium
volumes are interpolated. When transport is selected, all conductivity
temperatures must lie inside the QHA range (if QHA is enabled). If SSCHA is
also enabled, every conductivity temperature must be one of the SSCHA sampling
points. Thus ``kappa.temps`` may be a subset of ``sscha.temps``: SSCHA runs at
all its points, and transport runs only at selected points. For example, QHA
``[0, 1000, 10]``, SSCHA ``[300, 900, 100]``, and conductivity
``[300, 900, 300]`` are compatible. For QHA or SSCHA transport,
``kappa.temps`` must be written explicitly. If transport is disabled, ``kappa.temps``
does not constrain the QHA or SSCHA temperature range. FourPhonon transport
uses ``kappa.temps`` too; any legacy ``four-phonon.temps`` value must match it.

``plot`` configures later plotting. ``layout`` is ``separate``, ``combined``,
or ``both``; ``path`` selects the dispersion path; ``tau`` selects the
scattering-rate channel; ``temperature`` and ``kappa`` choose the displayed
transport data; ``dpi`` controls PNG resolution. With FC2 and matching cell
metadata, ``nepkappa plot`` makes dispersion, DOS, heat-capacity, and
group-velocity plots. With conductivity data it adds conductivity; linewidth
data also enables scattering rate, and mode conductivity enables cumulative
conductivity. RTA output can show
normal and Umklapp rates; Wigner output shows particle, coherence, and total
conductivity. FourPhonon output adds separate 3ph/4ph rates and compares the
available conductivity solutions. A bare FC2 array without cell metadata is
insufficient for dispersion plotting.

``output`` sets ``result-dir`` and optional ``progress``. Use a distinct result
directory for each material and parameter set so that outputs are not mixed.

Parallel execution
------------------

An optional ``parallel`` section keeps site-specific execution resources out
of the six scientific sections:

.. code-block:: yaml

   parallel:
     force-constant:
       backend: slurm
       jobs: 16
       submit: false
     kappa:
       backend: slurm
       submit: false

``parallel.force-constant`` controls displacement-force jobs.
``parallel.kappa`` controls phono3py LBTE or FourPhonon jobs according to the
engine. ``submit: false`` writes scripts for review without calling ``sbatch``.
MPI/OpenMP settings for FourPhonon may also be specified there. QHA-coupled
calculations should be submitted as a whole allocation; nested force-job arrays
are not supported. Force-job caches and Slurm collection state are kept under
``output.result-dir``.

Time-dependent BTE input
------------------------

For a new calculation, TD-BTE uses six sections. ``nepkappa run`` generates
matching FC2, FC3, and displacement metadata before building the collision
kernel and evolving the excitation. It does not calculate conductivity:

.. code-block:: yaml

   structure:
     poscar: examples/structures/3C-SiC/POSCAR_primitive
     relaxation: false
   calculator:
     name: nep
     nep-model: potentials/3C-SiC/nep_3C-SiC.txt
   force-constant:
     dim-fc2: [3, 3, 3]
     dim-fc3: [3, 3, 3]
   tdbte:
     mesh: [5, 5, 5]
     temperature: 300
     branches: [3, 4, 5]
     duration-ps: 10
   plot:
     layout: both
     dpi: 300
   output:
     result-dir: calculations/3C-SiC/tdbte-from-scratch

Only phono3py FC2/FC3 are supported in this route; QHA, SSCHA, FC4, and
FourPhonon settings do not feed the present TD-BTE solver. ``plot`` accepts
``layout`` and ``dpi`` for saved TD-BTE energy plots; run ``nepkappa plot``
after the calculation to generate them. ``nepkappa tdbte`` with the same
six-section file uses the already generated files in ``output.result-dir``
without recalculating forces.

For existing force constants, the compact two-section input remains available:

.. code-block:: yaml

   tdbte:
     force-constants: calculations/3C-SiC/fc
     mesh: [5, 5, 5]
     branches: [3, 4, 5]
   output:
     result-dir: calculations/3C-SiC/tdbte-reuse

The ``force-constants`` directory must contain mutually matching FC2, FC3,
and phono3py displacement metadata. A previously built ``kernel`` may be
used instead. A prior conductivity result is never a prerequisite.
``excitation``, ``max-step-ps``, and ``samples`` control the initial perturbation
and time integration. See :doc:`tdbte` for the numerical model and audits.
