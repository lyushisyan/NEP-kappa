Input files and parameters
==============================

YAML organization
---------------------

``nepkappa run`` requires six sections for its complete calculation workflow.
Phonon-property and thermal-transport calculations use:
``structure``, ``calculator``, ``force-constant``, ``kappa``, ``plot``, and
``output``. TD-BTE replaces ``kappa`` with ``tdbte``. All six section names
must appear; a section without options may be written as ``{}``. Optional
force-job execution resources may be placed in a seventh ``parallel`` section.
Direct commands accept only the sections needed by their stage. For example,
``kappa``/``output`` reuses existing force constants with ``nepkappa kappa``,
and ``tdbte``/``output`` reuses force constants or kernels with ``nepkappa tdbte``.
These compact inputs cannot be executed with ``run``. Do not combine ``kappa``
and ``tdbte`` in one input.

A compact 3C-SiC phonon/transport input is:

.. code-block:: yaml

   structure:
     poscar: examples/structures/3C-SiC/POSCAR_primitive
     relaxation: true
   calculator:
     name: nep
     model: potentials/3C-SiC/nep_3C-SiC.txt
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

Single-stage inputs
-----------------------

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

``run`` requires all six sections. A direct command accepts a reduced input
or a complete input, and runs only its named part:

* ``relax``: ``structure`` and ``calculator``.
* ``fc2`` / ``fc2fc3``: ``structure``, ``calculator``, ``force-constant``.
* Ordinary phono3py ``kappa``: ``kappa``; optionally set
  ``kappa.force-constants`` to the existing FC2/FC3 directory.
* ``tdbte``: ``tdbte``; optionally set ``tdbte.force-constants`` to the
  existing FC2/FC3 directory, or ``tdbte.kernel`` to an existing collision kernel.
* ``plot``: ``plot`` and ``output``; ``output`` alone uses default plot settings.
  Add ``kappa`` for mesh/temperature choices or ``structure`` for geometry
  normalization. With a plot-only input, saved ``tdbte/branch-energy.csv``
  selects TD-BTE plots; ``tdbte: {}`` also explicitly selects this route.
* ``report``: ``output.result-dir``, or pass a result directory directly.

Except for report inputs, ``output`` is optional and defaults to ``result``
below the launch directory. Explicit paths are recommended. Nested QHA,
SSCHA and FourPhonon settings still belong in ``force-constant``; coupled
transport that computes new forces also requires ``structure`` and ``calculator``.
Optional ``parallel`` settings can accompany a reduced input.
Use ``nepkappa info input.yaml --for <command>`` to check a stage input;
``info`` without ``--for`` checks the six-section ``run`` form. ``report`` is
also accepted as an ``info --for`` target.

For direct phono3py ``kappa`` and ``tdbte``, ``force-constants`` in the
corresponding section selects the input directory; ``output.result-dir``
selects where new results are written. Without an explicit source (or a
TD-BTE kernel), both commands read force constants from ``output.result-dir``.
Neither command computes new atomic forces when reusing these inputs.
The source must contain ``phono3py_disp.yaml``, ``fc2.hdf5`` and ``fc3.hdf5``;
missing files cause an error. Relative paths are resolved from the launch
directory, and ``~`` in the source path expands to the user's home directory.

For example, save this input and run ``nepkappa kappa kappa.yaml``:

.. code-block:: yaml

   kappa:
     force-constants: calculations/3C-SiC/fc
     method: rta
     mesh: [5, 5, 5]
     temps: [300]
   output:
     result-dir: calculations/3C-SiC/kappa-reuse

When the directories differ, ``kappa`` copies the three required inputs and
an optional ``BORN`` file into the result directory for phono3py and subsequent
plotting. It leaves the source unchanged and refuses conflicting destination
inputs; identical copies can be reused. This also applies to custom phono3py
commands and Slurm LBTE preparation. TD-BTE reads its source directly.
The explicit ``kappa.force-constants`` option is for ordinary phono3py reuse,
not FourPhonon or coupled QHA/SSCHA transport. Omit explicit force-constant
sources for ``run``, which generates and passes on its own force constants.

If ``structure.relaxation`` is enabled, run
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

structure: crystal structure
--------------------------------

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

calculator: energies and forces
-----------------------------------

``calculator`` selects the force and energy backend. Both NEP and MACE use
``model``: ``name: nep`` selects an NEP potential file, while ``name: mace``
selects a MACE checkpoint or a model identifier accepted by its foundation
factory. The default calculator is NEP, regardless of YAML key order.
The legacy NEP key ``nep-model`` is still accepted, but it cannot appear
together with ``model`` (even when both values are identical).
VASP uses ``name: vasp``, a ``vasp-command`` or ``vasp-path``, and
``potcar-path``. The executable and licensed POTCAR files must exist on the
execution host. Native MACE accepts ``model``, ``foundation``, ``device``, and
``dtype``. Another ASE-compatible calculator can use ``name: ase`` with
``factory``, ``kwargs``, and optional ``model-files``. An installed calculator
plugin may supply a short name. An input check does not launch these programs.

force-constant: force constants and temperature effects
-----------------------------------------------------------

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

For a QHA phonon overview, set
``force-constant.qha.equilibrium-phonons: true`` (default: false). After the
volume fit, ``nepkappa fc2 input.yaml`` also recalculates FC2 at every output
temperature's equilibrium volume, with the same internal-coordinate relaxation
setting. This adds force calculations; it does not interpolate the scanned FC2.
The files are stored under ``qha-phonons/T.../``, with a completed index in
``qha-phonons.yaml``. An equilibrium volume outside the scanned range is rejected.

Then ``nepkappa plot input.yaml`` creates four QHA panels: free energy versus
volume with fitted minima, equilibrium volume versus temperature, volumetric
thermal expansion, and temperature-dependent dispersion with DOS. It only reads
completed data and never starts a force calculation. ``plot.layout`` supports
``separate``, ``combined``, and ``both``; PNG and PDF files are written under
``plots/`` as ``qha_free_energy``, ``qha_volume``, ``qha_thermal_expansion``,
``qha_phonons``, and ``qha_combined``. Dispersion uses the existing ``plot.path``
and custom-path parameters, and DOS uses the mesh saved by the QHA phonon run.
Up to ten evenly spaced available temperature indices, including both endpoints,
are selected for the colored free-energy and phonon curves. The volume and
expansion panels retain all output temperatures; derivative-only auxiliary
temperatures are excluded. If the scan includes volume ratio 1, its phonons are
shown as a dashed static-volume reference. The curves describe volume-induced
frequency changes, not explicit SSCHA renormalization.

Pure-QHA results are recognized by ``qha-summary.yaml``, so a later plot-only
input can contain just ``plot`` and ``output``. QHA transport and SSCHA inputs
retain their existing temperature-specific plot routing. Existing QHA runs
without equilibrium FC2 must be rerun with the switch enabled before generating
this four-panel overview. Already generated QHA PDF files remain usable.

SSCHA uses the configured NEP/MACE/ASE calculator to refit an auxiliary
harmonic FC2 at each requested temperature. Common options are ``snapshots``,
``iterations``, ``transient``, ``sscha-mesh``, and ``initial-fc2``.
When QHA and SSCHA are combined, set
``force-constant.sscha.run-transport: true`` to use both the QHA volume and
the SSCHA FC2 for thermal transport. An explicit opposing QHA switch is rejected.
Enabling ``bubble`` saves separate diagonal on-shell
frequency shifts from the input FC3; those shifts do not update the exported
FC2 or conductivity. See :doc:`harmonic_temperature` for fixed-volume and
QHA+SSCHA inputs. The latter runs SSCHA at each QHA equilibrium volume. This sequential
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

kappa: scattering and thermal conductivity
----------------------------------------------

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

For a controlled three-phonon versus three-plus-four-phonon comparison, set
``force-constant.four-phonon.compare-3ph: true``. The same FourPhonon executable
runs twice in sequence, with identical FC2/FC3, structure, mesh, temperatures,
broadening, isotope and parallel settings. The first run disables four-phonon
scattering and writes under ``fourphonon/threephonon/``; the second uses FC4 and
writes under ``fourphonon/``. With ``method-3ph: lbte`` and ``method-4ph: rta``,
plotting overlays ``3ph RTA``, ``3ph LBTE``, ``3ph+4ph RTA`` and
``3ph LBTE + 4ph RTA``. The RTA/iterative tensors within a single run are not
the three-phonon/four-phonon pair. Any enabled isotope contribution is included
in both cases. ``compare-3ph`` defaults to false and currently requires ordinary
transport with generated CONTROL (not custom CONTROL, Wigner or only-harmonic).
Local runs and Slurm scripts execute both cases sequentially with the requested
resources, rather than doubling the MPI allocation. Use a fresh result directory;
existing comparison tensors are not overwritten. Collection checks that the
staged force constants, CONTROL settings and temperature lists match before
publishing the comparison summary.

``fc3-backend: thirdorder`` generates the ShengBTE FC3 directly from Thirdorder
displacements and the chosen calculator's forces. ``format: both`` preserves
this original ``FORCE_CONSTANTS_3RD`` for FourPhonon and additionally exports HDF5.
It requires a working ``thirdorder_vasp.py`` command (or ``thirdorder-command``);
the filename does not require VASP when NEP supplies the forces. The negative
``cutoff-fc3`` convention selects neighbor shells and must be converged separately
from the FC4 cutoff and transport mesh.

FC2, Thirdorder FC3, Fourthorder FC4 and generated CONTROL use the same
canonical unit-cell coordinates. Periodic coordinates are folded into [0, 1)
and floating-point noise within 1e-10 of a boundary is snapped to zero before
generating any displacements. Atom order and lattice vectors are retained.
This matters because equivalent atomic positions in different cells require
different IFC lattice-translation labels. Changing CONTROL alone does not
convert existing force constants. For reused ShengBTE FC2, adjacent phono3py
unit-cell metadata is checked against the generated CONTROL convention;
inconsistent legacy inputs are rejected. Regenerate a matched set in a new
result directory, or explicitly remap the IFCs using their original metadata.
An explicitly supplied CONTROL retains the user's external convention.

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

tdbte: nonequilibrium dynamics
----------------------------------

For a new calculation, TD-BTE uses six sections. ``nepkappa run`` generates
matching FC2, FC3, and displacement metadata before building the collision
kernel and evolving the excitation. It does not calculate conductivity:

.. code-block:: yaml

   structure:
     poscar: examples/structures/3C-SiC/POSCAR_primitive
     relaxation: false
   calculator:
     name: nep
     model: potentials/3C-SiC/nep_3C-SiC.txt
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
``layout``, ``dpi``, ``path`` and ``tdbte-times`` for the six-panel overview;
run ``nepkappa plot`` after the calculation. See :doc:`tdbte` for the
energy-colored snapshots and branch temperatures. ``nepkappa tdbte`` with the same
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

plot: figures and analysis
------------------------------

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

For phono3py results, ``plot.tau: phph`` plots the three-phonon sum of
Normal and Umklapp processes, while ``nu`` plots N and U separately and
requires both ``gamma_N`` and ``gamma_U``. ``total`` compares ph-ph,
isotope, and their sum when ``gamma_isotope`` is present; without isotope
data it shows only ph-ph. ``isotope`` selects just the isotope channel;
``all`` adds any available N/U channels to the total comparison. The legacy
``normal`` and ``umklapp`` selections remain supported.
Both linewidth contributions are converted to rates with ``4*pi*gamma``;
temperature-independent isotope linewidths are broadcast across temperatures.
Saved HDF5 data determine the available channels, even if the current YAML
omits ``isotope`` or sets it to false. If YAML explicitly enables isotope
scattering but the result lacks ``gamma_isotope``, ``total``/``all`` refuse
to present an incomplete total; use ``phph`` or supply the missing result.
This plotted total excludes boundary scattering. FourPhonon's 3ph and 4ph
rates are overlaid as separately labelled series in one panel, without
summing them or adding isotope scattering. The separate image is named
``scattering_rate_3ph_4ph.png``. With the four harmonic panels and one
conductivity panel, this gives a 2-by-3 combined figure; optional N/U or
mode-resolved data may add further panels.

``plot.scattering-rate-min`` sets a display-only lower cutoff in ps^-1,
for example ``1.0e-8`` to hide extremely small rates on logarithmic axes.
The default ``0`` retains all positive rates. Points equal to the cutoff
are retained. The cutoff applies separately to each plotted channel after
conversion to scattering rates, including FourPhonon rate panels, and to
both separate and combined figures. It does not change saved linewidths,
channel sums or thermal conductivity. State any cutoff in figure captions.

For ordinary non-Wigner LBTE files containing ``kappa_RTA``, both ``kappa``
(LBTE, solid line) and ``kappa_RTA`` (RTA, dashed line) are plotted together in
``kappa.png`` and the combined figure. The tensors must have matching shapes
and finite values, and receive the same geometry correction. With
``plot.kappa: all``, this comparison shows the mean of the three diagonal
components for each method; ``x``, ``y``, or ``z`` selects the corresponding
component for both. Without ``kappa_RTA``, the existing single-solution plot
is retained. Scattering rates and cumulative conductivity are unchanged;
this does not substitute RTA mode contributions for LBTE ones.

Standalone SSCHA generates FC2 only when both transport and bubble
postprocessing are disabled. Enabling either also generates FC3.

Without transport enabled, SSCHA plotting overlays all temperatures requested
by ``sscha.temps`` in dispersion, DOS, and group velocity panels. Each heat
capacity point uses the FC2 and volume of that same temperature. All requested
SSCHA results must be complete. ``plot.temperature`` does not select these
overview data. ``kappa.temps`` sets only the fixed-volume harmonic heat-capacity
reference grid (for example ``[100, 1000, 50]``); ``kappa.mesh`` sets the integration mesh.
The four individual PNGs and ``combined.png`` are written to the top-level
``output.result-dir/plots`` directory. This heat capacity is an auxiliary
harmonic-model estimate, not the temperature derivative of full SSCHA free energy.

When initial harmonic FC2 and matching metadata are available, the overview
also displays harmonic references. Fixed-volume SSCHA uses one gray reference;
QHA+SSCHA uses temperature-colored dashed reference curves (open markers for
group velocities) at each corresponding QHA volume. In both routes, the heat-capacity
reference is a gray dashed curve evaluated with one fixed-volume harmonic FC2,
while the SSCHA values are colored points without connecting lines. For QHA+SSCHA,
this curve uses the completed QHA scan point with volume ratio 1, not the
temperature-dependent equilibrium volumes. If that point or its FC2/metadata
is unavailable, the curve is omitted with a message; no expanded-volume reference
is substituted. Reference cells and atomic mappings must match for the spectral overlays.
Reference FC2 comes from the SSCHA run's base directory (or
``sscha.initial-fc2`` for fixed-volume SSCHA). Plotting never computes forces.
With transport enabled, ``plot.temperature`` retains the single-temperature
selection for the standard transport plotting route.

Wigner plots and reports use a shared conductivity reader. In phono3py LBTE
files, ``kappa`` can contain the particle solution also stored as
``kappa_intra_exact``; the total is that solution plus ``kappa_inter``.
The reader checks the stored convention and does not substitute
``kappa_intra_RTA`` for the LBTE particle solution.

output: result files
------------------------

``output`` sets ``result-dir`` and optional ``progress``. Use a distinct result
directory for each material and parameter set so that outputs are not mixed.

Parallel execution
----------------------

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

Custom phono3py command sampling
------------------------------------

``kappa.command`` replaces the complete phono3py command. Sampling can be
written once, inside the command:

.. code-block:: yaml

   kappa:
     engine: phono3py
     command: >-
       phono3py phono3py_disp.yaml --br --nu --mesh 31 31 31
       --tmin 100 --tmax 1000 --tstep 50
       --isotope --boundary-mfp 1000000 --tt smm19
   output:
     result-dir: results

NEP-kappa reads three positive mesh dimensions from ``--mesh`` and temperatures
from ``--ts`` (alias ``--temperature``) or the complete
``--tmin/--tmax/--tstep`` range. ``--ts`` accepts a single temperature or an
increasing, equally spaced list; it is interpreted as actual temperature
points, not the YAML range convention. Quoted values and ``--option=value``
forms are supported. The inferred settings also apply to ``info`` and ``plot``.
If ``kappa.mesh`` or ``kappa.temps`` is also provided, it must agree with the
command; conflicts are errors. Other calculation flags are not appended to
the custom command.

Settings hidden in external scripts or phono3py configuration files are not
read automatically; specify the corresponding YAML mesh and temperatures
when using those mechanisms. The command runs in ``output.result-dir``
without a shell, so shell pipes and redirection are not interpreted.
