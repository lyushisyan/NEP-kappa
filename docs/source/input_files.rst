Input Files
===========

NEP-kappa v1.1 uses YAML input files. Each YAML file describes one workflow by
grouping settings into the same stages used by the command line:

- ``structure``: input POSCAR
- ``calculator``: force backend, either NEP or VASP
- ``relaxation``: structure relaxation settings
- ``force-constant``: force-constant generation settings
- ``kappa``: thermal-conductivity settings passed to ``phono3py``
- ``plot``: plot layout, band path, relaxation-time channel, and kappa component
- ``output``: progress display and result directory

The ``compare`` command uses a smaller YAML file with ``reference``,
``candidate``, ``compare``, and ``plot`` sections.

Command behavior
----------------

.. code-block:: bash

   nepkappa relax input.yaml
   nepkappa fc input.yaml
   nepkappa kappa input.yaml
   nepkappa plot input.yaml
   nepkappa compare compare.yaml
   nepkappa run input.yaml
   nepkappa info input.yaml

- ``nepkappa relax`` relaxes the structure and writes ``POSCAR_relaxed`` to ``output.result_dir``.
- ``nepkappa fc`` generates ``phono3py_disp.yaml``, ``fc2.hdf5``, and ``fc3.hdf5``.
- ``nepkappa kappa`` computes thermal conductivity using existing ``phono3py_disp.yaml``, ``fc2.hdf5``, and ``fc3.hdf5``.
- ``nepkappa plot`` creates standard plots from ``fc2.hdf5`` and ``kappa-m*.hdf5``.
- ``nepkappa compare`` overlays DFT and NEP result directories in the same standard figures.
- ``nepkappa run`` executes ``relax``, ``fc``, and ``kappa`` in sequence.
- ``nepkappa info`` prints the parsed configuration without running a calculation.

When ``relaxation.enabled`` is ``true``, ``nepkappa fc`` reads
``POSCAR_relaxed`` from ``output.result_dir``. Run ``nepkappa relax`` first, or
use ``nepkappa run``.

Packaged examples
-----------------

The repository provides example POSCAR files and YAML files in ``examples/``:

- ``examples/POSCAR_bulk``: bulk Si structure
- ``examples/POSCAR_film``: Si film structure
- ``examples/1-bulk-nep-rta.yaml``: bulk, NEP forces, finite displacement, RTA
- ``examples/2-bulk-nep-hiphive-rta.yaml``: bulk, NEP forces, HiPhive, RTA
- ``examples/3-bulk-nep-lbte.yaml``: bulk, NEP forces, finite displacement, LBTE
- ``examples/4-bulk-nep-rta-wigner.yaml``: bulk, NEP forces, finite displacement, Wigner transport
- ``examples/5-bulk-vasp-rta.yaml``: bulk, VASP relaxation and VASP forces, finite displacement, RTA
- ``examples/6-bulk-vasp-hiphive-rta.yaml``: bulk, VASP relaxation and VASP forces, HiPhive, RTA
- ``examples/7-film-nep-rta.yaml``: film, NEP forces, finite displacement, RTA
- ``examples/8-film-nep-hiphive-rta.yaml``: film, NEP forces, HiPhive, RTA
- ``examples/compare.yaml``: DFT-vs-NEP comparison plotting

Example outputs are written to ``results/...`` and are not tracked by Git.

Minimal NEP example
-------------------

.. code-block:: yaml

   structure:
     poscar: examples/POSCAR_bulk
     dimensionality: 3

   calculator:
     name: nep
     nep_model: potentials/Si_Bulk_Fan.txt

   relaxation:
     enabled: true

   force-constant:
     dim-fc2: [3, 3, 3]
     dim-fc3: [3, 3, 3]
     use_hiphive: false

   kappa:
     mesh: [21, 21, 21]
     temps: [100, 1000, 50]
     method: rta
     isotope: false
     bfmp: 1.0e6
     wigner: false

   plot:
     layout: separate
     path: seekpath
     tau: total
     kappa: all
     temperature: 300
     dpi: 300

   output:
     progress: true
     result_dir: results/1-bulk-nep-rta

HiPhive example
---------------

.. code-block:: yaml

   structure:
     poscar: examples/POSCAR_film
     dimensionality: 2
     effective_thickness: 10.0
     vacuum_axis: z

   calculator:
     name: nep
     nep_model: potentials/Si_NWs_XuKe.txt

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
     result_dir: results/7-film-nep-hiphive-rta

VASP example
------------

.. code-block:: yaml

   structure:
     poscar: examples/POSCAR_bulk

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

   kappa:
     mesh: [21, 21, 21]
     temps: [100, 1000, 50]
     method: rta
     isotope: false
     bfmp: 1.0e6
     wigner: false

   output:
     progress: true
     result_dir: results/5-bulk-vasp-rta

``structure``
-------------

``poscar`` is the input POSCAR path. Relative paths are interpreted from the
directory where ``nepkappa`` is launched, usually the repository root.

.. code-block:: yaml

   structure:
     poscar: examples/POSCAR_bulk
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
--------------

For NEP calculations:

.. code-block:: yaml

   calculator:
     name: nep
     nep_model: potentials/Si_Bulk_Fan.txt

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

``relaxation``
--------------

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
------------------

Finite-displacement route:

.. code-block:: yaml

   force-constant:
     dim-fc2: [3, 3, 3]
     dim-fc3: [3, 3, 3]
     use_hiphive: false

``dim-fc2`` sets the phonon supercell used for FC2. ``dim-fc3`` sets the
supercell used for FC3. The deprecated ``dim`` key is still accepted for
compatibility and is interpreted as both ``dim-fc2`` and ``dim-fc3``.

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

``kappa``
---------

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
true`` enables Wigner transport through ``phono3py-wte``.

``plot``
--------

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
- ``combined``: write one 2-by-3 multi-panel ``combined.png`` figure
- ``both``: write separate PNG files and ``combined.png``

NEP-kappa always generates the six standard plots: ``dispersion``, ``dos``,
``heat_capacity``, ``group_velocity``, ``relaxation_time``, and ``kappa``.

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
- ``all``: plot total, N, and U channels together when available

``kappa`` controls the thermal-conductivity components shown in the kappa
figure:

- ``x``: plot only ``kappa_xx``
- ``y``: plot only ``kappa_yy``
- ``z``: plot only ``kappa_zz``
- ``all``: plot ``kappa_xx``, ``kappa_yy``, ``kappa_zz``, and their average

``temperature`` selects the target temperature for the relaxation-time plot.
NEP-kappa uses the closest temperature available in ``kappa-m*.hdf5``.

``output``
----------

.. code-block:: yaml

   output:
     progress: true
     result_dir: results/1-bulk-nep-rta

All generated files are written inside ``result_dir``. The terminal output is
also saved to ``run.log`` in that directory.

``plot`` output
---------------

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
- ``kappa.png``: selected thermal conductivity component or all diagonal components plus average
- ``combined.png``: 2-by-3 combined multi-panel figure when ``layout`` is ``combined`` or ``both``

DFT-vs-NEP comparison
---------------------

``nepkappa compare compare.yaml`` reads one completed DFT result directory and
one completed NEP result directory, then overlays them in the same six standard
figures. Each result directory must contain ``phono3py_disp.yaml``,
``fc2.hdf5``, and ``kappa-m*.hdf5``.

.. code-block:: yaml

   reference:
     dft_dir: results/dft
     label: DFT

   candidate:
     nep_dir: results/nep
     label: NEP

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
