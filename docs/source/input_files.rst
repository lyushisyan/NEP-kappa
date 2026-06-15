Input Files
===========

NEP-kappa v1.1 uses YAML input files. Each YAML file describes one workflow by
grouping settings into the same stages used by the command line:

- ``structure``: input POSCAR
- ``calculator``: force backend, either NEP or VASP
- ``relaxation``: structure relaxation settings
- ``force-constant``: force-constant generation settings
- ``kappa``: thermal-conductivity settings passed to ``phono3py``
- ``output``: progress display and result directory

Command behavior
----------------

.. code-block:: bash

   nepkappa relax input.yaml
   nepkappa fc input.yaml
   nepkappa kappa input.yaml
   nepkappa run input.yaml
   nepkappa info input.yaml

- ``nepkappa relax`` relaxes the structure and writes ``POSCAR_relaxed`` to ``output.result_dir``.
- ``nepkappa fc`` generates ``phono3py_disp.yaml``, ``fc2.hdf5``, and ``fc3.hdf5``.
- ``nepkappa kappa`` computes thermal conductivity using existing ``phono3py_disp.yaml``, ``fc2.hdf5``, and ``fc3.hdf5``.
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
- ``examples/7-film-nep-hiphive-rta.yaml``: film, NEP forces, HiPhive, RTA

Example outputs are written to ``results/...`` and are not tracked by Git.

Minimal NEP example
-------------------

.. code-block:: yaml

   structure:
     poscar: examples/POSCAR_bulk

   calculator:
     name: nep
     nep_model: potentials/Si_Bulk_Fan.txt

   relaxation:
     enabled: true

   force-constant:
     dim: [3, 3, 3]
     fc_calculator: symfc
     use_hiphive: false

   kappa:
     mesh: [21, 21, 21]
     temps: [100, 1000, 50]
     method: rta
     wigner: false

   output:
     progress: true
     result_dir: results/1-bulk-nep-rta

HiPhive example
---------------

.. code-block:: yaml

   structure:
     poscar: examples/POSCAR_film

   calculator:
     name: nep
     nep_model: potentials/Si_NWs_XuKe.txt

   relaxation:
     enabled: false

   force-constant:
     dim: [4, 4, 1]
     use_hiphive: true
     n_structures: 500
     rattle_std: 0.03
     min_dist: 2.2
     cutoffs: [5.0, 4.0]

   kappa:
     mesh: [21, 21, 1]
     temps: [100, 1000, 50]
     method: rta
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
     dim: [3, 3, 3]
     fc_calculator: symfc
     use_hiphive: false

   kappa:
     mesh: [21, 21, 21]
     temps: [100, 1000, 50]
     method: rta
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
     dim: [3, 3, 3]
     fc_calculator: symfc
     use_hiphive: false

HiPhive route:

.. code-block:: yaml

   force-constant:
     dim: [3, 3, 3]
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
     wigner: false

``temps`` accepts either one temperature or ``[tmin, tmax, tstep]``.
``method`` can be ``rta`` or ``lbte``. ``wigner: true`` enables Wigner transport
through ``phono3py-wte``.

``output``
----------

.. code-block:: yaml

   output:
     progress: true
     result_dir: results/1-bulk-nep-rta

All generated files are written inside ``result_dir``. The terminal output is
also saved to ``run.log`` in that directory.
