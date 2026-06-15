Input Files
===========

NEP-kappa v1.1 uses YAML input files. The YAML is organized by workflow stage:

- ``structure``: input structure
- ``calculator``: NEP or VASP force backend
- ``relaxation``: optional structure relaxation
- ``force-constant``: finite-displacement or HiPhive force constants
- ``kappa``: phono3py thermal-conductivity settings
- ``output``: progress display and result directory

Commands
--------

.. code-block:: bash

   nepkappa relax input.yaml
   nepkappa fc input.yaml
   nepkappa kappa input.yaml
   nepkappa run input.yaml
   nepkappa info input.yaml

- ``nepkappa relax`` prepares or relaxes the input structure.
- ``nepkappa fc`` generates ``fc2.hdf5``, ``fc3.hdf5``, and ``phono3py_disp.yaml``.
- ``nepkappa kappa`` computes thermal conductivity from existing force constants.
- ``nepkappa run`` executes ``relax``, ``fc``, and ``kappa`` in sequence.
- ``nepkappa info`` prints the parsed configuration without running.

NEP finite-displacement example
-------------------------------

.. code-block:: yaml

   structure:
     poscar: examples-input/POSCAR_bulk

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
     result_dir: examples-output/1-bulk-nep-rta

NEP HiPhive example
-------------------

.. code-block:: yaml

   structure:
     poscar: examples-input/POSCAR_film

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
     result_dir: examples-output/7-film-nep-hiphive-rta

VASP finite-displacement example
--------------------------------

.. code-block:: yaml

   structure:
     poscar: examples-input/POSCAR_bulk

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
     result_dir: examples-output/5-bulk-vasp-rta

Packaged examples
-----------------

- ``input_1-bulk-nep-rta.yaml``: bulk, NEP forces, finite displacement, RTA
- ``input_2-bulk-nep-hiphive-rta.yaml``: bulk, NEP forces, HiPhive, RTA
- ``input_3-bulk-nep-lbte.yaml``: bulk, NEP forces, finite displacement, LBTE
- ``input_4-bulk-nep-rta-wigner.yaml``: bulk, NEP forces, finite displacement, Wigner transport
- ``input_5-bulk-vasp-rta.yaml``: bulk, VASP relaxation and VASP forces, finite displacement, RTA
- ``input_6-bulk-vasp-hiphive-rta.yaml``: bulk, VASP relaxation and VASP forces, HiPhive, RTA
- ``input_7-film-nep-hiphive-rta.yaml``: film, NEP forces, HiPhive, RTA

Section reference
-----------------

``structure``
^^^^^^^^^^^^^

``poscar`` is the POSCAR-style structure file.

.. code-block:: yaml

   structure:
     poscar: examples-input/POSCAR_bulk

``calculator``
^^^^^^^^^^^^^^

Use ``name: nep`` for NEP forces:

.. code-block:: yaml

   calculator:
     name: nep
     nep_model: potentials/Si_Bulk_Fan.txt

Use ``name: vasp`` for VASP forces:

.. code-block:: yaml

   calculator:
     name: vasp
     vasp_command: mpirun -np 24 /path/to/vasp_std
     vasp_path: /path/to/vasp_std
     potcar_path: /path/to/potpaw_PBE.64

``potcar_path`` may point to a ready POTCAR file or to a POTCAR library
directory. For multi-element POSCAR files, NEP-kappa concatenates POTCAR chunks
in POSCAR element order.

``relaxation``
^^^^^^^^^^^^^^

For NEP runs, relaxation uses the NEP calculator. For VASP runs, NEP-kappa can
run staged VASP relaxation and pass each stage's relaxed structure to the next
stage.

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

``force-constant``
^^^^^^^^^^^^^^^^^^

Common options:

- ``dim``: supercell dimensions
- ``fc_calculator``: finite-displacement solver, normally ``symfc``
- ``use_hiphive``: use HiPhive instead of direct finite displacement
- ``n_structures``: number of rattled structures for HiPhive
- ``rattle_std``: random-displacement standard deviation for HiPhive
- ``min_dist``: minimum allowed distance in rattled structures
- ``cutoffs``: HiPhive cluster-space cutoffs
- ``workdir``: VASP force-calculation subdirectory
- ``vasp_kwargs``: important VASP single-point INCAR settings

NEP-kappa supplies default VASP single-point settings, so YAML files usually
only need calculation-specific values such as ``encut``, ``kspacing``, ``kpar``,
``ncore``, and ``lreal``.

``kappa``
^^^^^^^^^

Common options:

- ``mesh``: q-point mesh
- ``temps``: one temperature or ``[tmin, tmax, tstep]``
- ``method``: ``rta`` or ``lbte``
- ``wigner``: enable Wigner transport through ``phono3py-wte``

``output``
^^^^^^^^^^

``result_dir`` controls where all generated files are written. ``run.log`` is
saved in the same directory.

.. code-block:: yaml

   output:
     progress: true
     result_dir: examples-output/1-bulk-nep-rta
