Input Files
===========

**NEP-kappa** supports three input modes:

- YAML input-file mode, recommended for new runs
- legacy text input-file mode
- command-line mode

The legacy text option names are the same as command-line options.

Common commands:

.. code-block:: bash

   nepkappa info input.yaml
   nepkappa fc input.yaml
   nepkappa kappa input.yaml
   nepkappa run input.yaml

YAML input file format
----------------------

YAML input is recommended because it groups structure, force-constant, kappa,
and output settings clearly.

Example ``result_dir`` paths such as ``examples-output/...`` are local run
artifacts and are not tracked by the repository.

Bulk Si + finite displacement + phono3py

.. code-block:: yaml

   structure:
     poscar: examples-input/POSCAR_bulk
     nep_model: NEP/Si_Bulk_Fan.txt
     relax: false

   force_constants:
     dim: [3, 3, 3]
     use_hiphive: false

   kappa:
     mesh: [21, 21, 21]
     temps: [100, 1000, 50]
     method: rta
     wigner: false

   output:
     progress: true
     result_dir: examples-output/bulk-fd-rta

Film Si + HiPhive + phono3py

.. code-block:: yaml

   structure:
     poscar: examples-input/POSCAR_film
     nep_model: NEP/Si_NWs_XuKe.txt
     relax: false

   force_constants:
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
     result_dir: examples-output/film-hiphive-rta

Bulk Si + VASP forces + finite displacement

.. code-block:: yaml

   structure:
     poscar: examples-input/POSCAR_bulk
     relax: false

   calculator:
     name: vasp
     vasp_path: /root/software/vasp.6.4.3/bin/vasp_std
     potcar_path: /root/software/potpaw_PBE.64
     vasp_workdir: vasp-runs
     vasp_kwargs:
       encut: 650
       kspacing: 0.2
       kgamma: true
       kpar: 2
       ncore: 3

   force_constants:
     dim: [3, 3, 3]
     use_hiphive: false

   kappa:
     mesh: [21, 21, 21]
     temps: [100, 1000, 50]
     method: rta
     wigner: false

   output:
     progress: true
     result_dir: examples-output/bulk-vasp-rta

YAML aliases
------------

- ``structure.relax`` maps to ``--do_relax``
- ``calculator.name`` maps to ``--calculator``
- ``kappa`` and ``transport`` are both accepted for transport settings

Legacy text input format
------------------------

Text input files remain supported for compatibility with existing runs.

Example:

.. code-block:: text

   --poscar      examples-input/POSCAR_film
   --nep_model   NEP/Si_NWs_XuKe.txt
   --calculator  nep
   --do_relax    false
   --dim         4 4 1
   --use_hiphive false
   --mesh        21 21 1
   --temps       100 1000 50
   --method      rta
   --wigner      false
   --progress    true
   --result_dir  result

Rules
-----

- In YAML input files, use ``.yaml`` or ``.yml`` extensions.
- In legacy text input files, one option per line is recommended.
- Legacy text values are separated by spaces.
- In legacy text files, ``#`` starts a comment.

Parameter reference
-------------------

``--poscar``
^^^^^^^^^^^^

Path to the input structure file.

**Type:** string
**Default:** ``POSCAR``

Example:

.. code-block:: text

   --poscar Structure/POSCAR

``--nep_model``
^^^^^^^^^^^^^^^

Path to the NEP model file.

**Type:** string
**Required:** yes when ``--calculator nep``

Example:

.. code-block:: text

   --nep_model NEP/Si_NWs_XuKe.txt

``--calculator``
^^^^^^^^^^^^^^^^

Force backend used by ``nepkappa fc``.

**Type:** string
**Choices:** ``nep``, ``vasp``
**Default:** ``nep``

- ``nep``: use ``calorine.calculators.CPUNEP``
- ``vasp``: write VASP input files, run VASP, and read forces from
  ``vasprun.xml`` or ``OUTCAR`` in one directory per structure

``--vasp_path``
^^^^^^^^^^^^^^^

Path to the VASP executable.

**Type:** string
**Default:** auto-detect common paths such as
``/root/software/vasp.6.4.3/bin/vasp_std``

Example:

.. code-block:: text

   --vasp_path /root/software/vasp.6.4.3/bin/vasp_std

``--potcar_path``
^^^^^^^^^^^^^^^^^

Path to a POTCAR file or a potential-library directory.

**Type:** string
**Default:** auto-detect common paths such as
``/root/software/potpaw_PBE.64``

If a directory is given, NEP-kappa tries to assemble a POTCAR using the ordered
elements in the POSCAR. For multi-element POSCAR files, element POTCAR chunks
are concatenated in POSCAR element order.

``--vasp_command``
^^^^^^^^^^^^^^^^^^

Optional full command used to run VASP. This takes precedence over
``--vasp_path`` and is useful for MPI launchers.

**Type:** string
**Default:** use ``--vasp_path`` or auto-detection

Example:

.. code-block:: text

   --vasp_command "mpirun -np 64 /root/software/vasp.6.4.3/bin/vasp_std"

``--vasp_workdir``
^^^^^^^^^^^^^^^^^^

Subdirectory under ``result_dir`` for VASP force calculations.

**Type:** string
**Default:** ``vasp-runs``

NEP-kappa writes force calculations to paths such as
``result/vasp-runs/fc2/00001`` and ``result/vasp-runs/fc3/00001``.

``--vasp_kwargs``
^^^^^^^^^^^^^^^^^

JSON object with INCAR/KPOINTS options. Common VASP single-point defaults are
filled by NEP-kappa, so the YAML usually only needs project-specific settings
such as cutoff, k-point spacing, and parallel parameters.

**Type:** JSON object
**Default:** ``{}``

Example:

.. code-block:: text

   --vasp_kwargs '{"encut": 650, "kspacing": 0.2, "kgamma": true, "kpar": 2, "ncore": 3}'

``--do_relax``
^^^^^^^^^^^^^^

Whether to relax the structure before force-constant generation.

**Type:** boolean
**Default:** ``false``

Example:

.. code-block:: text

   --do_relax true

``--dim``
^^^^^^^^^

Supercell dimensions for force-constant calculations.

**Type:** three integers
**Default:** ``4 4 1``

Example:

.. code-block:: text

   --dim 4 4 1

``--use_hiphive``
^^^^^^^^^^^^^^^^^

Select the force-constant generation route.

**Type:** boolean
**Default:** ``false``

- ``false``: finite displacement
- ``true``: HiPhive

``--n_structures``
^^^^^^^^^^^^^^^^^^

Number of rattled structures used in HiPhive fitting.

**Type:** integer
**Default:** ``50``

Example:

.. code-block:: text

   --n_structures 100

``--rattle_std``
^^^^^^^^^^^^^^^^

Standard deviation of random displacements in HiPhive.

**Type:** float
**Default:** ``0.03``

Example:

.. code-block:: text

   --rattle_std 0.03

``--cutoffs``
^^^^^^^^^^^^^

Cutoff distances used in HiPhive cluster-space construction.

**Type:** one or more floats
**Default:** ``5.0``

Examples:

.. code-block:: text

   --cutoffs 5.0

.. code-block:: text

   --cutoffs 5.0 4.0

``--min_dist``
^^^^^^^^^^^^^^

Minimum allowed interatomic distance for rattled structures.

**Type:** float
**Default:** ``2.0``

Example:

.. code-block:: text

   --min_dist 2.0

``--mesh``
^^^^^^^^^^

q-point mesh used by phono3py.

**Type:** three integers
**Default:** ``21 21 1``

Example:

.. code-block:: text

   --mesh 21 21 1

``--temps``
^^^^^^^^^^^

Temperature setting.

Accepted formats:

Single temperature:

.. code-block:: text

   --temps 300

Temperature range:

.. code-block:: text

   --temps 100 1000 50

Rule:

- one value: single temperature
- three values: ``tmin tmax tstep``

``--method``
^^^^^^^^^^^^

Thermal conductivity solution method.

**Type:** string
**Choices:** ``lbte``, ``rta``
**Default:** ``lbte``

``--wigner``
^^^^^^^^^^^^

Whether to enable Wigner transport through the ``phono3py-wte`` plugin.
When this is ``true``, NEP-kappa passes ``--tt wte`` to phono3py.

**Type:** boolean
**Default:** ``false``

``--progress``
^^^^^^^^^^^^^^

Whether to show progress bars, ETA, and stage-level timing summaries.

**Type:** boolean
**Default:** ``true``

``--result_dir``
^^^^^^^^^^^^^^^^

Directory for generated outputs and ``run.log``.

**Type:** string
**Default:** ``result``

Generated files such as ``fc2.hdf5``, ``fc3.hdf5``, ``phono3py_disp.yaml``, and
the default ``kappa-m{mesh}.hdf5`` file are written under this directory.
