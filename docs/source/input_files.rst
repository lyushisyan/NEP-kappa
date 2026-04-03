Input Files
===========

**NEP-kappa** supports two input modes:

- input-file mode
- command-line mode

The option names are the same in both modes.

Input file format
-----------------

Example:

NEP + phono3py

.. code-block:: text

   --poscar      Structure/POSCAR
   --nep_model   NEP/Si_2025_Xuke.txt
   --do_relax    false
   --dim         4 4 1
   --fc2fc3      true
   --use_hiphive false
   --mesh        21 21 1
   --temps       100 1000 50
   --method      rta
   --wigner      true

NEP + HiPhive + phono3py

.. code-block:: text

   --poscar      Structure/POSCAR
   --nep_model   NEP/Si_2025_Xuke.txt
   --do_relax    false
   --dim         4 4 1
   --fc2fc3      true
   --use_hiphive false
   --mesh        21 21 1
   --temps       100 1000 50
   --method      rta
   --wigner      true
   --n_structures 500
   --rattle_std   0.03
   --min_dist     2.2
   --cutoffs      5.0 4.0

Rules
-----

- One option per line is recommended
- Values are separated by spaces
- ``#`` starts a comment

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
**Required:** yes

Example:

.. code-block:: text

   --nep_model NEP/Si_2025_Xuke.txt

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

``--fc2fc3``
^^^^^^^^^^^^

Whether to recompute force constants.

**Type:** boolean  
**Default:** ``false``

- ``true``: recompute ``fc2.hdf5`` and ``fc3.hdf5``
- ``false``: reuse existing force constants

``--method``
^^^^^^^^^^^^

Thermal conductivity solution method.

**Type:** string  
**Choices:** ``lbte``, ``rta``  
**Default:** ``lbte``

``--wigner``
^^^^^^^^^^^^

Whether to enable the Wigner option in phono3py.

**Type:** boolean  
**Default:** ``false``