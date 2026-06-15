Tutorial
========

This tutorial shows a minimal bulk Si calculation using the current YAML
framework.

Step 1: Choose the input structure
----------------------------------

Use the packaged bulk POSCAR:

.. code-block:: text

   examples-input/POSCAR_bulk

For your own calculation, replace this with any POSCAR path in the YAML
``structure.poscar`` field.

Step 2: Create a YAML input
---------------------------

Create ``input.yaml``:

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
     result_dir: result/bulk-nep-rta

Step 3: Run the workflow
------------------------

Run everything:

.. code-block:: bash

   nepkappa run input.yaml

Or run the stages separately:

.. code-block:: bash

   nepkappa relax input.yaml
   nepkappa fc input.yaml
   nepkappa kappa input.yaml

You can inspect the parsed settings first:

.. code-block:: bash

   nepkappa info input.yaml

Step 4: Check outputs
---------------------

The calculation writes all generated files to ``result/bulk-nep-rta``:

- ``run.log``
- ``POSCAR_relaxed``
- ``phono3py_disp.yaml``
- ``fc2.hdf5``
- ``fc3.hdf5``
- ``kappa-m*.hdf5``

For details on reading ``kappa-m*.hdf5`` files, see the
`phono3py HDF5 documentation <https://phonopy.github.io/phono3py/hdf5_howto.html>`_.
