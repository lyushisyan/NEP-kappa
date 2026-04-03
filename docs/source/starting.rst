Quick Start
===========

This page shows the minimal steps required to run **NEP-kappa**.

Prepare a structure
-------------------

Use the example structure builder:

.. code-block:: bash

   python Structure/build_and_relax_prim.py --model film --thick 8 --vac 15
   python Structure/build_and_relax_prim.py --model wire --thick 8 --vac 15
   python Structure/build_and_relax_prim.py --model bulk

This writes:

- ``Structure/POSCAR``

You may also use your own structure file and pass it later with ``--poscar``.

Run the workflow
----------------

Using an input file:

.. code-block:: bash

   python nepkappa.py input_1.txt

or

.. code-block:: bash

   python nepkappa.py input_2.txt

Using direct command-line arguments:

.. code-block:: bash

   python nepkappa.py \
     --poscar Structure/POSCAR \
     --nep_model NEP/Si_2025_Xuke.txt \
     --do_relax false \
     --dim 4 4 1 \
     --mesh 21 21 1 \
     --temps 100 1000 50 \
     --fc2fc3 true \
     --use_hiphive false \
     --method rta \
     --wigner true

.. note::

   For input-file syntax and parameter details, see :doc:`input_files`.

Typical outputs
---------------

Depending on the selected workflow, typical output files include:

- ``POSCAR_relaxed``
- ``phono3py_disp.yaml``
- ``hiphive_model.fcp``
- ``fc2.hdf5``
- ``fc3.hdf5``
- ``kappa-mXXXXX.hdf5``
