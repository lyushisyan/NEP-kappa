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
     --wigner true \
     --result_dir result \
     --output_name kappa

.. note::

   For input-file syntax and parameter details, see :doc:`input_files`.

Typical outputs
---------------

Depending on the selected workflow, typical output files include:

- ``result/run.log``
- ``result/POSCAR_relaxed``
- ``result/phono3py_disp.yaml``
- ``result/hiphive_model.fcp``
- ``result/fc2.hdf5``
- ``result/fc3.hdf5``
- ``result/kappa.hdf5`` or ``result/kappa-mXXXXX.hdf5``

By default, workflow outputs are written under ``result/``. The terminal output
is also saved to ``result/run.log``.
