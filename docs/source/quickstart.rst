Quickstart
==========

Environment
-----------

Install doc/runtime dependencies according to your environment, then ensure
`phono3py` is available from shell.

Prepare structure
-----------------

Example builder script:

.. code-block:: bash

   python Structure/build_and_relax_prim.py --model film --thick 8 --vac 15
   python Structure/build_and_relax_prim.py --model wire --thick 8 --vac 15
   python Structure/build_and_relax_prim.py --model bulk

Run workflow
------------

With input files:

.. code-block:: bash

   python nepkappa.py input_1.txt
   python nepkappa.py input_2.txt

With CLI arguments:

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

Key points
----------

- `--temps` accepts only 1 value or 3 values.
- `--fc2fc3 false` means reuse existing force constants.
