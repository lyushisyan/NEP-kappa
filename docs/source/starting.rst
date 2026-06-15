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

Using a recommended YAML input file for the full workflow:

.. code-block:: bash

   nepkappa run examples-input/input_bulk-rta.yaml

or

.. code-block:: bash

   nepkappa run examples-input/input_film-rta.yaml

You can also run the two stages separately:

.. code-block:: bash

   nepkappa fc examples-input/input_bulk-rta.yaml
   nepkappa kappa examples-input/input_bulk-rta.yaml

For VASP forces, edit ``examples-input/input_bulk-vasp-rta.yaml`` for your VASP
executable path, POTCAR path, and calculator settings, then run:

.. code-block:: bash

   nepkappa fc examples-input/input_bulk-vasp-rta.yaml
   nepkappa kappa examples-input/input_bulk-vasp-rta.yaml

You can inspect an input before running an expensive calculation:

.. code-block:: bash

   nepkappa info examples-input/input_bulk-rta.yaml

Using direct command-line arguments:

.. code-block:: bash

   nepkappa run \
     --poscar examples-input/POSCAR_film \
     --nep_model NEP/Si_NWs_XuKe.txt \
     --do_relax false \
     --dim 4 4 1 \
     --mesh 21 21 1 \
     --temps 100 1000 50 \
     --use_hiphive false \
     --method rta \
     --wigner false \
     --result_dir result

.. note::

   For input-file syntax and parameter details, see :doc:`input_files`.

Typical outputs
---------------

Generated example outputs are local run artifacts. They are written to the
configured ``result_dir`` such as ``examples-output/...`` and are not tracked
by the repository.

Depending on the selected workflow, typical output files include:

- ``result/run.log``
- ``result/POSCAR_relaxed``
- ``result/phono3py_disp.yaml``
- ``result/hiphive_model.fcp``
- ``result/fc2.hdf5``
- ``result/fc3.hdf5``
- ``result/vasp-runs/`` when using the VASP force backend
- ``result/kappa.hdf5`` or ``result/kappa-mXXXXX.hdf5``

By default, workflow outputs are written under ``result/``. The terminal output
is also saved to ``result/run.log``.
