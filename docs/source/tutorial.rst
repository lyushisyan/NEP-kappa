Tutorial
========

This is a minimal end-to-end example (film, finite displacement, RTA).

Step 1: Build a film structure
------------------------------

.. code-block:: bash

   python Structure/build_and_relax_prim.py --model film --thick 8 --vac 15

This writes:

- `Structure/POSCAR`

Step 2: Create an input file
----------------------------

Create `input.yaml`:

.. code-block:: yaml

   structure:
     poscar: examples-input/POSCAR_film
     nep_model: NEP/Si_NWs_XuKe.txt
     relax: false

   force_constants:
     dim: [4, 4, 1]
     use_hiphive: false

   kappa:
     mesh: [21, 21, 1]
     temps: [100, 1000, 50]
     method: rta
     wigner: false

   output:
     result_dir: result

Step 3: Run force constants and kappa
-------------------------------------

.. code-block:: bash

   nepkappa fc input.yaml
   nepkappa kappa input.yaml

For a one-command run, use:

.. code-block:: bash

   nepkappa run input.yaml

You can check the parsed settings first:

.. code-block:: bash

   nepkappa info input.yaml

Typical outputs:

Generated outputs are local run artifacts and are not tracked by the
repository.

- `result/run.log`
- `result/fc2.hdf5`
- `result/fc3.hdf5`
- `result/kappa-m*.hdf5`

For details on reading ``kappa-m*.hdf5`` files, please refer to the
`phono3py HDF5 documentation <https://phonopy.github.io/phono3py/hdf5_howto.html>`_.

Notes
-----

- For your own results, place input structures where convenient and update
  ``structure.poscar`` accordingly.
