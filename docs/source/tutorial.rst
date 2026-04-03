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

Create `input.txt`:

.. code-block:: text

   --poscar      Structure/POSCAR
   --nep_model   NEP/Si_2025_Xuke.txt
   --do_relax    false
   --dim         4 4 1
   --mesh        21 21 21
   --temps       100 1000 100
   --fc2fc3      true
   --use_hiphive false
   --method      rta
   --wigner      true

Step 3: Run workflow
--------------------

.. code-block:: bash

   python nepkappa.py input.txt

Typical outputs:

- `fc2.hdf5`
- `fc3.hdf5`
- `kappa-m11111.hdf5` (depends on mesh)

Step 4: Plot examples
---------------------

.. code-block:: bash

   python Example-Resuts/Film-1nm/plot_film.py

This generates:

- `Example-Resuts/Film-1nm/Si_film_1nm_4panel.pdf`

Notes
-----

- Example plotting scripts are templates for your own data.
- For your own results, copy and edit script paths accordingly.
