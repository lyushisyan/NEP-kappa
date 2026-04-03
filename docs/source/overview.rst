Overview
========

Workflow
--------

1. Build or prepare a structure.
2. Compute force constants via:

   - Finite Displacement (`--use_hiphive false`)
   - HiPhive fitting (`--use_hiphive true`)

3. Run `phono3py` to compute thermal conductivity.
4. Plot example results from the `Example-Resuts/` folder.

Repository layout
-----------------

- `nepkappa.py`: CLI entry point
- `workflow.py`: main workflow logic
- `Structure/build_and_relax_prim.py`: example structure builder script
- `Example-Resuts/`: example results and example plotting scripts

Notes
-----

- `Structure/build_and_relax_prim.py` is an example script.
- If you have other structures, place them directly under `Structure/` and pass them via `--poscar`.
- `Example-Resuts/` plotting scripts are examples that you can copy/modify.
