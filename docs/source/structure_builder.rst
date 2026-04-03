Structure Builder (Example)
===========================

The script `Structure/build_and_relax_prim.py` is an example implementation.

Supported models
----------------

- `film`: Si(001) slab
- `wire`: 1D Si nanowire
- `bulk`: Si primitive cell

Notes
-----

- Output is always written to `Structure/POSCAR`.
- Default NEP model: `NEP/Si_2025_Xuke.txt`.
- In `bulk` mode, `--thick` and `--vac` are ignored.
- For non-example structures, place your structure files in `Structure/`
  and pass the path to `--poscar` when running `nepkappa.py`.
