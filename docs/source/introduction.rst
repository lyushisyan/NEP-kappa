Introduction
============

**NEP-kappa** is a workflow for lattice thermal conductivity calculations.

It combines:

- **NEP** for structure relaxation and force evaluation
- **Finite Displacement** or **HiPhive** for force-constant generation
- **phono3py** for thermal conductivity calculations

Workflow
--------

The standard workflow is:

1. Prepare a structure
2. Relax the structure with NEP (optional)
3. Generate force constants
4. Run ``phono3py`` to compute thermal conductivity
5. Plot example results

Supported routes
----------------

Two force-constant routes are supported:

- **Finite Displacement**
- **HiPhive fitting**

Repository layout
-----------------

- ``nepkappa.py``:
  workflow entry point
- ``workflow.py``:
  core workflow implementation
- ``Structure/build_and_relax_prim.py``:
  example structure builder for Si systems
- ``input_1.txt``:
  example input for finite displacement
- ``input_2.txt``:
  example input for HiPhive
- ``Example-Resuts/``:
  example results and plotting scripts

Notes
-----

- The structure builder is an example script for Si systems.
- For other materials, prepare your own structure file and pass it with ``--poscar``.
- Plotting scripts in ``Example-Resuts/`` are examples only.