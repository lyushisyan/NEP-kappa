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

Supported routes
----------------

Two force-constant routes are supported:

- **Finite Displacement**
- **HiPhive fitting**

Repository layout
-----------------

- ``pyproject.toml``:
  package metadata and install configuration
- ``src/nepkappa/``:
  installable Python package
- ``Structure/build_and_relax_prim.py``:
  example structure builder for Si systems
- ``examples-input/input_bulk-rta.yaml``:
  bulk Si example using finite displacement
- ``examples-input/input_film-rta.yaml``:
  film Si example using HiPhive
- ``examples-input/``:
  structures used by the packaged YAML examples

Notes
-----

- The structure builder is an example script for Si systems.
- For other materials, prepare your own structure file and pass it with ``--poscar``.
- Example inputs are kept under ``examples-input/`` so the package root stays clean.
