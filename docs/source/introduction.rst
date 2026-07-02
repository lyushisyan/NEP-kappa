Introduction
============

**NEP-kappa** is an installable workflow package for lattice thermal
conductivity calculations. It is designed around the same calculation stages
used in production runs:

1. optional structure relaxation
2. force-constant generation
3. thermal conductivity calculation with ``phono3py``
4. plotting phonon and thermal-transport results from existing HDF5 outputs

The commands are:

- ``nepkappa relax input.yaml`` relaxes the structure
- ``nepkappa fc input.yaml`` generates ``fc2.hdf5``, ``fc3.hdf5``, and ``phono3py_disp.yaml``
- ``nepkappa kappa input.yaml`` computes thermal conductivity from existing force constants
- ``nepkappa plot input.yaml`` plots phonon and thermal-transport results
- ``nepkappa run input.yaml`` runs ``relax``, ``fc``, and ``kappa`` in sequence
- ``nepkappa info input.yaml`` prints the parsed configuration without running

Calculator backends
-------------------

Two force backends are supported:

- ``nep``: use a NEP model through ``calorine.calculators.CPUNEP``
- ``vasp``: write VASP inputs, run VASP, and read forces from VASP outputs

Force-constant routes
---------------------

NEP-kappa can generate force constants by:

- finite displacement with ``phono3py``
- random structure generation and fitting with ``HiPhive``

Plotting outputs
----------------

``nepkappa plot`` reads ``phono3py_disp.yaml``, ``fc2.hdf5``, and
``kappa-m*.hdf5`` from ``output.result_dir``. It generates phonon dispersion,
DOS, volume heat capacity, group velocity in km/s, relaxation time, and thermal
conductivity figures. The plot layout, high-symmetry path, relaxation-time
channel, and thermal-conductivity component are controlled by the YAML
``plot`` section.

Repository layout
-----------------

- ``pyproject.toml``: package metadata and dependencies
- ``src/nepkappa/``: Python package and command-line implementation
- ``examples/``: repository-provided POSCAR and YAML examples
- ``potentials/``: example NEP model files
- ``docs/``: Sphinx documentation
- ``tests/``: lightweight configuration and CLI tests

Run outputs are written to the YAML ``output.result_dir`` value, commonly
``results/...``. These generated files are local artifacts and are not
tracked by the repository.
