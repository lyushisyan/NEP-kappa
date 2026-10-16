Requirements
==============

The core dependencies below are installed by ``pip`` from ``pyproject.toml``.
See :doc:`installation` for environment setup.

- python>=3.10
- numpy
- scipy
- ase
- h5py
- matplotlib
- phonopy>=4.1
- phono3py>=4.0.1
- hiphive
- trainstation
- seekpath
- calorine
- tqdm
- PyYAML

Additional requirements
-------------------------

- a valid NEP model file is required when using the ``nep`` force backend
- VASP force calculations require a working VASP executable and POTCAR file or potential library
- VASP examples contain placeholder paths that must be configured before use
- bundled examples use paths relative to the repository root; custom inputs
  can use paths relative to another launch directory
- Slurm LBTE mode requires ``sbatch`` and a result directory on a shared filesystem
- generic ASE/plugin calculators require their provider package to be installed
  separately; those optional machine-learning frameworks are not core dependencies
- native MACE force evaluation requires the ``mace`` optional dependency:
  ``python -m pip install -e '.[mace]'``
- Phonopy SSCHA requires an ASE calculator that supplies both total energy and
  forces; the NEP backend through Calorine satisfies this requirement

- Thirdorder FC3 requires ``thirdorder_vasp.py`` or a configured executable
- Fourthorder FC4 requires ``Fourthorder_vasp.py`` or a configured executable
- four-phonon transport requires FourPhonon; its Wigner route needs Wigner_Park
- TD-BTE uses Numba acceleration automatically; Numba is installed through
  HiPhive's dependencies, so no additional installation step is needed
- documentation builds use Sphinx and the Read the Docs theme:
  ``python -m pip install -e '.[docs]'``

The numerical backend and MPI libraries used to compile external executables
must be available on the calculation host. Python installation does not install
VASP, Thirdorder, Fourthorder or FourPhonon.
