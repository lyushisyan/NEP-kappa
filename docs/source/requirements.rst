Requirements
==============

**NEP-kappa** requires the following Python packages:

- python>=3.9
- numpy
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
