Requirements
============

**NEP-kappa** requires the following Python packages:

- python>=3.9
- numpy
- ase
- h5py
- matplotlib
- phonopy
- phono3py>=4.0.1
- hiphive
- trainstation
- seekpath
- calorine
- tqdm
- PyYAML

Addition
-----------------------

- a valid NEP model file is required when using the ``nep`` force backend
- VASP force calculations require a working VASP executable and POTCAR file or potential library
- calculations should be run from the repository root directory
