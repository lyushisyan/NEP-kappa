Installation
==============

Install in an isolated Python environment
-------------------------------------------

The package declares Python >=3.9; Python 3.11 is a practical starting choice.
Dependency availability may vary with the platform.

.. code-block:: bash

   git clone https://github.com/lyushisyan/NEP-kappa.git
   cd NEP-kappa
   python -m venv .venv
   source .venv/bin/activate
   python -m pip install -e .
   nepkappa --version
   nepkappa validate examples/nep-rta.yaml

The editable installation uses this checkout directly. After local source
updates, the command sees them; reinstall when dependency declarations change.
Python libraries for the standard NEP/Phonopy/phono3py workflow are installed
from ``pyproject.toml``.

Choose additional backends
----------------------------

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Route
     - Additional setup
   * - NEP
     - Matching model file; Calorine CPUNEP supports CPU inference.
   * - MACE
     - Install ``python -m pip install -e '.[mace]'`` and supply a checkpoint
       or configured foundation model. CPU is supported.
   * - VASP
     - Configure your licensed executable, MPI launcher, and POTCAR library.
       These are not distributed by this project.
   * - Thirdorder FC3
     - Make ``thirdorder_vasp.py`` available or set ``thirdorder-command``.
   * - FC4 / four-phonon transport
     - Install Fourthorder and FourPhonon separately; configure
       ``Fourthorder_vasp.py`` and the FourPhonon executable.
   * - Phonopy SSCHA
     - The installed Phonopy/symfc stack and an energy/force ASE calculator.
       The current route does not use ALAMODE.
   * - Slurm
     - A shared result directory, scheduler commands, and the same Python
       environment and backend executables on compute nodes.

Executable discovery can be checked with ``command -v phono3py`` and, when
needed, ``command -v thirdorder_vasp.py`` or ``command -v Fourthorder_vasp.py``.
A parser check alone does not test these programs or run a force calculation.

Next: :doc:`starting`. For testing and documentation builds, see
:doc:`development`.
