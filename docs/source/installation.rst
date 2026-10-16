Installation and environment setup
======================================

System requirements
-----------------------

NEP-kappa requires Python 3.10 or later, matching the minimum required by its
Phonopy and phono3py dependencies. Python 3.11 is used in this example.
Linux is suitable for local and cluster calculations; external executables and
MPI must match the target machine. See :doc:`requirements` for dependencies.

Create a Python environment
-------------------------------

.. code-block:: bash

   conda create -n nepkappa python=3.11
   conda activate nepkappa

Alternatively, create and activate a virtual environment with
``python -m venv .venv`` and ``source .venv/bin/activate``.

Obtain the source
---------------------

.. code-block:: bash

   git clone https://github.com/lyushisyan/NEP-kappa.git
   cd NEP-kappa

Install and check
---------------------

.. code-block:: bash

   python -m pip install -e .
   python -m pip check
   nepkappa --version
   nepkappa info examples/nep-rta-wigner-3ph.yaml

The editable installation uses this checkout directly. Python dependencies
for the standard NEP/Phonopy/phono3py workflow are installed from
``pyproject.toml``. ``info`` checks configuration without starting a calculation.

Choose additional backends
------------------------------

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

Environment management and updates
--------------------------------------

Activate the environment in each new shell before using NEP-kappa. For a Git
checkout, review local changes before updating:

.. code-block:: bash

   conda activate nepkappa
   cd /path/to/NEP-kappa
   git status
   git pull --ff-only
   python -m pip install -e .
   python -m pip check
   nepkappa --version

Use the same environment and external programs on cluster compute nodes.
Continue with :doc:`starting`; development and documentation build instructions
are in :doc:`development`.
