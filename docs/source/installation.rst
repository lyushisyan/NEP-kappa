Installation and environment setup
==================================

NEP-kappa supports Python 3.10 or later; Python 3.11 is recommended.
Python dependencies, including MACE, are managed in ``pyproject.toml``.
Configure external solvers according to the calculation you plan to run.

Quick installation
------------------------

Create an isolated Python environment with Conda:

.. code-block:: bash

   conda create -n nepkappa python=3.11
   conda activate nepkappa

Download the source and install:

.. code-block:: bash

   git clone https://github.com/lyushisyan/NEP-kappa.git
   cd NEP-kappa
   python -m pip install -e .

This installs NEP-kappa in editable mode together with the MACE backend,
which is convenient for updating or developing the source.

Check the installation:

.. code-block:: bash

   python -m pip check
   nepkappa --version
   nepkappa info examples/nep-rta-wigner-3ph.yaml

If these commands succeed, the package is installed and can parse the example
input. ``nepkappa info`` checks the configuration without starting a calculation.

.. note::

   You can also create an environment with Python's built-in ``venv``:

   .. code-block:: bash

      python -m venv .venv
      source .venv/bin/activate

Calculation backends and additional dependencies
------------------------------------------------

The installation above includes NEP, MACE, Phonopy, phono3py and HiPhive
dependencies. Configure the following external programs as needed.

.. list-table::
   :header-rows: 1
   :widths: 25 75
   :width: 100%

   * - Backend
     - Additional requirements
   * - **VASP**
     - Install a licensed VASP executable and configure MPI, POTCAR and the runtime environment.
   * - **Thirdorder / Fourthorder**
     - Configure ``thirdorder_vasp.py`` and ``Fourthorder_vasp.py``.
   * - **FourPhonon**
     - Install FourPhonon separately; four-phonon Wigner calculations also require Wigner_Park.

.. note::

   - External programs must be built for the operating system, compiler and MPI environment in use.
   - Install the programs required by your calculation route.

Updating the software
------------------------

Enter your existing NEP-kappa checkout, obtain the latest source and update:

.. code-block:: bash

   conda activate nepkappa
   cd /path/to/NEP-kappa
   git status
   git pull --ff-only
   python -m pip install -e .

Check the updated environment:

.. code-block:: bash

   python -m pip check
   nepkappa --version

Review local changes before updating to avoid conflicts with remote changes.

Continue with :doc:`starting` to run your first example.
