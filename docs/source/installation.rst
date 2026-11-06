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
``pyproject.toml``, including Calorine >= 3.5 for both CPU and GPU use.
``info`` checks configuration without starting a calculation.

Choose additional backends
------------------------------

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Route
     - Additional setup
   * - NEP
     - Matching model file; Calorine CPUNEP supports CPU inference.
   * - NEP GPU batches
     - Provide GPUMD's ``nep`` executable on a compatible GPU node.
       Python dependencies are included in the standard installation. See below.
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

Configure optional NEP GPU batches
-----------------------------------

On a local machine without a GPU, keep the normal installation above and the
default ``calculator.nep-force-backend: cpu``. CPUNEP does not need either the
``nep`` or ``gpumd`` executable. GPU batch inputs can still be prepared and
checked locally with ``nepkappa info``.

On the GPU server, activate the intended Python environment and run these
commands from the NEP-kappa checkout:

.. code-block:: bash

   python -m pip install -e .
   python -m pip check
   python -c "import calorine; from calorine.tools import batch_predict_properties; print(calorine.__version__)"

CPU and GPU use the same Python installation, which requires Calorine >= 3.5.
Source installations of Calorine also require a C++17 compiler. Installing
the Python package does not install GPUMD or a GPU driver/runtime.

Separately, load your cluster's GPUMD module or install the ``nep`` executable
following the `GPUMD installation instructions <https://gpumd.org/installation.html>`_.
Use a build compatible with the compute node's GPU and driver/runtime.
The GPUMD project provides both ``nep`` and ``gpumd``; this batch-prediction
route invokes only ``nep``, in prediction mode, and does not run molecular
dynamics with ``gpumd``.

.. code-block:: bash

   command -v nep

For NVIDIA nodes, run ``nvidia-smi`` inside the allocated GPU job to check
device visibility. A working executable path or Python import alone does
not verify GPU availability. In the input, select:

.. code-block:: yaml

   calculator:
     name: nep
     model: /path/to/nep.txt
     nep-force-backend: gpu-batch
     nep-batch-size: 128
     nep-command: /path/to/nep

``nep.txt`` is the trained model; ``nep-command`` is the executable path.
It can be omitted when ``nep`` is on PATH. These paths must be visible on
compute nodes. Request a GPU through the scheduler before running
``nepkappa fc2`` or ``nepkappa fc2fc3``; input validation does not allocate
one. See :doc:`input_files` for Slurm force-array settings and supported
workflow combinations. Missing dependencies or failed GPU predictions raise
an error instead of silently switching to CPU.

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
The same installation command checks the minimum Calorine version for both
CPU and GPU use.
Continue with :doc:`starting`; development and documentation build instructions
are in :doc:`development`.
