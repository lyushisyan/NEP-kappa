Installation
============

It is recommended to install **NEP-kappa** in a dedicated Conda environment.

Create a Conda environment
--------------------------

.. code-block:: bash

   conda create -n nepkappa-env python=3.10 -y
   conda activate nepkappa-env

Install dependencies
--------------------

Install the package from the repository root. This also installs the Python
dependencies declared by NEP-kappa.

.. code-block:: bash

   git clone https://github.com/lyushisyan/NEP-kappa.git
   cd NEP-kappa
   python -m pip install -e .

If an older ``phono3py`` is already installed in the environment, upgrade the
editable install:

.. code-block:: bash

   python -m pip install --upgrade -e .

Verify installation
-------------------

Make sure the command-line entry point and ``phono3py`` are available:

.. code-block:: bash

   nepkappa --help
   python -c "import phono3py; print(phono3py.__version__)"

If the command runs successfully, the environment is ready for NEP-kappa.

Development checks
------------------

For development, install the optional test dependencies and run the tests:

.. code-block:: bash

   python -m pip install -e ".[dev]"
   python -m pytest
