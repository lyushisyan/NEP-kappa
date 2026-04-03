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

Install the required Python packages:

.. code-block:: bash

   pip install numpy h5py ase matplotlib seekpath
   pip install phonopy phono3py hiphive trainstation calorine

Clone the repository
--------------------

.. code-block:: bash

   git clone https://github.com/lyushisyan/NEP-kappa.git
   cd NEP-kappa

Verify installation
-------------------

Make sure ``phono3py`` is available:

.. code-block:: bash

   phono3py --version

If the command runs successfully, the environment is ready for NEP-kappa.