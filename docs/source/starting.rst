Quick Start
===========

Install the package
-------------------

From the repository root:

.. code-block:: bash

   python -m pip install -e .

Check that the command is available:

.. code-block:: bash

   nepkappa --help

Inspect an example
------------------

The example inputs are in ``examples-input/``. Before starting a calculation,
inspect the parsed settings:

.. code-block:: bash

   nepkappa info examples-input/input_1-bulk-nep-rta.yaml

Run a complete workflow
-----------------------

For the default bulk NEP finite-displacement RTA example:

.. code-block:: bash

   nepkappa run examples-input/input_1-bulk-nep-rta.yaml

Run stages separately
---------------------

The same workflow can be split into stages:

.. code-block:: bash

   nepkappa relax examples-input/input_1-bulk-nep-rta.yaml
   nepkappa fc examples-input/input_1-bulk-nep-rta.yaml
   nepkappa kappa examples-input/input_1-bulk-nep-rta.yaml

If ``relaxation.enabled`` is ``false``, the ``relax`` stage just keeps the input
structure and records the preparation step.

Available examples
------------------

- ``examples-input/input_1-bulk-nep-rta.yaml``: bulk, NEP forces, finite displacement, RTA
- ``examples-input/input_2-bulk-nep-hiphive-rta.yaml``: bulk, NEP forces, HiPhive, RTA
- ``examples-input/input_3-bulk-nep-lbte.yaml``: bulk, NEP forces, finite displacement, LBTE
- ``examples-input/input_4-bulk-nep-rta-wigner.yaml``: bulk, NEP forces, finite displacement, Wigner transport
- ``examples-input/input_5-bulk-vasp-rta.yaml``: bulk, VASP relaxation and VASP forces, finite displacement, RTA
- ``examples-input/input_6-bulk-vasp-hiphive-rta.yaml``: bulk, VASP relaxation and VASP forces, HiPhive, RTA
- ``examples-input/input_7-film-nep-hiphive-rta.yaml``: film, NEP forces, HiPhive, RTA

For VASP examples, edit the VASP executable, MPI command, and POTCAR path before
running on your own machine or cluster.

Typical outputs
---------------

Generated files are written under ``output.result_dir`` from the YAML file.
Typical outputs include:

- ``run.log``
- ``POSCAR_relaxed`` when relaxation is enabled
- ``phono3py_disp.yaml``
- ``fc2.hdf5``
- ``fc3.hdf5``
- ``hiphive_model.fcp`` for HiPhive runs
- ``vasp-runs/`` and ``vasp-relax/`` for VASP runs
- ``kappa-mXXXXX.hdf5`` from ``phono3py``
