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

The example inputs are in ``examples/``. Before starting a calculation,
inspect the parsed settings:

.. code-block:: bash

   nepkappa info examples/1-bulk-nep-rta.yaml

Run a complete workflow
-----------------------

For the default bulk NEP finite-displacement RTA example:

.. code-block:: bash

   nepkappa run examples/1-bulk-nep-rta.yaml

Run stages separately
---------------------

The same workflow can be split into stages:

.. code-block:: bash

   nepkappa relax examples/1-bulk-nep-rta.yaml
   nepkappa fc examples/1-bulk-nep-rta.yaml
   nepkappa kappa examples/1-bulk-nep-rta.yaml
   nepkappa plot examples/1-bulk-nep-rta.yaml

If ``relaxation.enabled`` is ``false``, the ``relax`` stage copies the input
structure to ``POSCAR_relaxed``.

For VASP examples, edit the VASP executable, MPI command, and POTCAR path before
running on your own machine or cluster.

Compare completed DFT and NEP results
-------------------------------------

After both result directories contain ``phono3py_disp.yaml``, ``fc2.hdf5``,
and ``kappa-m*.hdf5``, compare them with:

.. code-block:: bash

   nepkappa compare examples/compare.yaml

Comparison figures are written under ``compare.compare_dir/plots`` from the
comparison YAML file.

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
- ``plots/`` with dispersion, DOS, volume heat capacity, group velocity in km/s, relaxation time, and kappa figures
- ``comparison/plots/`` for DFT-vs-NEP comparison figures
