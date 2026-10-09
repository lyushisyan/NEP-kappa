Quick start
=============

Check the bundled 3C-SiC calculation
------------------------------------

After :doc:`installation`, run these commands from the repository root:

.. code-block:: bash

   nepkappa info examples/nep-rta-wigner-3ph.yaml

To run this example in an appropriate compute allocation, execute
``nepkappa run examples/nep-rta-wigner-3ph.yaml``. It relaxes 3C-SiC,
generates FC2/FC3, and computes three-phonon RTA and SMM19 Wigner transport.
Plots and a report are generated separately. Results appear in
``calculations/example-runs/3c-sic-nep-rta-wigner-3ph/``:

- ``run.log``: progress and errors.
- ``fc2.hdf5``, ``fc3.hdf5``, ``phono3py_disp.yaml``: force constants and metadata.
- ``kappa-m*.hdf5``: conductivity and mode data.
- ``plots/``: figures.
- ``report.md`` and ``report.yaml``: result summaries.

The template demonstrates the workflow. Converge supercell size, q mesh, and
force-constant settings before using a result quantitatively.

Create your own input
-----------------------

.. code-block:: bash

   cp examples/nep-rta-wigner-3ph.yaml input.yaml
   nepkappa info input.yaml
   nepkappa run input.yaml

Edit the copied file's structure, calculator, numerical settings, and output
directory for your material. Provide a structure and potential for the same
material. The six-section static format is described in :doc:`input_files`.

Ordinary YAML paths resolve from the **launch directory**, not the YAML file's
directory. To work elsewhere, use suitable relative paths or absolute paths.
Choose another workflow
-------------------------

Select a template from :doc:`examples` and run it with the same ``run`` command.
The QHA, SSCHA, and four-phonon switches are inside ``force-constant``.
For time-dependent BTE, start from ``examples/tdbte.yaml``; it uses only
``tdbte`` and ``output`` sections.

For analysis or restarts, use a stage command instead:

.. code-block:: bash

   nepkappa info input.yaml --for kappa
   nepkappa stage kappa input.yaml

This recalculates transport from existing matching FC2, FC3, and phono3py
metadata. See :doc:`tutorial` for stage prerequisites, QHA/SSCHA, and Slurm.

``info`` checks supported settings; it does not test model accuracy or
external executables. ``nepkappa status input.yaml`` reads saved submission
state and, where available, Slurm queue state. For a traceback, put ``--debug``
before the command.
