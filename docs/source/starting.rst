Quick start
=============

Run the bundled Si calculation
--------------------------------

After :doc:`installation`, run these commands from the repository root:

.. code-block:: bash

   nepkappa validate examples/nep-rta.yaml
   nepkappa run examples/nep-rta.yaml
   nepkappa plot examples/nep-rta.yaml
   nepkappa report examples/nep-rta.yaml

The calculation relaxes bulk Si, generates FC2/FC3, and computes three-phonon
RTA conductivity. Plots and a report are generated separately.
Results appear in ``calculations/example-runs/nep-rta/``:

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

   nepkappa init input.yaml
   nepkappa validate input.yaml
   nepkappa info input.yaml
   nepkappa run input.yaml

The initializer asks for the calculation type, structure, calculator, numerical
settings, and output directory. Provide a structure and potential for the same
material. Existing input files are preserved unless ``--force`` is supplied.

Ordinary YAML paths resolve from the **launch directory**, not the YAML file's
directory. To work elsewhere, use suitable relative paths or absolute paths.
Convergence study ``base`` and ``study.directory`` are exceptions: they resolve
from the study YAML's directory.

For scripts, create a non-interactive input:

.. code-block:: bash

   nepkappa init input.yaml --non-interactive \
     --preset three-phonon --calculator nep \
     --structure examples/structures/Si/POSCAR_bulk \
     --model potentials/Si/Si_Bulk_Fan.txt \
     --dim 3 3 3 --mesh 21 21 21

Choose another workflow
-------------------------

Select a template from :doc:`examples` and run it with the same ``run`` command.
The presets are ``three-phonon``, ``four-phonon``, ``qha``, ``scph``, and
``qha-sscha``. ``scph`` selects the implemented Phonopy stochastic SSCHA route.

For analysis or restarts, use a stage command instead:

.. code-block:: bash

   nepkappa validate input.yaml --for kappa
   nepkappa kappa input.yaml

This recalculates transport from existing matching FC2, FC3, and phono3py
metadata. See :doc:`tutorial` for stage prerequisites, QHA/SSCHA, and Slurm.

Validation checks supported settings; it does not test model accuracy or
external executables. ``nepkappa status input.yaml`` reads saved submission
state and, where available, Slurm queue state. For a traceback, put ``--debug``
before the command.
