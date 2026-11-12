Calculation workflows and commands
==================================

NEP-kappa uses YAML to configure structures, atomic-force calculations,
force constants, lattice thermal transport and nonequilibrium phonon dynamics.

Use ``nepkappa run`` for a complete workflow, or direct commands to execute
selected stages. This page introduces the input organization, common settings
and the scope of each module. New users should first read :doc:`starting`.

Complete workflows
------------------------

A complete YAML input has six main sections:

.. list-table:: YAML configuration sections
   :header-rows: 1
   :widths: 30 70

   * - Section
     - Purpose
   * - ``structure``
     - Crystal structure and relaxation.
   * - ``calculator``
     - Interatomic model, energies and atomic forces.
   * - ``force-constant``
     - Generation of second-, third- and fourth-order force constants.
   * - ``kappa`` / ``tdbte``
     - Lattice thermal transport or nonequilibrium phonon dynamics.
   * - ``plot``
     - Plot contents, formats and output settings.
   * - ``output``
     - Result-directory settings.

``kappa`` and ``tdbte`` select separate routes and cannot be combined in the
same complete input. Sections without explicit options may be written as ``{}``.

Save the input as ``input.yaml`` and run:

.. code-block:: bash

   nepkappa info input.yaml
   nepkappa run input.yaml

``info`` checks the configuration. ``run`` organizes and executes stages
according to the YAML settings. For example, a three-phonon transport workflow
relaxes the structure, generates FC2/FC3 and calculates conductivity, saving
the results of each stage for later use.

Direct stage calculations
-------------------------

Use direct commands when only part of a calculation is needed.
For example, after generating FC2/FC3, ``nepkappa kappa`` can calculate
conductivity at different temperatures or q meshes without recalculating forces.

Such an input needs only relevant sections such as ``kappa`` and ``output``.
``kappa.force-constants`` selects the existing force-constant directory and
``output.result-dir`` selects the new result directory.

.. code-block:: bash

   nepkappa info kappa.yaml --for kappa
   nepkappa kappa kappa.yaml

The first command checks the direct-stage input; the second calculates conductivity.

Other direct commands are:

.. list-table:: Direct commands
   :header-rows: 1
   :widths: 28 72

   * - Command
     - Purpose and required sections
   * - ``relax``
     - Relax the structure; requires ``structure`` and ``calculator``.
   * - ``fc2``
     - Generate FC2; requires ``structure``, ``calculator`` and ``force-constant``.
   * - ``fc2fc3``
     - Generate FC2/FC3 using the same sections.
   * - ``fc2fc3fc4``
     - Generate FC2/FC3/FC4 using the same sections.
   * - ``kappa``
     - Calculate conductivity from existing force constants; requires ``kappa``.
   * - ``tdbte``
     - Solve the time-dependent BTE using existing force constants or a collision kernel; requires ``tdbte``.
   * - ``plot``
     - Plot saved results using ``plot`` and ``output``.
   * - ``report``
     - Generate a report; specify ``output.result-dir`` or pass a result directory directly.

Direct commands accept compact inputs or the relevant settings in a complete YAML input.

.. note::

   - Relative paths are resolved from the launch directory, rather than the YAML file's directory.
   - Use a separate ``output.result-dir`` for each calculation to keep results organized.
   - Check compact inputs with ``nepkappa info input.yaml --for <command>``.
   - If relaxation is enabled, run ``relax`` before a direct force-constant command.
