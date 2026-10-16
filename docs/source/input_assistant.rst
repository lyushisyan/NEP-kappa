Input assistant
====================

The ``nepkappa-input`` skill provides instructions for preparing and checking
YAML inputs with an AI assistant. It is stored in
``.agents/skills/nepkappa-input/`` and uses the parser, examples, and documentation
from this checkout. Local validation requires terminal access and an installed
NEP-kappa environment.
For phonon properties and thermal transport, it prefers the six-section
input with QHA, SSCHA, and four-phonon switches inside ``force-constant``.
Nonequilibrium phonon dynamics uses six sections with ``tdbte`` in place of
``kappa`` when force constants must be generated. A compact ``tdbte``/``output``
input reuses existing artifacts. Older stage-first and preset inputs are rejected.

Using the skill
-----------------

In Codex, open this repository and ask, for example:

.. code-block:: text

   Use $nepkappa-input to prepare a three-phonon RTA input for bulk Si.
   Use examples/structures/Si/POSCAR_bulk and potentials/Si/Si_Bulk_Fan.txt, at 300 K.
   Write si-rta.yaml and validate it. Do not start the calculation.

For 3C-SiC:

.. code-block:: text

   Use $nepkappa-input to prepare a three-phonon NEP input for 3C-SiC.
   Use calculations/3C-SiC/POSCAR and potentials/3C-SiC/nep_3C-SiC.txt.
   Calculate 300–1000 K in 100 K steps. Check the files and validate the input.
   Do not start the calculation.

Replace the SiC paths with your own files. Specify whether the request is to
prepare an input, generate plots, or run a calculation.

Plotting from existing results:

.. code-block:: text

   Use $nepkappa-input to inspect calculations/3C-SiC/results.
   Prepare an input for dispersion, DOS, volumetric heat capacity and group
   velocity, with separate and combined figures. Reuse the existing FC2
   and metadata. Validate the input without generating plots yet.

VASP setup:

.. code-block:: text

   Use $nepkappa-input to check VASP and POTCAR locations on my server,
   then inspect the launch command, MPI environment and Slurm settings
   in input.yaml. Correct input paths if needed. Do not submit a calculation.

Remote file checks require access to that server. For TD-BTE inputs, supply
either matching ``phono3py_disp.yaml``, ``fc2.hdf5``, and ``fc3.hdf5`` in one
directory, or a completed reusable kernel described in :doc:`tdbte`.

Availability in assistants
--------------------------

Codex discovers the repository skill under ``.agents/skills/``. If needed,
start a new session in the checkout. In another assistant that supports the
Agent Skills format, install or link the entire ``nepkappa-input`` directory
in that assistant's skill location. For example, from the repository root:

.. code-block:: bash

   # Claude Code project installation; only needed when not already installed.
   mkdir -p .claude/skills
   ln -s ../../.agents/skills/nepkappa-input .claude/skills/nepkappa-input

   # Cursor project installation; only needed when not already installed.
   mkdir -p .cursor/skills
   ln -s ../../.agents/skills/nepkappa-input .cursor/skills/nepkappa-input

Create these links only if the destination does not already exist. Other
assistants may use different skill directories. ``agents/openai.yaml`` contains
optional OpenAI-specific display metadata.

A copied skill still needs access to the matching checkout's documentation
and examples. Tell the assistant where NEP-kappa is located if it is not
working inside this repository. An assistant without filesystem access can
provide YAML text, but cannot verify local files or run the local validator.

Expected output
-----------------

The assistant delivers the YAML file, launch directory, intended command,
assumptions, and validation status. For a complete workflow it uses:

.. code-block:: bash

   nepkappa info input.yaml --for run

Stage-only inputs use the corresponding target, such as ``--for kappa``.
``info`` checks the configuration, not model accuracy or numerical convergence.
File and environment checks are reported separately. Preparing an input does
not submit a calculation.

Maintaining the skill
-----------------------

Update the skill alongside parser and example changes. Validate generated
inputs against the same checkout.
