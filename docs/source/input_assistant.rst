AI Input Assistant
====================

The repository includes the ``nepkappa-input`` Agent Skill in
``.agents/skills/nepkappa-input/``. It helps an AI assistant create, edit,
explain, and validate NEP-kappa YAML inputs from a calculation description.
It covers complete workflows, individual stages, plotting, result comparison,
and convergence-study inputs, and responds in the user's language.

It starts with the descriptive templates listed in :doc:`examples`, including
QHA, Phonopy SSCHA, QHA+SSCHA, and combined three/four-phonon transport.
It checks model species and identity instead of trusting filenames alone.

The skill uses this checkout's parser, input documentation, and examples.
It does not require its own model API key or a separate agent service; the
user supplies a compatible AI assistant. Automatic validation also needs
terminal access and a matching NEP-kappa installation.

Using the skill
-----------------

In Codex, open this repository and ask, for example:

.. code-block:: text

   Use $nepkappa-input to prepare a three-phonon RTA input for bulk Si.
   Use examples/structures/Si/POSCAR_bulk and potentials/Si/Si_Bulk_Fan.txt, at 300 K.
   Write si-rta.yaml and validate it. Do not start the calculation.

Chinese requests work as well:

.. code-block:: text

   使用 nepkappa-input，帮我给 3C-SiC 写一个 NEP 三声子热导率输入文件。
   结构在 calculations/3C-SiC/POSCAR，势文件在 potentials/3C-SiC/nep_3C-SiC.txt。
   温度范围是 300–1000 K，间隔 100 K。先检查文件，缺少的信息请问我。

The SiC paths above illustrate user-supplied paths; provide your actual files.
The skill must not replace missing SiC inputs with bundled Si examples.
The general NEP89 model is stored separately at
``potentials/nep89_20250409.txt``. Covering the structure's elements does not
alone prove accuracy for its phase or temperature range.

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

These links share the same maintained skill. Do not replace an existing
installation blindly. Other hosts may have different installation and
invocation conventions; follow their Agent Skills documentation. The core
instructions use ordinary files and CLI commands; ``agents/openai.yaml`` is
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

   nepkappa validate input.yaml --for run
   nepkappa info input.yaml --for run

Stage-only inputs use the corresponding target, such as ``--for kappa``.
Comparison and convergence inputs use their separate Python parsers as
described in the skill; they are not supported ``validate --for`` targets.

Parser validation, file/environment checks, and numerical convergence are
reported separately. A passing parser does not establish that model files
exist, a cluster is configured correctly, or a q mesh is converged. Preparing
input files alone does not start calculations or submit Slurm jobs.

Maintaining the skill
-----------------------

Update the skill when input keys, workflow behavior, or examples change.
Keep exact defaults and schemas in the software rather than duplicating them
in the skill. Verify representative generated inputs with the matching
installed validator, including stage-only inputs and invalid configurations.
