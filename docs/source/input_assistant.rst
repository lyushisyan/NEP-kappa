Input assistant
====================

The ``nepkappa-input`` skill provides instructions for preparing and checking
YAML inputs with an AI assistant. It is stored in
``.agents/skills/nepkappa-input/`` and uses the parser, examples, and documentation
from this checkout. Local validation requires terminal access and an installed
NEP-kappa environment.
For static calculations, it prefers the six-section input with QHA, SSCHA,
and four-phonon switches inside ``force-constant``. Dynamic TD-BTE uses
``tdbte`` and ``output``. Existing stage-first and preset inputs remain
compatible.

Using the skill
-----------------

In Codex, open this repository and ask, for example:

.. code-block:: text

   Use $nepkappa-input to prepare a three-phonon RTA input for bulk Si.
   Use examples/structures/Si/POSCAR_bulk and potentials/Si/Si_Bulk_Fan.txt, at 300 K.
   Write si-rta.yaml and validate it. Do not start the calculation.

For 3C-SiC:

.. code-block:: text

   使用 nepkappa-input，帮我给 3C-SiC 写一个 NEP 三声子热导率输入文件。
   结构在 calculations/3C-SiC/POSCAR，势文件在 potentials/3C-SiC/nep_3C-SiC.txt。
   温度范围是 300–1000 K，间隔 100 K。先检查文件，缺少的信息请问我。

Replace the SiC paths with your own files. Specify whether the request is to
prepare an input, generate plots, or run a calculation.

Plotting from existing results:

.. code-block:: text

   使用 $nepkappa-input，检查 calculations/3C-SiC/results 里现有的数据。
   我只想画色散、DOS、体积热容和群速度，单图和组合图都要。
   如果只有 FC2，请不要要求 FC3 或重新计算力常数。
   先生成并验证绘图输入，告诉我缺什么；暂时不要绘图。

VASP setup:

.. code-block:: text

   使用 $nepkappa-input，帮我检查目标服务器上的 VASP 和 POTCAR 在哪里，
   检查现有 input.yaml 的启动命令、MPI 环境与 Slurm 配置。
   如果需要改路径，请修改输入文件，不要改源码，也不要提交计算。

Remote file checks require access to that server. For TD-BTE inputs, also supply
the existing kernel NPZ and JSON metadata described in :doc:`tdbte`.

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

   nepkappa validate input.yaml --for run
   nepkappa info input.yaml --for run

Stage-only inputs use the corresponding target, such as ``--for kappa``.
Comparison and convergence inputs use their separate Python parsers as
described in the skill; they are not supported ``validate --for`` targets.

Validation checks the configuration, not model accuracy or numerical convergence.
File and environment checks are reported separately. Preparing an input does
not submit a calculation.

Maintaining the skill
-----------------------

Update the skill alongside parser and example changes. Validate generated
inputs against the same checkout.
