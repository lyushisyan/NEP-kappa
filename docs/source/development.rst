Development
====================================

Repository contents
---------------------

Keep source code, documentation, public examples, and the input-assistant skill
in version control. ``tests/`` and ``benchmarks/`` are maintained locally and
ignored by Git. They are not included in a fresh GitHub checkout.

``calculations/`` is the local research workspace. Its guide is public, while
material runs, site-specific submission scripts, plots, force jobs, and archived
inputs are ignored, along with local manuscript and release-worktree notes.

Run checks
------------

From the repository root, build the public documentation:

.. code-block:: bash

   python -m pip install -e '.[docs]'
   python -m sphinx -b html -W --keep-going docs/source docs/build/html

``.github/workflows/checks.yml`` runs the strict documentation build on GitHub
pushes and pull requests using Python 3.11. It does not depend on local tests
or reference data.

If the local ``tests/`` and ``benchmarks/`` directories are available, also run:

.. code-block:: bash

   python -m pip install -e '.[dev]'
   python -m pytest -q

Tests use temporary result directories, calculator/scheduler mocks,
small numerical checks, and frozen data. They do not need the local
research archive, VASP, SSH, or a Slurm allocation. This suite does not replace
production convergence studies or running every external solver.

Source layout
---------------

- ``cli.py`` and ``command_registry.py``: public commands.
- ``config.py`` and ``config_models.py``: YAML parsing and typed configuration.
- ``application.py`` and ``stages/``: workflow orchestration and execution.
- ``adapters/``, ``calculators.py``, ``transport.py``, ``fourphonon.py``:
  calculator and solver interfaces.
- ``qha.py``, ``sscha.py``, ``qha_sscha.py``: temperature-dependent workflows.
- ``bubble.py``, ``approximations.py``: on-shell corrections and method scope.
- ``tdbte.py``: experimental artifact-based population dynamics and audits.
- ``artifacts.py``, ``provenance.py``, ``run_state.py``: cache and result identity.
- ``scheduler.py`` and ``slurm.py``: scheduler interaction.
- ``plot.py``, ``report.py``, ``convergence.py``: analysis.

``workflow.py`` retains compatibility methods used by the stage implementations.

Maintaining examples and the skill
------------------------------------

Prefer one example per distinct workflow and explain small variants in
``examples/README.md``. Do not duplicate a full YAML just to change RTA to LBTE.
Use public input paths and placeholder executable paths; specify
``submit: false`` in Slurm examples.

When keys or templates change, update the input reference and
``.agents/skills/nepkappa-input/`` together. The example-catalog test validates
workflow schemas, structures and model paths, and separate analysis schemas.
Keep compatibility-only inputs under ``tests/fixtures/``.

Baseline snapshots
--------------------

The local ``benchmarks/data/`` includes compact data and original BAs source summaries
needed to check their hashes without private calculation directories.
These checks verify hashes and recorded values without rerunning DFT.
See :doc:`release_baseline` for their scope.
