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
- ``tdbte.py``: population dynamics and audits; ``tdbte_builder.py`` builds native
  FC2/FC3 energy-shell chunks, ``tdbte_storage.py`` validates and streams them,
  and ``tdbte_quadrature.py`` owns the surface integration rule.
- ``artifacts.py``, ``provenance.py``, ``run_state.py``: cache and result identity.
- ``scheduler.py`` and ``slurm.py``: scheduler interaction.
- ``plot.py`` and ``report.py``: result analysis and reporting.

``workflow.py`` retains the shared context used by the stage implementations.

Maintaining examples and the skill
------------------------------------

The public catalog intentionally contains nine 3C-SiC static
workflows and one TD-BTE input. Explain their method differences in
``examples/README.md``. Use public input paths and placeholder executable
paths. Site-specific Slurm settings belong in an external batch script or
optional ``parallel`` settings.

When keys or templates change, update the input reference and
``.agents/skills/nepkappa-input/`` together. The example-catalog test validates
the nine static inputs, the separate TD-BTE schema, and their public assets.
Keep internal fixture inputs in ``tests/fixtures/`` aligned with the public
six-section static or two-section dynamic shapes.

Baseline snapshots
--------------------

The local ``benchmarks/data/`` includes compact data and original BAs source summaries
needed to check their hashes without private calculation directories.
These checks verify hashes and recorded values without rerunning DFT.
See :doc:`release_baseline` for their scope.
