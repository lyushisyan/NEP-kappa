Architecture
============

The CLI reads a YAML configuration, resolves a workflow plan, and invokes
stage handlers through the application layer. Calculators supply energies,
forces, or stresses; the force-constant and transport stages consume those
results and record artifacts and provenance.

.. code-block:: text

   CLI / initializer
          |
   configuration + workflow plan
          |
   application / stage runner
          |
          +-- structure relaxation
          +-- FC2 / FC3 / FC4 generation
          +-- QHA / Phonopy SSCHA / QHA+SSCHA
          +-- phono3py / FourPhonon transport
          +-- plotting / reports

   Shared services: calculators, scheduler, artifacts, provenance, run state

Module boundaries
-----------------

``command_registry.py`` defines the public command catalog. ``cli.py`` owns
terminal handling and concise error reporting; ``initializer.py`` writes
validated starter inputs. ``config.py`` validates YAML and exposes typed views
through ``config_models.py`` while preserving historical attributes.

``application.py`` constructs stage implementations and executes plans.
``stages/`` contains common lifecycle hooks and force/structure stages.
``workflow.py`` provides shared context and compatibility methods; it is
still a larger module, as is the centralized parser.

``calculators.py`` supplies NEP, MACE, ASE, and plugin interfaces.
``adapters/vasp.py`` owns VASP-specific inputs, execution, and output parsing.
``transport.py`` and ``fourphonon.py`` isolate transport interfaces.
Temperature workflows are in ``qha.py``, ``sscha.py``, and ``qha_sscha.py``.

``scheduler.py`` and ``slurm.py`` manage batch submission. ``run_state.py``
records deferred execution; ``status`` can combine it with live queue state.
``artifacts.py`` stores force-job outputs and cache identities, and
``provenance.py`` records input hashes, versions, timings, and outputs.

Maintaining compatibility
-------------------------

Scientific behavior and existing accepted YAML keys should stay stable during
structural cleanup. Introduce a new stage behind the application interface and
add parser, stage, and numerical checks appropriate to its behavior.
Keep backend-specific details out of CLI handlers. Update the public example
catalog, input reference, and skill when the user-facing interface changes.

The current tests cover unit behavior, orchestration, small numerical
regressions, and frozen scientific snapshots. Full production reproducibility
and convergence need separate calculations. See :doc:`development` and
:doc:`release_baseline`.
