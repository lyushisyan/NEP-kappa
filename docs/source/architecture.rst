Architecture
============

The user-facing calculation input has two shapes. Static inputs contain six
sections and enable QHA, SSCHA, or four-phonon work inside
``force-constant``. Dynamic TD-BTE inputs contain ``tdbte`` and ``output``
and reuse completed FC2/FC3. The CLI compiles either shape into an ordered
plan and invokes the scientific stage handlers. Calculators supply energies,
forces, or stresses for static stages; dynamic TD-BTE reads existing matching
force constants and metadata.

.. code-block:: text

   Static: structure + calculator + force-constant + kappa + plot + output
                  |
          structure -> forces -> FC2/FC3 -> QHA/SSCHA/FC4 -> 3ph/4ph
                                                          |
                                                     plots/report

   Dynamic: tdbte + output -> existing FC2/FC3 -> kernel -> populations(t)

   Both shapes -> application / ordered stage runner

   Shared services: calculators, scheduler, artifacts, provenance, run state

Implementation behind the stages
--------------------------------

``command_registry.py`` defines the public command catalog. ``cli.py`` owns
terminal handling and concise error reporting. ``config.py`` validates YAML and exposes typed views
through ``config_models.py`` while preserving historical attributes.

``config.py`` also expands static nested switches and infers the separate
dynamic TD-BTE route. ``stage_plan.py`` validates older explicit stage choices
and compiles them to existing step names and backend settings.
``application.py`` constructs stage implementations
and executes those steps in order.
``stages/`` contains common lifecycle hooks and force/structure stages.
``workflow.py`` provides shared context and compatibility methods; it is
still a larger module, as is the centralized parser.

``calculators.py`` supplies NEP, MACE, ASE, and plugin interfaces.
``adapters/vasp.py`` owns VASP-specific inputs, execution, and output parsing.
``transport.py`` and ``fourphonon.py`` isolate transport interfaces.
Temperature workflows are in ``qha.py``, ``sscha.py``, and ``qha_sscha.py``.
``bubble.py`` writes separate diagonal on-shell frequency corrections;
``approximations.py`` records interpretation limits. ``tdbte.py`` orchestrates
kernel preparation and audited population dynamics. ``tdbte_builder.py``
evaluates native vertices from existing FC2/FC3; ``tdbte_quadrature.py`` integrates
linear shells; ``tdbte_storage.py`` validates and streams checksummed chunks.
``tdbte_accumulate.py`` optionally compiles block reductions. This route does not
load a force calculator. Physical-rate and convergence validation remain separate
from its numerical conservation checks.

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
