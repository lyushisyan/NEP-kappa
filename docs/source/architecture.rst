Calculation stages and routes
=============================

The two calculation routes are **phonon properties and thermal transport**
and **nonequilibrium phonon dynamics (TD-BTE)**. These names describe the
requested results: QHA and SSCHA in the first route are temperature dependent,
while the current TD-BTE route follows energy redistribution in time without
spatial drift. The routes share force-constant
preparation but neither route must run before the other.

Both from-scratch routes use six sections. Phonon properties and transport use
``structure``/``calculator``/``force-constant``/``kappa``/``plot``/``output``;
TD-BTE substitutes ``tdbte`` for ``kappa``. QHA, SSCHA, and four-phonon
options add steps only in the transport route. A compact ``tdbte``/``output``
input can instead reuse completed, mutually matching phono3py artifacts.
``nepkappa run`` derives the ordered steps from the selected input shape.

.. code-block:: text

   Preparation: structure + calculator -> displacements and forces
                                      -> FC2 (+ FC3 when needed) and metadata
                                                        |
                    +-----------------------------------+--------------------+
                    |                                                        |
       Phonon properties and thermal transport                  Nonequilibrium TD-BTE
       (six-section input)                                      (six sections or reuse)
       |                                                        |
       +-> FC2 -> dispersion, DOS, heat capacity                +-> 3ph kernel
       +-> FC2 + FC3 -> 3ph RTA/LBTE/Wigner kappa(T)             +-> n(q, branch, t)
       +-> optional QHA, SSCHA, FC4 -> selected results

   ``nepkappa run`` on a six-section TD-BTE input prepares matching FC2/FC3
   and metadata first; ``nepkappa fc2fc3`` can also prepare them separately.
   A calculated kappa(T) is not a TD-BTE input or prerequisite.

This is a data-dependency diagram, not a claim that all options use one fixed
order. QHA scans volumes and may regenerate force constants at each equilibrium
volume; SSCHA refits an auxiliary FC2; four-phonon transport additionally
requires FC4. The present TD-BTE builder consumes matching
``phono3py_disp.yaml``, ``fc2.hdf5``, and ``fc3.hdf5`` (or a previously built
compatible kernel). It does not read a ``kappa-m*.hdf5`` conductivity result.
The TD-BTE solver currently includes three-phonon collisions; QHA, SSCHA,
four-phonon and Wigner options are not automatically transferred to it.

Implementation behind the stages
--------------------------------

``command_registry.py`` defines the public command catalog. ``cli.py`` owns
terminal handling and concise error reporting. ``config.py`` validates the two
public YAML shapes and exposes typed views through ``config_models.py``.

``config.py`` also expands nested phonon/transport switches and infers the
full or artifact-reuse TD-BTE route.
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
``tdbte_accumulate.py`` optionally compiles block reductions. The TD-BTE stage
itself does not load a force calculator; the six-section ``run`` plan does so in
the preceding FC2/FC3 stage. Physical-rate and convergence validation remain separate
from its numerical conservation checks.

``scheduler.py`` and ``slurm.py`` manage batch submission. ``run_state.py``
records deferred execution; Slurm tools can inspect live queue state.
``artifacts.py`` stores force-job outputs and cache identities, and
``provenance.py`` records input hashes, versions, timings, and outputs.

Maintaining the input contract
------------------------------

Keep the six-section phonon/transport schema, six-section TD-BTE schema, and
two-section TD-BTE artifact-reuse schema explicit.
Introduce a new stage behind the application interface and add parser, stage,
and numerical checks appropriate to its behavior.
Keep backend-specific details out of CLI handlers. Update the public example
catalog, input reference, and skill when the user-facing interface changes.

The current tests cover unit behavior, orchestration, small numerical
regressions, and frozen scientific snapshots. Full production reproducibility
and convergence need separate calculations. See :doc:`development` and
:doc:`release_baseline`.
