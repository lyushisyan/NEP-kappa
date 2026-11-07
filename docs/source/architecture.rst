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
options add steps only in the phonon-properties/transport route. A compact ``tdbte``/``output``
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

The dependency graph for a public calculation command is:

.. code-block:: text

   cli.py: command, terminal/log streams, provenance
      |
      +-> config.py: read YAML once, adapt values to argparse
      |      +-> config_schema.py: accepted keys and scalar conversion
      |      +-> config_plan.py: expand scientific switches, select steps
      |      +-> config_validation.py: validate the selected operation
      |      +-> config_models.py: flat WorkflowConfig and typed snapshots
      |
      +-> config_display.py: show active settings
      |
      +-> application.py: order stages, time them, defer/resume batch work
                +-> stages/structure.py and stages/force_constants.py
                +-> transport.py / fourphonon.py / tdbte.py
                +-> qha.py / sscha.py / qha_transport.py / qha_sscha.py
                +-> plotting / reports
                          |
                     workflow.py: shared context (no plan runner)
                          +-> calculators.py / adapters/vasp.py
                          +-> artifacts.py / slurm.py

   force_worker.py -> config.py -> workflow.py -> cached force calculation

``command_registry.py`` defines public command names and input-selected
operations. ``application.py`` is the scientific orchestration entry point.
It selects finite displacement, HiPhive or Thirdorder directly, without a
second callback-based FC2/FC3 dispatcher. Coupled temperature workflows still
own their physical volume/temperature loops and call the application stage
runner for each case. These loops are scientific algorithms, not redundant
command dispatchers.

``CommandExecution.status`` is the application/CLI outcome contract for
successful commands: ``complete`` means the requested calculations finished,
while ``prepared`` and ``submitted`` mean scientific work is still pending.
``run_stage`` retains the status returned by force, RTA/LBTE or FourPhonon scheduler
manifests; synchronous SSCHA transport shares the same execution record. A plan
stops when that record indicates deferred work. The CLI preserves exit code
zero for successful preparation/submission, but reports Prepared/Submitted and
writes the same status to provenance instead of claiming completion.
Exceptions still produce a failed command and nonzero exit code. Scheduler
submission manifests separately track the subsequent job lifecycle; finishing
a batch job does not retroactively rewrite the original command's provenance.

SSCHA and QHA-SSCHA reject configured nested Slurm RTA/LBTE in input validation
and again at their workflow entry points, before temperature calculations or
submission. An unexpected deferred result from nested transport raises an error
before later temperatures or workflow completion summaries. Both temperature
loops remain sequential; direct Slurm RTA/LBTE retains its deferred status behavior.
Pure-QHA transport also rejects nested Slurm transport before its temperature
loop and rejects unexpected prepared/submitted results from a synchronous call.

``workflow.py`` retains the shared structure, calculator instances, artifact
paths, cache identities, force-job preparation and export helpers. This context
is necessary for the stages and Slurm workers. Its former ``run*``, stage
forwarders and duplicate timing methods have been removed. ``stages/base.py``
keeps only prerequisite validation and execution; the unused prepare, collect
and status hooks have been removed. VASP operations use the existing adapter
directly rather than going through one-line context proxies.

Input parsing and configuration ownership
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The input path is YAML -> expanded settings -> argparse tokens ->
``WorkflowConfig``. Token conversion remains a deliberate compatibility
boundary: argparse still owns existing scalar coercions, boolean spellings,
list cardinalities, choices and defaults. There is no second type-conversion
framework. The redundant token-order table has been removed. ``fc2`` and
``kappa`` now read YAML and invoke argparse once, preserving the preflight
checks before validating their input-selected operation.

Strict YAML loading (including duplicate-key rejection) lives in
``config_schema.py`` and is shared by ``config.py`` and ``report.py``.
Configuration validation can still use the report result-directory resolver;
the report module no longer imports the parser. QHA equilibrium-volume loading,
temperature-range checks and primitive-volume scaling live in ``qha.py``;
both ``qha_transport.py`` and ``qha_sscha.py`` import them directly.

Planning, validation and display have separate modules. The flat configuration
is the canonical mutable state. ``WorkflowSections`` is an on-demand typed
snapshot, not a second parsing/validation pass. It is deliberately uncached:
QHA and SSCHA update copied configurations for individual volumes and
temperatures. Application code reads simple plan/output fields directly
instead of constructing every typed section for a single field.

Timing text is shared through ``output.py``. CLI total duration, stage timing
and activity feedback keep their existing output formats.

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

Keep six-section inputs mandatory for ``run``. Direct commands require only
the sections used by their stage; preserve command-specific validation and
compatibility with complete inputs. TD-BTE artifact reuse uses ``tdbte``, not
a compact ``run`` input.
Introduce a new stage behind the application interface and add parser, stage,
and numerical checks appropriate to its behavior.
Keep backend-specific details out of CLI handlers. Update the public example
catalog, input reference, and skill when the user-facing interface changes.

The public ``tests/public/`` suite covers parsing, stage selection, deferred
execution, force-worker context use, cache identity and small numerical
invariants without private datasets, a VASP executable or a Slurm cluster.
The separate local suite covers additional numerical regressions and frozen
scientific snapshots. Full production reproducibility
and convergence need separate calculations. See :doc:`development` and
:doc:`release_baseline`.

Refactoring and cache identity
-----------------------------------

Force-job fingerprint payloads retain their existing schema, inputs and hashing
rules. One input is the SHA-256 of ``workflow.py`` itself. Editing that file
therefore invalidates old force caches even when the numerical algorithm is
unchanged; a refactor must not substitute a historic hash to force cache reuse.
Provenance continues to record actual input, source/software and artifact
identities. This refactor does not change public YAML or scientific output
formats, units, mode ordering or numerical rejection/conservation policies.

Force-cache publication invalidates ``job-input.sha256`` before replacing
``forces.npy``. The JSON input record is installed next, followed by an atomic
replacement of the digest marker as the final commit point. An interrupted
update before that point has no valid marker and must recalculate. Invalid
new force arrays are rejected before the old force file and metadata are
changed. Thirdorder/Fourthorder jobs likewise invalidate their marker before
a backend can overwrite ``vasprun.xml``. This protects single-writer jobs from
process interruption; it does not introduce concurrent-writer locking or a
power-loss durability protocol. Identity schemas and final file formats are
unchanged, and no old identities are migrated or retroactively certified.

Force-cache audit and proposed replacement (not implemented)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The present schema-2 key is built by ``workflow._job_input_payload`` and
checked by ``ForceArtifactStore`` for finite-displacement/HiPhive forces;
Thirdorder/Fourthorder stages use the same fingerprint for their XML jobs.
It includes displacement POSCAR bytes, backend, calculator inputs, the package
version list and the entire ``workflow.py`` hash. NEP includes model bytes;
VASP includes command identity, merged default/user settings and POTCAR sources;
ASE/MACE includes factory source, kwargs and discovered/declared model files.

The dependency audit finds two opposing problems:

* Unrelated edits to ``workflow.py`` invalidate every force job. Versions of
  plotting, fitting and transport packages also affect force keys, although
  those operations do not calculate the force at a fixed geometry.
* Changes to ``adapters/vasp.py``, the force-handling code in ``calculators.py``,
  or a same-version rebuild of Calorine/native calculator libraries need not
  invalidate the key. Hashing a factory's source file does not cover imported
  helpers or native libraries. ``_command_cache_identity`` scans command tokens
  backwards for a file/executable; a file argument can be mistaken for the
  executable, and launcher/script/engine dependencies are not all covered.

In addition, VASP settings in the key precede INCAR/KPOINTS serialization, so
changes to those writers or to force parsing are missed. Structure matching
checks species/order, cell and positions; POSCAR does not represent all ASE
state (for example charges, magnetic moments or arbitrary calculator metadata).
The digest marker is compared directly, rather than authenticating the JSON
record or force-file contents. Finite force corruption is therefore outside
the current shape/finiteness checks. This is an input cache, not a checksum
audit of the force output.

A future fingerprint revision should use an explicit, small dependency list
per supported force backend, without automatic import-graph discovery:

1. **Actual force input:** ordered species, cell, positions and periodicity,
   plus supported force-affecting atom state, with exact deterministic
   serialization and no tolerance-based rounding. Include effective calculator
   parameters and content hashes of the model/checkpoint/parameter/POTCAR files
   actually loaded. For VASP, identify the generated scientific INCAR/KPOINTS
   settings and ordered POTCAR contents. Keep input paths and original files
   in provenance, separate from content identity.
2. **Force implementation:** explicitly include the workflow force-evaluation
   and order-job/XML-writing functions, backend construction/evaluation,
   VASP input writers and output parser, and any unit/order conversions they
   call. For in-repository Python, a normalized AST digest excluding comments,
   docstrings and location fields can avoid cosmetic invalidation; include
   referenced constants and helper implementations, not just function bodies.
   Backend modules can be hashed conservatively as whole normalized modules.
   Dependency-list changes must be reviewed with numerical code changes.
   This is a conservative source identity, not proof of semantic equivalence.
3. **Numerical engine:** include relevant provider versions and content hashes
   of provider implementation modules and native binaries (not just the Python
   factory). For NEP this includes CPUNEP and its compiled implementation;
   for ASE/MACE it includes the selected provider and numerical runtime;
   for VASP it includes the resolved engine binary and applicable wrapper
   scripts, with command tokens retained separately. Record precision/device
   and reproducibility-affecting execution settings. An opaque plugin or
   wrapper with unidentifiable dependencies cannot certify persistent reuse.

Displacement generation and FC fitting do not belong in a fixed-geometry
force identity: their changed geometry is already an input, and fitting
consumes forces later. Thirdorder/Fourthorder tool identities remain provenance
for generation/reaping; XML force encoding still belongs in the force contract.
Full source hashes, broad environment versions, and run provenance should
remain available independently of the narrower reuse key. A new record should
also checksum its force artifact and verify it before reuse.

Adoption requires a new fingerprint schema and explicit invalidation. Old or
incomplete identities must never be translated into new valid markers. Keep
old results intact and recompute in a fresh result directory; only successful,
validated new force calculations may acquire the new identity. No such change
or migration is made here. The three dead private-method deletions still
change the current ``workflow.py`` hash and therefore invalidate existing
schema-2 force caches on their next use.

The removed Python stage-forwarding methods were internal entry points. Local
tests now call the application runner, scientific stage or VASP adapter that
owns the operation. Archived calculation packages contain their own snapshots;
those snapshots are not rewritten as part of a source refactor.
