Development
====================================

Repository contents
---------------------

Keep source code, documentation, public examples, and the input-assistant skill
in version control. ``tests/public/`` is a small public test suite; two opt-in
GPU validation and benchmark scripts are published in ``tests/gpu/``. Other
``tests/`` contents and ``benchmarks/`` remain local and ignored by Git.

``calculations/`` is the local research workspace. Its guide is public, while
material runs, site-specific submission scripts, plots, force jobs, and archived
inputs are ignored, along with local manuscript and release-worktree notes.
Generated artifacts in ``output/`` and ``outputs/``, local manual copies and
archives, and the scratch ``texput.log`` are also ignored.

Run checks
------------

From the repository root, run the public tests and build the documentation:

.. code-block:: bash

   python -m pip install -e '.[dev,docs]'
   RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest tests/public -q
   python -m sphinx -b html -W --keep-going docs/source docs/build/html

``.github/workflows/checks.yml`` runs the public tests and strict documentation
build on GitHub pushes and pull requests using Python 3.11. The public tests
use checked-in examples and synthetic temporary data. They require no private
models/datasets, VASP executable, network service or Slurm cluster. Large
scientific benchmarks remain separate from this CI job.

If the local ``tests/`` and ``benchmarks/`` directories are available, also run:

.. code-block:: bash

   python -m pip install -e '.[dev]'
   RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q

Keep the collection scope consistent when comparing regression totals.
``pytest tests/public`` runs the public suite; ``pytest tests`` also runs local
tests in that directory, but excludes research tests under ``calculations/``.
The root-level command above includes those research tests when present.
Use ``python -m pytest --collect-only -q`` to compare complete node IDs,
including parametrized cases, before interpreting a changed test count.

Most tests use temporary result directories, calculator/scheduler mocks,
small numerical checks, and frozen data. Some TD-BTE research checks also
import scripts from ``calculations/3C-SiC/tdbte-pilot/``; the full local suite
requires those scripts in addition to ``tests/`` and ``benchmarks/``.
These checks do not launch VASP, SSH, or Slurm jobs. They do not replace
production convergence studies or running every external solver.

The focused cache/status cleanup was validated on 2026-11-04 with 779 passing
repository-wide cases (65 public), compared with 760 (46 public) before it.
The increase of 19 comprises two additional direct-FC2/full-run Slurm
combinations, five cache interruption/recovery cases, nine execution-status
cases (including nested SSCHA transport), two strict-YAML cases and one shared
QHA-helper case. The original two Slurm parameter IDs were expanded into four;
all other prior node IDs remain, including the 23 research cases outside
``tests/``. No cases were skipped or deselected.

The subsequent synchronous-SSCHA correction passed 794 cases (80 public),
with no failures or skips. Three tests that mocked the entire SSCHA ``run()``
method were replaced by 18 cases: four CLI rejections, four workflow-entry
rejections, eight real-temperature-loop checks and two direct Slurm LBTE
acceptance checks. All other prior node IDs remain. Nested Slurm LBTE with
active SSCHA transport is now intentionally rejected; see :doc:`input_files`.

The optional NEP GPU-batch path adds 24 public cases: 818 tests pass in total,
including 104 public tests. All 794 previous node IDs remain, with no skips or
deselection. Tests cover batch-size validation, unsupported combinations,
cache identity and interruption/recovery, invalid GPU output, missing
dependencies, real HiPhive fitting, and real FC2/FC3 fitting plus Slurm worker
collection (including cutoff-excluded pairs). GPU prediction is mocked in this
public suite; it does not establish GPU numerical accuracy or performance.
The opt-in ``tests/gpu/validate_nep_batch.py`` script compares actual CPUNEP and
GPU forces and FC2/FC3 on a GPU allocation, records versions/errors/timing,
and verifies cache reuse. Its output is a separate scientific validation artifact.

GPU validation on 2026-11-06 completed successfully on a Tesla V100-SXM2 32GB
(Slurm exit code 0), using Calorine 3.5, GPUMD 5.9.1, Phonopy/phono3py 4.1.0,
phonors 0.2.1 and NumPy 1.26.4. The public 3C-SiC model and a 2x2x2 primitive
supercell (16 atoms) produced 116 displacement force jobs. Maximum CPU/GPU
differences were 3.6113e-5 eV/Angstrom for forces, 3.1573e-4 eV/Angstrom^2
for FC2 and 1.9350e-3 eV/Angstrom^3 for FC3. These passed the preselected
absolute tolerances of 1e-4, 5e-3 and 0.2, respectively. A second GPU workflow
reused every cached force with the predictor disabled.

The finite-displacement stage took 1.89 seconds on CPU and 17.84 seconds with
GPU batches in this small case; GPU was slower. The script's first end-to-end
CPU time (205.12 seconds) also includes cold dependency initialization and must
not be compared with the later GPU time as a speedup benchmark. These are single
observations including I/O and fitting, not isolated force-kernel timings.
This validates the tested model and geometry only; it does not establish
converged transport results or performance for large systems. Slurm force-array
worker/collection integration is covered by the public mocked tests; the GPU
validation ran the direct FC2/FC3 workflow inside one real Slurm allocation.

RTA Slurm validation
----------------------

The distributed RTA path adds 22 public cases covering single temperatures and
temperature ranges, ordinary and Wigner transport, dependency chaining,
prepared/submitted status, automatic Rayon settings, and rejection of nested
Slurm transport in QHA/SSCHA temperature loops. The former local test requiring
RTA Slurm rejection now checks acceptance.

Two integration cases use the checked-in 3C-SiC NEP model to generate FC2/FC3
in a 2x2x2 primitive supercell, then run real phono3py 4.1.0 on a 3x3x3 mesh
at 200, 300 and 400 K. The generated prepare, array and collect shell scripts
are executed locally with simulated Slurm environment variables. Ordinary RTA
and SMM19 Wigner results match serial calculations with both relative and
absolute tolerances of 1e-10, including gamma, N/U, isotope and conductivity
datasets. A deliberately stale full-grid result is ignored, missing worker
data fails collection without replacing the published result, and reports
summarize only the published result. These checks do not submit real Slurm
jobs or establish cluster performance or convergence of the small mesh.

Source layout
---------------

- ``cli.py`` and ``command_registry.py``: public commands and terminal handling.
- ``config.py``: YAML reading, argparse coercion/defaults and parse entry points.
- ``config_schema.py``, ``config_plan.py``, ``config_validation.py`` and
  ``config_display.py``: accepted keys, stage selection, scientific checks and
  human-readable presentation, respectively.
- ``config_models.py``: typed snapshots of the mutable flat configuration.
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

``workflow.py`` retains the shared context used by stages and force workers.
It no longer runs plans or owns stage timing. See :doc:`architecture` for the
execution and configuration dependency graph.

Maintaining examples and the skill
------------------------------------

The public catalog intentionally contains nine 3C-SiC phonon/transport
workflows and one TD-BTE input. Explain their method differences in
``examples/README.md``. Use public input paths and placeholder executable
paths. Site-specific Slurm settings belong in an external batch script or
optional ``parallel`` settings.

When keys or templates change, update the input reference and
``.agents/skills/nepkappa-input/`` together. The example-catalog test validates
the nine phonon/transport inputs, six-section TD-BTE input, and their public assets.
Keep internal fixture inputs in ``tests/fixtures/`` aligned with the public
six-section ``run`` inputs and reduced direct-command inputs. Test both
acceptance of stage inputs and rejection of incomplete ``run`` inputs.

Documentation structure
-----------------------

The numbered user guide follows the seven chapters of the manual:
``introduction``, ``installation``, ``starting``, ``input_files``,
``harmonic_temperature``, ``thermal_transport`` and ``tdbte``. Keep this order
in ``index.rst``. Supplementary recipes, troubleshooting and developer notes
remain outside the numbered chapters.

Downloadable chapter inputs live in ``docs/source/inputs/`` and are displayed
with ``literalinclude`` so that the shown and downloaded YAML stay identical.
They complement the repository's ten public workflow templates. Check a changed
input with ``nepkappa info <file> --for <stage>`` for both its calculation stage
and ``plot``; these checks do not start a calculation. Use portable executable
placeholders and matching MPI/Slurm resource counts in cluster instructions.

Baseline snapshots
--------------------

The local ``benchmarks/data/`` includes compact data and original BAs source summaries
needed to check their hashes without private calculation directories.
These checks verify hashes and recorded values without rerunning DFT.
The archived baseline scope and local check commands are documented below.

.. toctree::
   :maxdepth: 1

   release_baseline
