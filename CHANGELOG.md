# Release notes

## Unreleased

- Expose ten direct commands: ``info``, ``run``, ``relax``, ``fc2``,
  ``fc2fc3``, ``qha``, ``kappa``, ``tdbte``, ``plot``, and ``report``.
  The input selects the transport backend and temperature-coupled route.
- Keep ``nepkappa info`` as the sole input-inspection command; it checks and
  displays the selected workflow. Remove the ``init`` and ``validate`` CLI
  commands; copy an example YAML to start a new input.
- Remove the ``stage`` and ``status`` command entries. Earlier direct-stage,
  comparison, and convergence commands beyond the ten listed above are no
  longer accepted.
- Add a stage-first ``workflow.stages`` interface with per-stage method and
  feature choices. Static inputs now use six sections with QHA, SSCHA, and
  four-phonon switches inside ``force-constant``. TD-BTE uses a separate ``tdbte``/``output``
  input. Older preset and stage-first inputs remain supported.
- Add an optional FourPhonon Wigner_Park RTA path for combined three/four-phonon
  population, coherence, and total conductivity, with tensor consistency checks.
- Expose independent `kappa.method-3ph` and `kappa.method-4ph` choices for
  FourPhonon static inputs, and add 3C-SiC Wigner RTA, 3ph LBTE/4ph RTA,
  and 3ph LBTE/4ph LBTE examples.
- Limit the public `examples/` YAML catalog to nine 3C-SiC static
  workflows and the separate TD-BTE input; update documentation links.
- Add explicit `kappa.engine` selection, an optional top-level `parallel` section,
  and 3C-SiC calculation inputs, including a bundled primitive structure.
  QHA-coupled 3ph RTA now regenerates
  FC2/FC3 and conductivity at each target temperature's QHA volume.
- Plot RTA Normal and Umklapp rates together, Wigner particle/coherence/total
  conductivity, and separate three-/four-phonon scattering rates. Compare only
  conductivity solutions actually present in phono3py or FourPhonon outputs.
- Build time-dependent BTE energy-shell kernels directly from matching FC2/FC3,
  phono3py metadata, and a q mesh.
- Stream checksummed event chunks during integration, with optional Numba
  acceleration and compatibility with existing single-file kernels.
- Record kernel provenance, build failures, and numerical audits in reports.
- Remove the TD-BTE `experimental` input switch; retain physical-validation
  limitations and conservation, equilibrium, and entropy checks.
- Update YAML examples, documentation, and the input-preparation Skill.

Local verification: 459 tests passed with a fresh Numba cache; strict
documentation build and a real SiC small-grid TD-BTE build-and-solve smoke test
passed. The FourPhonon Wigner_Park interface has not yet been checked against a
completed external-solver run. Absolute-rate validation and mesh convergence
remain separate from these implementation checks.

## 2.0.1 — 2026-09-23

- Unified FC3 cutoff input as `cutoff-fc3`, retaining the old phono3py alias.
- Added FC2-only plotting and improved result reports.
- Clarified SSCHA and QHA-volume approximations; added optional diagonal
  on-shell bubble frequency corrections (without updating conductivity).
- Added experimental time-dependent BTE from a prebuilt energy-shell kernel.
  Physical-rate validation and automatic kernel construction remain incomplete.
- Updated documentation and the input-assistant Skill, including VASP setup.

Local verification: 348 tests passed and the strict documentation build passed.
Research kernels, local tests, and benchmark datasets are not part of this release.

## 2.0.0 — 2026-09-17

NEP-kappa 2.0 unifies phonon and lattice thermal-transport workflows using
first-principles and machine-learning force providers.

### Workflows

- Workflow presets and an interactive `nepkappa init` input generator.
- NEP, VASP, native MACE, external ASE factories, and calculator plugins.
- FC2/FC3 through finite displacement or HiPhive, Thirdorder FC3, and
  Fourthorder FC4 generation.
- phono3py RTA, iterative LBTE, and SMM19 Wigner transport.
- FourPhonon three-plus-four-phonon transport and a template computing both
  pure three-phonon and combined three-plus-four-phonon results.
- Isotropic QHA, fixed-volume Phonopy stochastic SSCHA, and QHA+SSCHA with
  independent three/four-phonon transport switches.
- Multi-model comparison, parameter convergence studies, and result reports.

### Execution and reproducibility

- Stage-specific configuration validation, focused inspection, and concise
  errors with optional debug tracebacks.
- Slurm force arrays, distributed LBTE, FourPhonon submission, and job status.
- Input-fingerprinted force caches, provenance, and recorded workflow state.
- Separate application, stage, calculator, scheduler, and artifact interfaces.

### Getting started and migration

- Descriptive example filenames replace numbered filenames. See
  [the catalog](examples/README.md) for the consolidated templates.
- Structures are in `examples/structures/`; models are in `potentials/`.
- Material calculations and outputs belong in the local `calculations/`
  workspace, whose contents are ignored except for its README.
- `tests/` and `benchmarks/` are maintained locally and excluded from Git.
  GitHub automation builds documentation without depending on these files.
- Updated quick-start, workflow, parameter, and development documentation.
- Repository-local `nepkappa-input` skill for preparing and validating YAML.

The `scph` command selects Phonopy stochastic SSCHA, not ALAMODE SCPH.
Temperature-renormalized transport uses the force-constant approximation
documented for each route. Example settings require material-specific
convergence checks. VASP, Thirdorder, Fourthorder, FourPhonon, and optional
machine-learning backends need their respective external installations.
