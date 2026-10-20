# Release notes

## 2.0.1 — 2026-10-20 (reissued)

This release replaces the September 23, 2026 `v2.0.1` release and tag.

### Highlights

- Fix the shared coordinate and lattice-translation conventions for FC2, FC3,
  and FC4 generation and export.
- Compare three-phonon and three-plus-four-phonon conductivity with the same
  FourPhonon backend, with independent RTA/LBTE settings and overlaid plots.
- Add FourPhonon Wigner RTA transport with particle, coherence, and total
  conductivity outputs.
- Build TD-BTE collision kernels automatically from FC2/FC3, stream event
  chunks during integration, and plot energy evolution in a six-panel figure.
- Unify CLI commands and YAML inputs; simplify terminal output while retaining
  detailed logs, and align documentation and examples with the user manual.

### Migration

The public commands are `info`, `run`, `relax`, `fc2`, `fc2fc3`, `qha`, `kappa`,
`tdbte`, `plot`, and `report`. Legacy commands and input schemas have been
removed. Start from the updated examples and check inputs with `nepkappa info`.
Regenerate matched force constants when older files use inconsistent coordinate
or lattice-translation conventions.

Verification: 691 local tests passed and the strict Sphinx documentation build
passed. TD-BTE absolute-rate validation and numerical convergence remain
necessary for quantitative physical interpretation.

### Earlier changes included from the September 23 snapshot

- Unified FC3 cutoff input as `cutoff-fc3`, retaining the old phono3py alias.
- Added FC2-only plotting and improved result reports.
- Clarified SSCHA and QHA-volume approximations; added optional diagonal
  on-shell bubble frequency corrections (without updating conductivity).
- Added experimental time-dependent BTE from a prebuilt energy-shell kernel.
  Physical-rate validation and automatic kernel construction remain incomplete.
- Updated documentation and the input-assistant Skill, including VASP setup.

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
