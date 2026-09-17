# Release notes

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
