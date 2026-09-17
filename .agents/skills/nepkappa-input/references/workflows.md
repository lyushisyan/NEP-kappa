# Workflow selection and constraints

Paths in this reference are relative to the target NEP-kappa checkout. These
are routing notes for the current repository, not a second schema; consult its
parser and input documentation for exact accepted fields and version changes.

## Starting examples

| Request | Starting point in `examples/` | Intended command |
| --- | --- | --- |
| NEP bulk RTA / serial LBTE | `nep-rta.yaml`; choose `kappa.method` | `run` |
| HiPhive fitting | `nep-hiphive.yaml` | `run` |
| Wigner transport | `wigner.yaml` | `run` |
| VASP / MACE | `vasp-rta.yaml` / `mace-rta.yaml` | `run` |
| Film geometry | `film.yaml`; combine fitting settings if requested | `run` |
| Slurm LBTE / force jobs | `slurm-lbte.yaml` / `slurm-vasp.yaml` | `kappa` / `fc2fc3`, or full `run` |
| FC4 generation only | `vasp-fc4.yaml` | `fc4` |
| Thirdorder FC3 | `thirdorder.yaml`; adapt calculator for VASP | `fc2fc3` |
| Isotropic QHA | `qha.yaml` | `run` or `qha` |
| Fixed-volume Phonopy SSCHA | `sscha.yaml` | `run` or `scph` |
| QHA-volume Phonopy SSCHA | `qha-sscha.yaml` | `run` or `qha-sscha` |
| Both 3ph and 3ph+4ph | `bas-three-four-phonon.yaml` | `run` |
| FourPhonon / existing FCs | `fourphonon.yaml` | `run` / `kappa4` |
| Two or more model results | `compare.yaml` | `compare` |
| q-mesh convergence | `converge-qmesh.yaml` | `converge` |

For existing FC2/FC3 transport, adapt only the necessary settings from a
transport example and use `kappa`. For harmonic-only generation use `fc2`.
For plotting existing outputs use `plot`. Do not let a copied preset change a
request to reuse existing results into one that regenerates force constants.

## Cross-section checks

- Normal full workflows use presets `three-phonon`, `four-phonon`, `qha`,
  `scph`, or `qha-sscha`. Explicit stage lists require
  `workflow.preset: custom`. Resolve the
  actual plan in the target parser; the SCPH preset expands to `fc2fc3 -> scph`.
- `kappa.temps` uses `[T]` or `[T_min, T_max, T_step]`, not an arbitrary list of
  three temperatures. Check the step is positive and the range is ordered even
  if an installed parser only checks length. Read the relevant section before
  applying temperature conventions to QHA, SCPH, or FourPhonon.
- `kappa` needs `phono3py_disp.yaml`, `fc2.hdf5`, and `fc3.hdf5` in the result
  directory. FC2/FC3 export must be `phono3py` or `both` if generation precedes
  this stage. FourPhonon uses ShengBTE-format force constants; FC4 generation
  alone does not calculate four-phonon thermal conductivity.
- With relaxation enabled, stage-only `fc2`/`fc2fc3` expects
  `POSCAR_relaxed` in the result directory. It is a future output when `relax`
  precedes these stages in the full plan.
- External ASE/plugin calculators, including MACE, support relaxation through
  `stages/structure.py`. Check the calculator's energy/force/stress capabilities
  before proposing cell relaxation. For an
  already-relaxed structure, use `relaxation.enabled: false`. Do not transplant
  VASP-specific relaxation settings to other calculators.
- Films need user-established effective thickness (angstrom); wires need
  effective area (angstrom squared). Keep supercells and q meshes consistent
  with the identified periodic axes (normally 1 along a vacuum direction).
  Geometry normalization currently affects plots, not original kappa HDF5.
- HiPhive cutoffs and displacements are material-dependent. Negative
  Thirdorder/Fourthorder cutoffs follow neighbor-shell conventions.
- Native phono3py FC3 accepts `pair-cutoff-fc3` in Angstrom. Thirdorder uses
  `cutoff-fc3` (negative neighbor shell or positive nm), while HiPhive uses
  `cutoffs` in Angstrom. Native finite-displacement FC2 has no independent
  real-space cutoff; its interaction range is controlled by `dim-fc2`.
- QHA needs at least five strictly increasing volume ratios and a sensible
  starting structure near its zero-temperature equilibrium volume. Neither a
  passing parser nor a template volume range establishes dynamical stability.
- Phonopy SSCHA needs a harmonic FC2 and an ASE calculator providing energy
  and forces; when `run-transport` is enabled it additionally needs matching
  phono3py metadata and FC3.
  The command-driven `calculator.name: vasp` backend is not supported by this
  SSCHA implementation; QHA itself supports VASP.
- The `scph` command selects Phonopy stochastic SSCHA, not the separate
  ALAMODE perturbative SCPH method. State the actual approximation in reports.
- `qha-sscha.three-phonon` and `qha-sscha.four-phonon` are independent switches.
  The latter needs Fourthorder and FourPhonon. Outside the coupled route, use
  the BAs custom plan `[fc2fc3, kappa, fc4, kappa4]` for both channels. `kappa4`
  means combined 3ph+4ph transport, not conductivity from 4ph scattering alone.
- `qha-sscha` interpolates each requested SSCHA temperature inside the
  completed QHA range, then regenerates matching FC2 and, for transport, FC3.
  It does not support nested force-job Slurm arrays; submit the whole run as
  one outer batch job.
- Automatically generated FourPhonon CONTROL files disable non-analytic
  corrections. Polar-material calculations that need those corrections
  require an appropriate custom CONTROL with dielectric/Born-charge data;
  do not invent this data. Espresso harmonic input also needs a matching
  custom CONTROL. Report this separately from YAML parser acceptance.
- Slurm settings are stage-specific: force jobs, LBTE, and FourPhonon have
  distinct `parallel` sections. Use the corresponding example; check MPI
  processes and OpenMP threads agree with the requested allocation. Ask for
  missing site-specific settings.

## Comparison inputs

Comparison YAML has `datasets`, `compare`, and `plot` sections (with optional
geometry). `datasets` lists labels and result directories. Use completed
results with `phono3py_disp.yaml`, `fc2.hdf5`, and `kappa-m*.hdf5`; inspect the
files when accessible. It is not a normal workflow YAML.

Neither `validate --for compare` nor `info --for compare` is supported.
Use the installed package's read-only parser with the matching Python:

```bash
python -c 'import sys; from nepkappa.config import parse_compare_args; parse_compare_args(sys.argv[1]); print("Comparison configuration parsed successfully")' compare.yaml
```

Do not invoke `nepkappa compare` just to validate: it writes plots and results.

## Convergence-study inputs

Use `base`, `parameter`, `values`, and `study`, following the convergence
example. `base` and `study.directory` resolve relative to the **study YAML's
directory**. Paths inside the base workflow still use the launch directory.
Make this explicit when the input and study files are in different folders.

The swept dotted parameter must already exist in the base YAML. There must
be at least two distinct values. Set `study.execute: false` for input-only
preparation unless the user explicitly requests execution.

Neither `validate --for converge` nor `info --for converge` is supported.
Use the matching Python environment:

```bash
python -c 'import sys; from nepkappa.convergence import parse_convergence_args; parse_convergence_args(sys.argv[1]); print("Study configuration parsed successfully")' converge.yaml
```

This parser checks the study structure and base-file existence, but does not
validate the dotted path or all generated workflows. Read the base YAML and
verify the path exists; follow the parser's hyphen/underscore normalization.
For full input validation, substitute each value into a temporary copy of the
base and run workflow validation for the intended command from the same
launch directory. Do not alter the original base during these checks.

Do not use `nepkappa converge` as a read-only validator: even with
`execute: false`, it creates case inputs and analysis outputs. Run it only
when the user's request includes preparing the study outputs or execution.
The last configured value is the comparison reference; do not present it as
the exact physical answer or claim convergence from incomplete calculations.
