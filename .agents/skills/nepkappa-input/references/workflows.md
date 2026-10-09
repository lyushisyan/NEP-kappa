# Workflow selection and constraints

Paths in this reference are relative to the target NEP-kappa checkout. These
are routing notes for the current repository, not a second schema; consult its
parser and input documentation for exact accepted fields and version changes.

## Starting examples

| Request | Starting point in `examples/` | Intended command |
| --- | --- | --- |
| VASP harmonic FC2 and phonon plots | `vasp-rta-3ph.yaml`; configure executable and POTCAR | `stage relax`, `stage fc2`, `plot` |
| VASP QHA only | Adapt a six-section input with `force-constant.qha` settings | `stage qha` |
| VASP 3ph RTA | `vasp-rta-3ph.yaml`; configure executable and POTCAR | `run` |
| NEP 3ph RTA + Wigner | `nep-rta-wigner-3ph.yaml` | `run` |
| NEP 3ph LBTE + Wigner | `nep-lbte-wigner-3ph.yaml` | `run` |
| NEP 3ph+4ph RTA | `nep-rta-3ph-4ph.yaml` | `run` |
| QHA-volume 3ph RTA | `nep-qha-rta-3ph.yaml`; regenerate FC2/FC3 at each QHA volume | `run` |
| QHA-volume SSCHA 3ph RTA | `nep-qha-sscha-rta-3ph.yaml` | `run` |
| 3ph+4ph RTA with Wigner coherence | `nep-rta-wigner-3ph-4ph.yaml`; Wigner_Park executable | `run` |
| 3ph LBTE + 4ph RTA | `nep-lbte-3ph-rta-4ph.yaml` | `run` |
| 3ph LBTE + 4ph LBTE | `nep-lbte-3ph-lbte-4ph.yaml` | `run` |
| Existing FC2 / transport plots | Minimal plotting sections; see `references/analysis.md` | `plot` |
| Time-dependent BTE | `tdbte.yaml`; existing FC2/FC3 + metadata, or a kernel | `stage tdbte` or `run` |

The public catalog is limited to these nine static inputs and TD-BTE. Other
supported calculations can use a copied six-section input with changed settings;
do not refer to removed example filenames. Site-specific Slurm resources belong
in an external batch script or optional top-level `parallel` settings.

For existing FC2/FC3 transport, adapt only the necessary settings from a
transport example and use `stage kappa`. For harmonic-only generation use `stage fc2`.
For plotting existing outputs use `plot`. Do not let a copied preset change a
request to reuse existing results into one that regenerates force constants.

## Cross-section checks

- New static inputs use six top-level sections: `structure`, `calculator`,
  `force-constant`, `kappa`, `plot`, and `output`. Put QHA, SSCHA, and four-phonon
  switches and their options inside `force-constant`; 4PH means FC4 generation
  plus conductivity. TD-BTE uses separate `tdbte` and `output` sections. Older
  `workflow.stages`, presets, and custom steps remain valid for compatibility.
  Do not mix nested static switches with an explicit workflow. QHA alone does
  not automatically feed corrected FC2 into transport; use
  `kappa.qha-volumes: true` for a new FC2/FC3/RTA calculation at each
  temperature's QHA equilibrium volume. Set `kappa.engine` explicitly in new
  examples. Put optional parallel settings in top-level `parallel` mappings.
- For FourPhonon static inputs, set both `kappa.method-3ph` and
  `kappa.method-4ph`. Supported pairs are `rta/rta`, `lbte/rta`, and
  `lbte/lbte`; `rta/lbte` is unavailable. Wigner_Park requires `rta/rta`
  and a Wigner_Park executable. Do not also set `kappa.method` or a
  conflicting `force-constant.four-phonon.solver`.
- `kappa.temps` uses `[T]` or `[T_min, T_max, T_step]`, not an arbitrary list of
  three temperatures. Check the step is positive and the range is ordered even
  if an installed parser only checks length. Read the relevant section before
  applying temperature conventions to QHA, SCPH, or FourPhonon.
- Phono3py `kappa` needs `phono3py_disp.yaml`, `fc2.hdf5`, and `fc3.hdf5` in
  the result directory. FC2/FC3 export must be `phono3py` or `both` for that
  route. A FourPhonon input makes the same `nepkappa kappa` command use
  ShengBTE-format force constants; FC4 generation alone does not calculate
  four-phonon thermal conductivity.
- With relaxation enabled, direct `fc2`/`fc2fc3` expects
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
- Use `cutoff-fc3` for native phono3py FC3 (positive Angstrom; omit for no cutoff).
  The old `pair-cutoff-fc3` is only a compatibility alias, not a key for new inputs.
  Thirdorder also uses
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
- In the static form, `force-constant.sscha.run-transport` and
  `force-constant.four-phonon.enabled` choose the three- and four-phonon routes
  within QHA+SSCHA. The latter needs Fourthorder and FourPhonon. Outside the coupled route, use
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
