---
name: nepkappa-input
description: Prepare, edit, explain, and validate NEP-kappa YAML inputs for phonons, thermal transport, QHA, SSCHA, QHA+SSCHA, existing-result plotting (including FC2-only), comparison, convergence studies, and time-dependent BTE. Use when turning calculation or plotting requirements into inputs or diagnosing input errors; preparing an input does not start calculations.
---

# NEP-kappa input assistant

Turn the user's calculation requirements into a minimal, version-compatible
input. Respond in their language and preserve their chosen material, calculator,
accuracy settings, and execution scope.

## Find the matching specification

Locate the checkout supplied by the user or open in the workspace. When this
skill lives at `.agents/skills/nepkappa-input`, the repository is three parents
above it. Verify `src/nepkappa/config.py` and `examples/`; a copied skill may
need a separately supplied repository path.

Check `nepkappa --version` and, when relevant, `nepkappa --help` and
`nepkappa validate --help`. Match the installation to the sources:

```bash
python -c 'import nepkappa; print(nepkappa.__version__); print(nepkappa.__file__)'
```

Use `config.py` and `command_registry.py` for accepted YAML keys and targets;
`docs/source/input_files.rst` for units; `examples/README.md` for templates.
Internal `config_models.py` field names are not necessarily YAML keys.
Read the relevant stage implementation when runtime behavior is unclear.

Read [workflow guidance](references/workflows.md) for route selection,
cutoffs, temperature conventions, and analysis-specific validation. Treat
that reference as guidance; the target version's parser is authoritative.
For plotting or time-dependent BTE, also read
[artifact-reuse guidance](references/analysis.md). Do not infer development
feature availability from the version number alone; inspect the actual checkout.
For VASP location, executable changes, POTCAR setup, or MPI/Slurm environment
questions, read [VASP environment guidance](references/vasp-environment.md).
Help check and configure the target host, not just YAML syntax.

## Establish inputs and choose one route

Inspect provided files before asking for missing information. Establish the
structure, species, geometry, calculator/model, temperatures, supercells,
q mesh, result directory, and intended launch environment **only as needed for
the requested stage**. First distinguish input preparation, result inspection,
plot generation, and calculation execution. For existing-result
requests, inspect available artifacts and choose a stage such as `kappa` or
`plot`; do not turn reuse into force-constant generation.
Do not require a potential or FC3 for FC2-only plotting, or a calculator for
TD-BTE using existing force constants or a kernel. Inspect available artifacts
before promising plots.

Use one minimal example. For static calculations, start from the six-section input
in `examples/nep-rta-wigner-3ph.yaml`: structure, calculator, force-constant, kappa, plot,
and output. `structure.relaxation` controls optimization. QHA, SSCHA, and
four-phonon switches belong inside `force-constant` and select their full
static workflows; read `docs/source/input_files.rst` for supported combinations.
The four-phonon switch runs both FC4 and conductivity. For dynamic TD-BTE,
use the separate `tdbte` and `output` shape in `examples/tdbte.yaml`.
Omit unneeded parameter lines rather than inventing a new schema.
Set `kappa.engine: phono3py` for three-phonon RTA/LBTE or
`kappa.engine: fourphonon` with `force-constant.four-phonon: true`.
For FourPhonon, use separate `kappa.method-3ph` and `kappa.method-4ph`
settings and check supported pairings in `docs/source/input_files.rst`.
Place optional Slurm/MPI/OpenMP settings in top-level `parallel.force-constant`
or `parallel.kappa` rather than lengthening the scientific sections. For QHA
thermal expansion to affect 3ph RTA, use `kappa.qha-volumes: true`; this runs
FC2/FC3 and transport anew at each target QHA volume.
For a VASP harmonic calculation, use `nepkappa fc2` on a suitable VASP
six-section input, then `nepkappa plot` after FC2 exists. Run `nepkappa relax`
first if `structure.relaxation: true`; `nepkappa fc2fc3` computes both force
constant orders without invoking transport.
Advanced stage-first `workflow.stages` inputs remain valid even though the
public catalog has nine static workflows. Existing `workflow.preset` and
`workflow.steps` inputs remain valid; do not combine either with
`workflow.stages`. Small variants such as RTA versus LBTE should edit the
transport method, not create duplicate sections. Put user calculations under
`calculations/<material-or-study>/`
unless they choose another location. `examples/` is the public template catalog.

Model identity is not established by its filename or directory. For readable
NEP models inspect the header and ensure the structure's species are covered.
`potentials/3C-SiC/nep_3C-SiC.txt` is material-specific;
`potentials/nep89_20250409.txt` is the general 89-element model. Confirm origin
or hashes when selecting among copies; use `potentials/README.md` and adjacent
model notes. Species compatibility alone does not establish phase accuracy.
Do not replace a missing model with a renamed or unrelated potential.

Films require physical thickness and wires require cross-sectional area; these
cannot be inferred from vacuum lengths alone. Remote paths not checked on their
host are unverified. Do not copy personal cluster accounts, module commands,
partitions, or resource allocations from historical inputs.

## Write and validate

Preserve explicit numerical choices. Identify unconverged defaults as starting
values, and keep essential unresolved fields visible in a labeled draft.
Use a distinct output directory when changing the material or study.

Ordinary input paths resolve from the **launch directory**, not the YAML
location. Convergence `base` and `study.directory` resolve from the study YAML.
State the launch directory in the handoff.

Use `submit: false` / `study.execute: false` for input-only drafts, while
preserving separately authorized execution and explicit user settings.
For ordinary inputs, validate the intended command without computing:

```bash
nepkappa validate input.yaml --for run
nepkappa info input.yaml --for run
```

Replace `run` with the actual stage (`fc2fc3`, `kappa`, `kappa4`, `qha`,
`scph`, `qha-sscha`, `tdbte`, or `plot`) for stage-only inputs. Comparison/convergence
use separate read-only parsers shown in the workflow reference. `report` and
`status` inspect results, not alternate YAML schemas.

Repair the reported cause and repeat failed checks without dropping requested
features. Check required files separately, distinguishing future stage outputs
from existing prerequisites. Do not require a model for transport-only reuse.
If the CLI is unavailable, state which checks could not run; do not install
large simulation dependencies solely to validate a draft.

Preparing an input does not itself authorize calculations. Do not invoke `run`,
`converge`, `compare`, `plot`, `tdbte`, `sbatch`, or load a calculator as a
validation shortcut. Adding a `plot` section does not add a plotting stage;
give the explicit stage command or enable `analysis.method: plot` in an
authorized stage-first plan.
When calculation execution was requested, continue within that scope after
input checks.

## Handoff

Return the input link, launch directory, command, and material assumptions.
Separate parser validation, checked file/environment prerequisites, and remaining
numerical convergence work. A valid YAML is not a verified scientific result.
Keep the handoff short: what was prepared, where to launch, one next command,
what was checked, and what is still missing. Explain settings in ordinary
language; do not make users learn all presets or commands before helping them.
