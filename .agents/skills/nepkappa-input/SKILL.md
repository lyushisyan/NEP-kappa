---
name: nepkappa-input
description: Prepare, edit, explain, and validate NEP-kappa YAML inputs for phonons, thermal transport, QHA, Phonopy SSCHA, QHA+SSCHA, comparison, and convergence studies. Use when translating calculation requirements into inputs; preparing an input alone does not start calculations.
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

## Establish inputs and choose one route

Inspect provided files before asking for missing information. Establish the
structure, species, geometry, calculator/model, temperatures, supercells,
q mesh, result directory, and intended launch environment. For existing-result
requests, inspect available artifacts and choose a stage such as `kappa` or
`plot`; do not turn reuse into force-constant generation.

Use one minimal example. For a full workflow prefer a preset; small variants
such as RTA versus LBTE should edit the relevant key, not create duplicate
sections. Put user calculations under `calculations/<material-or-study>/`
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
`scph`, `qha-sscha`, or `plot`) for stage-only inputs. Comparison/convergence
use separate read-only parsers shown in the workflow reference. `report` and
`status` inspect results, not alternate YAML schemas.

Repair the reported cause and repeat failed checks without dropping requested
features. Check required files separately, distinguishing future stage outputs
from existing prerequisites. Do not require a model for transport-only reuse.
If the CLI is unavailable, state which checks could not run; do not install
large simulation dependencies solely to validate a draft.

Preparing an input does not itself authorize calculations. Do not invoke `run`,
`converge`, `compare`, `sbatch`, or load a calculator as a validation shortcut.
When calculation execution was requested, continue within that scope after
input checks.

## Handoff

Return the input link, launch directory, command, and material assumptions.
Separate parser validation, checked file/environment prerequisites, and remaining
numerical convergence work. A valid YAML is not a verified scientific result.
