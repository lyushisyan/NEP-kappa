# VASP executable and environment setup

Use for VASP path discovery, moving calculations between hosts, and preparing
VASP inputs. Inspect `src/nepkappa/adapters/vasp.py`, the calculator section of
`docs/source/input_files.rst`, and `examples/vasp-rta-3ph.yaml` in the matching
checkout. VASP and its licensed POTCAR library are external installations,
not bundled with NEP-kappa.

## Locate on the intended host

Establish the target host first. A local `command -v vasp_std` says nothing
about a remote server. Inspect the user's existing input and batch scripts;
historical paths are candidates, not evidence of current availability. With
authorized access, use read-only checks on that host:

```bash
hostname
command -v vasp_std
command -v mpirun
command -v srun
```

If the site uses environment modules, inspect `module list` and available
VASP/compiler/MPI modules in its initialized shell. A missing PATH entry does
not prove VASP is uninstalled. Check user-provided or site-documented locations
with `ls -l` / `test -x`; avoid searching the entire filesystem. Do not execute
VASP (even with `--version`) as a discovery shortcut on a login node. Do not
install, upload licensed software, change shell startup files, or submit jobs
as an implicit consequence of a path question.

## Configure explicitly

Resolution priority is `calculator.vasp_command`, then `calculator.vasp_path`,
then automatic discovery. Current discovery checks historical `/root/software`
locations before PATH. Those defaults are not portable server configuration.
Prefer an explicit command for reproducible MPI execution:

```yaml
calculator:
  name: vasp
  vasp_command: "srun /absolute/path/to/vasp_std"
  potcar_path: /absolute/path/to/potpaw_PBE
```

This is a template, not a claim that this launcher works at every site. Use
the site's supported launcher and MPI library for the actual VASP build.
If both command and path exist, changing only `vasp_path` has no effect;
update the binary inside `vasp_command` too, or remove the redundant path.
Setting only `vasp_path` does not add MPI parallelism.

The adapter uses `shlex.split`, not a shell. Do not put `module load ... &&`,
`export`, redirections, `$VARIABLE`, or `~` expansion in `vasp_command`.
Use absolute paths (quote paths with spaces inside the command). Environment
setup belongs in the launch shell/batch script or the force jobs'
`force-constant.parallel.preamble`. That preamble does not configure every
other workflow stage; outer relaxation/QHA execution needs its own environment
and allocation. Check generated scripts and avoid nested MPI launchers.

Match launcher ranks to Slurm `ntasks`, and OpenMP threads to `cpus-per-task`.
Never copy a historical partition/account/core count unverified. Changing an
executable path is not permission to change ENCUT, k points, KPAR, or NCORE.

## POTCAR and validation boundaries

`potcar_path` accepts a ready combined POTCAR or a library directory. Check
readability, species order, functional, and potential variants/version. For a
library, inspect which per-element files the current adapter will select; do
not silently substitute `_sv`/`_pv` or another functional. A supplied combined
POTCAR must already match the structure; parser acceptance does not check this.
Inspect limited header metadata, not full licensed contents. Never commit
POTCAR data into public examples or distribute it without authorization.

Do not call `resolve_potcar` merely to inspect: library assembly writes files.
Run `nepkappa info input.yaml --for <stage>` for syntax and plan checks only.
Report separately: configured command, binary existence/executable bit,
environment/MPI checks, POTCAR checks, and whether a compute-node smoke test
has actually run. Remote or compute-node access not exercised remains unverified.
