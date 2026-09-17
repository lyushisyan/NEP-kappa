# Examples

Choose one template, copy it to your own input, and run from the repository
root. These are workflow demonstrations, not converged production settings.
Outputs are stored under `calculations/example-runs/`.

## Common workflows

| Input | What it demonstrates | Preparation |
| --- | --- | --- |
| [nep-rta.yaml](nep-rta.yaml) | Bulk Si three-phonon RTA | Bundled Si structure and NEP |
| [vasp-rta.yaml](vasp-rta.yaml) | DFT relaxation and forces | Set VASP and POTCAR paths |
| [mace-rta.yaml](mace-rta.yaml) | MACE on CPU | Install MACE; set a checkpoint |
| [wigner.yaml](wigner.yaml) | Wigner transport | Bundled Si demonstration |
| [qha.yaml](qha.yaml) | Isotropic thermal expansion | Check volume range and stability |
| [sscha.yaml](sscha.yaml) | Fixed-volume Phonopy SSCHA + transport | Converge sampling and iterations |
| [qha-sscha.yaml](qha-sscha.yaml) | QHA-volume SSCHA; independent 3ph/4ph switches | FourPhonon needed when 4ph enabled |
| [bas-three-four-phonon.yaml](bas-three-four-phonon.yaml) | BAs 3ph and 3ph+4ph | Thirdorder, Fourthorder, FourPhonon |

For any ordinary workflow:

```bash
cp examples/nep-rta.yaml input.yaml
nepkappa validate input.yaml
nepkappa info input.yaml
nepkappa run input.yaml
```

Set a distinct `output.result_dir` before changing the material or physical
settings. Use `nepkappa plot input.yaml` for completed phono3py results and
`nepkappa report input.yaml` for a summary.

## Advanced templates

| Input | Purpose |
| --- | --- |
| [nep-hiphive.yaml](nep-hiphive.yaml) | Fit FC2/FC3 with HiPhive |
| [film.yaml](film.yaml) | Film geometry and effective thickness |
| [thirdorder.yaml](thirdorder.yaml) | FC2/FC3 for ShengBTE; `fc2fc3` only |
| [vasp-fc4.yaml](vasp-fc4.yaml) | VASP/Fourthorder FC4; `fc4` only |
| [slurm-vasp.yaml](slurm-vasp.yaml) | VASP force arrays and continuation |
| [slurm-lbte.yaml](slurm-lbte.yaml) | Distributed LBTE |
| [fourphonon.yaml](fourphonon.yaml) | FourPhonon with Slurm transport |

Slurm templates need your site's environment and resources. All start with
`submit: false`. A workflow can still execute preceding local stages.
For existing IFCs, use `nepkappa info examples/fourphonon.yaml --for kappa4`
and the `kappa4` command to reuse them.

## Analysis inputs

These use separate schemas and commands:

- [compare.yaml](compare.yaml): `nepkappa compare examples/compare.yaml`;
  keep two or more datasets pointing to completed calculations.
- [converge-qmesh.yaml](converge-qmesh.yaml):
  `nepkappa converge examples/converge-qmesh.yaml`; `study.execute: false`
  prepares cases and analysis without running simulations. Its `base` and
  study directory resolve relative to the study YAML.

## Small changes do not need another template

Edit the existing section of the copied YAML; do not append duplicate keys.

| Variant | Change |
| --- | --- |
| Serial LBTE | In `nep-rta.yaml`, set `kappa.method: lbte` |
| Film + HiPhive | Add fitting settings from `nep-hiphive.yaml` to `film.yaml`; choose valid cutoffs |
| VASP + HiPhive | Combine the VASP calculator with the fitting section |
| VASP + Thirdorder | Replace the NEP calculator in `thirdorder.yaml` with configured VASP settings |
| Two-model comparison | Keep two entries in `compare.yaml` |

Numbered filenames now have descriptive names. Duplicate variants have been
consolidated; the removed local inputs are retained in
`calculations/legacy-examples/`.
Structures are in `structures/Si/` and `structures/BAs/`. Models are maintained
in [../potentials/](../potentials/README.md).
