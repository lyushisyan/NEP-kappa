# Reusable interatomic potentials

Potential files are grouped by material. A model being present here means only
that it is available to reproduce a workflow; it does not establish accuracy
outside its training domain.

| File or directory | Model species | Typical use |
| --- | --- | --- |
| `nep89_20250409.txt` | 89 elements | general-purpose NEP89 foundation model |
| `Si/` | Si | bulk and nanostructure examples |
| `3C-SiC/nep_3C-SiC.txt` | C, Si | material-specific 3C-SiC NEP |
| `BAs/` | B, As | BAs phonon and four-phonon calculations |
| `Cs2InAgCl6/` | Cs, In, Ag, Cl | double-perovskite/Wigner calculations |

Before a production calculation, verify the model species, training reference,
energy/force convention, intended phase, and validation errors. Large MACE
checkpoints are intentionally not duplicated here; calculation inputs should
reference the separately managed checkpoint or a named foundation model.

The bundled NEP89 copy has header `nep4_zbl 89` and SHA-256
`75168ece02e840e4a32644f982b78d43cba697f5b64b4c8134ab66c7a8c28be1`.
It was restored from the downloaded copy deployed on dongfang at
`/public1/home/wanggang/apps/NEP-kappa/potentials/nep89_20250409.txt`.

The material-specific 3C-SiC copy has header `nep4 2 C Si` and SHA-256
`ded60f644be08a5755f707e0888cae47f816d277569949d30270cfa79c0fc5c8`.
It was restored from the calculation copy deployed on dongfang at
`/public1/home/wanggang/apps/NEP-kappa/datasets/3C-SiC/nep_3C-SiC.txt`.
NEP89 and this model are distinct: coverage of C and Si alone does not make
a general model a material-specific model.
