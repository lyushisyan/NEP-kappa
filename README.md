# NEP-kappa

[![Documentation](https://readthedocs.org/projects/nep-kappa/badge/?version=latest)](https://nep-kappa.readthedocs.io/en/latest/)

**Phonons and lattice thermal transport with first-principles and machine-learning potentials.**

NEP-kappa connects NEP, VASP, MACE, and ASE calculators to force-constant generation,
Phonopy/phono3py, and FourPhonon. Describe a calculation in YAML and run it with
`nepkappa run input.yaml`.

Current version: **2.0.0**. See the [release notes](CHANGELOG.md).
[Quick start](docs/source/starting.rst) · [Examples](examples/README.md) ·
[Input reference](docs/source/input_files.rst) · [中文上手](#中文上手)

## Start here

The package requires Python 3.9 or newer; Python 3.11 is a practical starting
environment. From a terminal:

```bash
git clone https://github.com/lyushisyan/NEP-kappa.git
cd NEP-kappa
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
nepkappa --version
```

Run the bundled bulk-Si example from the **repository root**:

```bash
nepkappa validate examples/nep-rta.yaml
nepkappa run examples/nep-rta.yaml
nepkappa plot examples/nep-rta.yaml
nepkappa report examples/nep-rta.yaml
```

Results go to `calculations/example-runs/nep-rta/`. This demonstrates the workflow;
its supercell and q mesh are starting settings, not a convergence claim.
Validation checks the configuration, not every external program or the accuracy
of a potential.

For your own material, create a separate input:

```bash
nepkappa init input.yaml
nepkappa validate input.yaml
nepkappa run input.yaml
```

Use your material's structure and matching potential. Ordinary input paths are
relative to the directory where you launch the command.

## Choose a calculation

| Goal | Start from | Backend / prerequisite |
| --- | --- | --- |
| Three-phonon RTA or LBTE | [nep-rta.yaml](examples/nep-rta.yaml) | NEP + phono3py; change `kappa.method` for LBTE |
| DFT forces and transport | [vasp-rta.yaml](examples/vasp-rta.yaml) | Your VASP executable and licensed POTCAR library |
| MACE forces and transport | [mace-rta.yaml](examples/mace-rta.yaml) | `pip install -e '.[mace]'` and a checkpoint |
| Wigner transport | [wigner.yaml](examples/wigner.yaml) | phono3py SMM19 |
| Thermal expansion (QHA) | [qha.yaml](examples/qha.yaml) | Phonopy volume scan |
| Fixed-volume renormalization | [sscha.yaml](examples/sscha.yaml) | Phonopy stochastic SSCHA; CLI name `scph` |
| QHA + SSCHA and optional transport | [qha-sscha.yaml](examples/qha-sscha.yaml) | SSCHA at QHA equilibrium volumes |
| Both 3ph and 3ph+4ph conductivity | [bas-three-four-phonon.yaml](examples/bas-three-four-phonon.yaml) | Thirdorder, Fourthorder, FourPhonon |
| Compare models / check convergence | [Example catalog](examples/README.md) | Completed results / a base input |

The catalog also covers HiPhive fitting, films, FC3/FC4 generation, and Slurm.
Optional executables are installed separately; see
[installation](docs/source/installation.rst). SSCHA here is the Phonopy route
implemented in this package. Transport with renormalized FC2 does not by itself
include every higher-order anharmonic correction.

## The commands most users need

| Command | Purpose |
| --- | --- |
| `nepkappa init input.yaml` | Create an input interactively |
| `nepkappa validate input.yaml` | Check input syntax and supported settings |
| `nepkappa run input.yaml` | Run the selected workflow |
| `nepkappa status input.yaml` | Inspect stored job state and the Slurm queue |
| `nepkappa plot input.yaml` | Plot completed phono3py outputs |
| `nepkappa report input.yaml` | Write a result summary |

Stage commands allow restarts and reuse of existing force constants.
See `nepkappa --help` and the [workflow guide](docs/source/tutorial.rst).
Slurm templates start with `submit: false`; configure the cluster and inspect
generated scripts before enabling submission.

## Input assistant

The repository includes [nepkappa-input](.agents/skills/nepkappa-input/SKILL.md).
In a compatible assistant opened in this checkout, ask:

> Use $nepkappa-input to prepare bulk Si NEP RTA at 300 K, using
> examples/structures/Si/POSCAR_bulk and potentials/Si/Si_Bulk_Fan.txt.
> Write input.yaml and validate it without starting the calculation.

It checks the installed schema, model species and paths, and distinguishes
configuration validity from numerical convergence.
[Skill usage](docs/source/input_assistant.rst)

## Repository layout

| Directory | Contents |
| --- | --- |
| `src/nepkappa/` | Installable source code |
| `examples/` | Reusable input templates and small structures |
| `potentials/` | Models and provenance notes |
| `docs/` | User and developer documentation |
| `tests/` | Local automated checks; ignored by Git |
| `benchmarks/` | Local frozen reference data; ignored by Git |
| `calculations/` | Local research inputs and results; ignored except its guide |
| `.agents/skills/` | Maintained input-assistant skill |

Follow the [development guide](docs/source/development.rst) for local testing
and documentation builds. Tests and reference data are maintained locally and
are not included in a fresh GitHub checkout.

## 中文上手

NEP-kappa 用一个 YAML 输入文件组织声子和晶格热输运计算。
安装后在项目根目录运行：

```bash
nepkappa init input.yaml       # 交互式生成自己的输入文件
nepkappa validate input.yaml   # 检查配置
nepkappa run input.yaml        # 执行所选流程
nepkappa report input.yaml     # 汇总已有结果
```

初次体验可用上面的 Si 算例。计算其他材料时必须替换为对应的结构和势函数。
常规路径相对于命令启动目录；结果统一放在 `calculations/`。
QHA、SSCHA、QHA+SSCHA、三/四声子、Wigner 和 Slurm 的入口见
[示例目录](examples/README.md)，参数含义见 [输入参考](docs/source/input_files.rst)。
`tests`、`benchmarks` 和实际计算结果保留在本地，不随代码提交。

## Citation

F. Yin et al., *Accelerated phonon transport calculations for nanostructures:
Combining neuroevolution potentials and compressed sensing*,
Journal of Applied Physics **139**, 135103 (2026).
[DOI: 10.1063/5.0324012](https://doi.org/10.1063/5.0324012)

Also cite the electronic-structure, potential, and transport methods used in
your calculation; see [references](docs/source/reference.rst).
