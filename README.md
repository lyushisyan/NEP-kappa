# NEP-kappa

[![Docs](https://readthedocs.org/projects/nep-kappa/badge/?version=latest)](https://nep-kappa.readthedocs.io/en/latest/)

[English](#english-version) | [中文](#中文版)


## English Version

NEP-kappa is a workflow for lattice thermal conductivity calculations:

1. Build/relax a structure with NEP
2. Generate force constants (Finite Displacement or HiPhive)
3. Run `phono3py` for thermal conductivity (`kappa`)

### Publication
If you use **NEP-kappa** in your research, please cite the following references.

[1] F. Yin, et al.,
*Accelerated phonon transport calculations for nanostructures: Combining neuroevolution potentials and compressed sensing*,
Journal of Applied Physics **139**, 135103 (2026). https://doi.org/10.1063/5.0324012


### Repository Layout

- `pyproject.toml`: package metadata and install configuration
- `src/nepkappa/`: installable Python package
- `examples-input/input_bulk-rta.yaml`: bulk Si example using finite displacement
- `examples-input/input_film-rta.yaml`: film Si example using HiPhive
- `examples-input/input_bulk-vasp-rta.yaml`: bulk Si example using VASP forces
- `examples-input/`: structures used by packaged YAML inputs
- `tests/`: basic parser and CLI tests
- `Structure/build_and_relax_prim.py`: example script to build and relax Si structure
- `Structure/POSCAR`: output structure file (always overwritten by build script)
- `NEP/Si_Bulk_Fan.txt`: example bulk Si NEP model
- `NEP/Si_NWs_XuKe.txt`: example film/nanowire Si NEP model

### Requirements

Recommended conda env (e.g. `nepkappa-env`) with:

- `numpy`
- `h5py`
- `ase`
- `calorine`
- `hiphive`
- `trainstation`
- `phonopy`
- `phono3py>=4.0.1`
- `seekpath`
- `matplotlib`
- `tqdm`
- `PyYAML`

Install from the repository root:

```bash
python -m pip install -e .
nepkappa --help
```

### 1) Build Structure

Build script supports three Si models:

- `film`: Si(001) slab
- `wire`: 1D Si nanowire
- `bulk`: Si primitive cell (no vacuum; ignores `--thick` and `--vac`)

Examples:

```bash
python Structure/build_and_relax_prim.py --model film --thick 8 --vac 15
python Structure/build_and_relax_prim.py --model wire --thick 8 --vac 15
python Structure/build_and_relax_prim.py --model bulk
```

Notes:

- The script always writes output to `Structure/POSCAR`.
- Example NEP models are provided in `NEP/`.
- This build script is an example implementation for Si models.
- If you have other structures, place them directly in `Structure/` and pass them via `--poscar`.

### 2) Run Workflow

Run the full workflow with recommended YAML input files:

```bash
nepkappa run examples-input/input_bulk-rta.yaml   # Bulk, finite displacement
nepkappa run examples-input/input_film-rta.yaml   # Film, HiPhive
nepkappa run examples-input/input_bulk-vasp-rta.yaml   # Bulk, VASP forces
```

For step-by-step runs, first generate force constants and then compute kappa:

```bash
nepkappa fc examples-input/input_bulk-rta.yaml
nepkappa kappa examples-input/input_bulk-rta.yaml
```

For VASP force calculations, edit `examples-input/input_bulk-vasp-rta.yaml` to match
your VASP executable, POTCAR path, and key INCAR/KPOINTS settings. Common static
single-point defaults are filled by NEP-kappa, then run:

```bash
nepkappa fc examples-input/input_bulk-vasp-rta.yaml
nepkappa kappa examples-input/input_bulk-vasp-rta.yaml
```

You can also run with CLI arguments:

```bash
nepkappa run \
  --poscar examples-input/POSCAR_film \
  --nep_model NEP/Si_NWs_XuKe.txt \
  --do_relax false \
  --dim 4 4 1 \
  --mesh 21 21 1 \
  --temps 100 1000 50 \
  --use_hiphive false \
  --method rta \
  --wigner false \
  --progress true \
  --result_dir result
```

Useful software-style commands:

```bash
nepkappa run examples-input/input_bulk-rta.yaml
nepkappa fc examples-input/input_bulk-rta.yaml
nepkappa kappa examples-input/input_bulk-rta.yaml
nepkappa info examples-input/input_bulk-rta.yaml
```

### Key Arguments

- `--calculator`: `nep` or `vasp`; controls the force backend used by `nepkappa fc`
- `--nep_model`: NEP model path; required when `--calculator nep`
- `--vasp_path`: path to the VASP executable, e.g. `/root/software/vasp.6.4.3/bin/vasp_std`
- `--potcar_path`: path to a POTCAR file or a potential-library directory
- `--vasp_command`: optional full VASP command, e.g. `mpirun -np 64 /path/to/vasp_std`
- `--vasp_workdir`: subdirectory under `result_dir` for per-structure VASP runs
- `--vasp_kwargs`: JSON object with INCAR/KPOINTS options

When `--potcar_path` is a directory, NEP-kappa assembles `POTCAR` by
concatenating element POTCAR files in the element order of the POSCAR.

- `--use_hiphive`: `false` for finite displacement, `true` for HiPhive
- `--temps`: must be either 1 value (`Ts`) or 3 values (`tmin tmax tstep`)
- `--method`: `lbte` or `rta`
- `--wigner`: enable Wigner transport through `phono3py-wte` (`--tt wte`); requires the plugin
- `--progress`: show progress bars, ETA, and timing summaries
- `--result_dir`: directory for generated files and `run.log`

Command behavior:

- `nepkappa fc`: generate `fc2.hdf5`, `fc3.hdf5`, and `phono3py_disp.yaml`
- `nepkappa kappa`: compute thermal conductivity from existing force constants
- `nepkappa run`: run both stages, equivalent to `nepkappa fc` followed by `nepkappa kappa`
- `nepkappa info`: print parsed settings without running

### Typical Outputs

Generated example outputs are local run artifacts. They are written to the
configured `result_dir` such as `examples-output/...` and are not tracked by
the repository.

- `result/run.log`
- `result/POSCAR_relaxed` (if `--do_relax true`)
- `result/phono3py_disp.yaml`
- `result/fc2.hdf5` / `result/fc3.hdf5`
- `result/vasp-runs/` (when `--calculator vasp`)
- `result/kappa-m{mesh}.hdf5`

### Contact info

If you have any questions, please email sxliu98@gmail.com or yinfei0426@outlook.com

---

## 中文版

NEP-kappa 是一个晶格热导率计算流程：

1. 用 NEP 构建/弛豫结构
2. 生成力常数（Finite Displacement 或 HiPhive）
3. 用 `phono3py` 计算热导率 `kappa`

### 论文引用

如果你在科研工作中使用了 **NEP-kappa**，请引用以下文献：

[1] F. Yin, et al.,
*Accelerated phonon transport calculations for nanostructures: Combining neuroevolution potentials and compressed sensing*,
Journal of Applied Physics **139**, 135103 (2026). https://doi.org/10.1063/5.0324012

### 仓库结构

- `pyproject.toml`：软件包元数据和安装配置
- `src/nepkappa/`：可安装的 Python 软件包
- `examples-input/input_bulk-rta.yaml`：bulk Si 有限位移示例输入
- `examples-input/input_film-rta.yaml`：film Si HiPhive 示例输入
- `examples-input/input_bulk-vasp-rta.yaml`：使用 VASP 力的 bulk Si 示例输入
- `examples-input/`：YAML 示例使用的结构文件
- `tests/`：基础解析和命令行测试
- `Structure/build_and_relax_prim.py`：构建并弛豫 Si 结构的示例脚本
- `Structure/POSCAR`：结构输出文件（build 脚本会固定覆盖）
- `NEP/Si_Bulk_Fan.txt`：bulk Si 示例势函数
- `NEP/Si_NWs_XuKe.txt`：film/nanowire Si 示例势函数

### 环境依赖

建议使用 conda 环境（如 `nepkappa-env`），并安装：

- `numpy`
- `h5py`
- `ase`
- `calorine`
- `hiphive`
- `trainstation`
- `phonopy`
- `phono3py>=4.0.1`
- `seekpath`
- `matplotlib`
- `tqdm`
- `PyYAML`

在仓库根目录安装：

```bash
python -m pip install -e .
nepkappa --help
```

### 1) 构建结构

当前支持三种 Si 模型：

- `film`：Si(001) 薄膜
- `wire`：1D Si 纳米线
- `bulk`：Si 原胞（无真空；会忽略 `--thick` 和 `--vac`）

示例：

```bash
python Structure/build_and_relax_prim.py --model film --thick 8 --vac 15
python Structure/build_and_relax_prim.py --model wire --thick 8 --vac 15
python Structure/build_and_relax_prim.py --model bulk
```

说明：

- 结构固定输出到 `Structure/POSCAR`
- 示例势函数位于 `NEP/` 目录
- 这个 build 脚本只是 Si 结构示例实现
- 如果你有其他结构，直接放到 `Structure/` 目录，并通过 `--poscar` 指定

### 2) 运行主流程

使用推荐的 YAML 输入文件运行完整流程：

```bash
nepkappa run examples-input/input_bulk-rta.yaml   # Bulk，有限位移
nepkappa run examples-input/input_film-rta.yaml   # Film，HiPhive
nepkappa run examples-input/input_bulk-vasp-rta.yaml   # Bulk，VASP 力
```

如果想分步运行，可以先生成力常数，再计算热导率：

```bash
nepkappa fc examples-input/input_bulk-rta.yaml
nepkappa kappa examples-input/input_bulk-rta.yaml
```

如果使用 VASP 计算力，请先按你的 VASP 可执行文件、POTCAR 路径和关键
INCAR/KPOINTS 参数修改 `examples-input/input_bulk-vasp-rta.yaml`。常用静态计算参数会由程序默认补齐，然后运行：

```bash
nepkappa fc examples-input/input_bulk-vasp-rta.yaml
nepkappa kappa examples-input/input_bulk-vasp-rta.yaml
```

或直接命令行：

```bash
nepkappa run \
  --poscar examples-input/POSCAR_film \
  --nep_model NEP/Si_NWs_XuKe.txt \
  --do_relax false \
  --dim 4 4 1 \
  --mesh 21 21 1 \
  --temps 100 1000 50 \
  --use_hiphive false \
  --method rta \
  --wigner false \
  --progress true \
  --result_dir result
```

常用软件式命令：

```bash
nepkappa run examples-input/input_bulk-rta.yaml
nepkappa fc examples-input/input_bulk-rta.yaml
nepkappa kappa examples-input/input_bulk-rta.yaml
nepkappa info examples-input/input_bulk-rta.yaml
```

### 关键参数

- `--calculator`：`nep` 或 `vasp`，控制 `nepkappa fc` 使用的力计算后端
- `--nep_model`：NEP 势函数路径；当 `--calculator nep` 时必需
- `--vasp_path`：VASP 可执行文件路径，例如 `/root/software/vasp.6.4.3/bin/vasp_std`
- `--potcar_path`：POTCAR 文件或势函数库目录路径
- `--vasp_command`：可选的完整 VASP 命令，例如 `mpirun -np 64 /path/to/vasp_std`
- `--vasp_workdir`：`result_dir` 下保存各个 VASP 计算目录的子目录
- `--vasp_kwargs`：INCAR/KPOINTS 相关 JSON 参数

当 `--potcar_path` 是目录时，NEP-kappa 会按照 POSCAR 中的元素顺序自动拼接多元素 `POTCAR`。

- `--use_hiphive`：`false` 为有限位移；`true` 为 HiPhive
- `--temps`：只能传 1 个值（单温）或 3 个值（`tmin tmax tstep`）
- `--method`：`lbte` 或 `rta`
- `--wigner`：通过 `phono3py-wte` 插件启用 Wigner 输运（`--tt wte`）；需要额外安装插件
- `--progress`：显示进度条、预计剩余时间和阶段耗时
- `--result_dir`：保存结果文件和 `run.log` 的目录

命令行为：

- `nepkappa fc`：生成 `fc2.hdf5`、`fc3.hdf5` 和 `phono3py_disp.yaml`
- `nepkappa kappa`：使用已有力常数计算热导率
- `nepkappa run`：连续运行两步，等价于先 `nepkappa fc` 再 `nepkappa kappa`
- `nepkappa info`：只打印解析后的配置，不运行计算

### 常见输出

示例计算结果是本地运行产物，会写入配置中的 `result_dir`，例如
`examples-output/...`，仓库不会跟踪这些结果文件。

- `result/run.log`
- `result/POSCAR_relaxed`（当 `--do_relax true`）
- `result/phono3py_disp.yaml`
- `result/fc2.hdf5` / `result/fc3.hdf5`
- `result/vasp-runs/`（当 `--calculator vasp`）
- `result/kappa-m{mesh}.hdf5`

### 联系我们

如果您有任何问题，请联系 sxliu98@gmail.com 或者 yinfei0426@outlook.com

也可加入 QQ 群交流：

![QQ Group](assets/qrcode.png)
