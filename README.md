# NEP-kappa

[English](#english-version) | [中文](#中文版)

## English Version


NEP-kappa is an automated workflow for lattice thermal conductivity calculations:

1. Compute forces with NEP potential (`calorine`)
2. Build force constants via Finite Displacement or HiPhive
3. Run `phono3py` to compute `kappa`
4. Plot NEP results with `plot.py`

### Repository Layout

- `main.py`: workflow entry point (CLI args or input file)
- `workflow.py`: core workflow implementation
- `input_1.txt`: example input for Finite Displacement
- `input_2.txt`: example input for HiPhive
- `Structure/build_and_relax_prim.py`: build and relax structures, output `POSCAR`
- `Structure/POSCAR`: example structure
- `NEP/Si_2025_Xuke.txt`: example Si potential
- `plot.py`: plot NEP dispersion / DOS / group velocity / relaxation time

### Requirements

Use a conda environment (e.g. `ph3-env`) with:

- `numpy`
- `h5py`
- `ase`
- `calorine`
- `hiphive`
- `trainstation`
- `phonopy`
- `phono3py`
- `seekpath`
- `matplotlib`

Also ensure `phono3py` command is available in your shell.

### Quick Start

1) Prepare structure and potential

- Structure: `--poscar` (e.g. `Structure/POSCAR`)
- Potential: `--nep_model` (e.g. `NEP/Si_2025_Xuke.txt`)

To generate a structure automatically:

```bash
python Structure/build_and_relax_prim.py --model film --thick 8 --vac 15
```

2) Run the workflow

```bash
python main.py input_1.txt   # Finite Displacement
python main.py input_2.txt   # HiPhive
```

Or use CLI arguments directly:

```bash
python main.py \
  --poscar Structure/POSCAR \
  --nep_model NEP/Si_2025_Xuke.txt \
  --do_relax false \
  --dim 4 4 1 \
  --mesh 21 21 1 \
  --temps 100 1000 50 \
  --fc2fc3 true \
  --use_hiphive false \
  --method rta \
  --wigner true
```

### Two Force-Constant Paths

- `--use_hiphive false`: Finite Displacement
- `--use_hiphive true`: HiPhive

### Key Arguments

- `--fc2fc3`: `true` to recompute FC2/FC3, `false` to reuse existing FCs
- `--temps`: accepts only 1 value (single temperature) or 3 values (`tmin tmax tstep`)
- `--method`: `lbte` or `rta`
- `--wigner`: enable `phono3py --wigner`

### Outputs

Typical outputs:

- `POSCAR_relaxed` (if `--do_relax true`)
- `phono3py_disp.yaml`
- `fc2.hdf5` / `fc3.hdf5`
- `FORCE_CONSTANTS` / `FORCE_CONSTANTS_3RD`
- `kappa-m{mesh}.hdf5`

### Plotting

```bash
python plot.py \
  --poscar Structure/POSCAR \
  --model NEP/Si_2025_Xuke.txt \
  --h5 kappa-m21211.hdf5 \
  --supercell 4 4 1 \
  --output NEP_results.pdf
```

Optional arguments:

- `--dos-mesh`
- `--freq-max`
- `--temp`
- `--truncate-label`

Show all options:

```bash
python plot.py --help
```

If not provided, defaults are:

- `NEP/POSCAR`
- `NEP/Si_2025_Xuke.txt`
- `NEP/kappa-m212121.hdf5`

### Troubleshooting

- `Error: Could not import 'workflow.py'`
  - Usually missing dependencies (`phono3py`, `hiphive`, `trainstation`)

- `phono3py failed with return code ...`
  - Check whether `fc2.hdf5` / `fc3.hdf5` are generated
  - Check whether `--dim` / `--mesh` are too large

- `--temps` validation error
  - Only 1 or 3 numeric values are supported

---

## 中文版


NEP-kappa 是一个用于晶格热导率计算的自动化流程，核心链路为：

1. 使用 NEP 势（`calorine`）计算受力
2. 通过 Finite Displacement 或 HiPhive 生成力常数
3. 调用 `phono3py` 计算热导率 `kappa`
4. 使用 `plot.py` 绘制 NEP 结果

### 目录说明

- `main.py`：流程入口（支持命令行参数或输入文件）
- `workflow.py`：核心计算逻辑
- `input_1.txt`：Finite Displacement 示例输入
- `input_2.txt`：HiPhive 示例输入
- `Structure/build_and_relax_prim.py`：生成并弛豫结构，输出 `POSCAR`
- `Structure/POSCAR`：示例结构
- `NEP/Si_2025_Xuke.txt`：示例 Si 势函数
- `plot.py`：绘制 NEP 的色散 / DOS / 群速度 / 弛豫时间

### 环境依赖

建议使用 `conda` 环境（如 `ph3-env`），并安装：

- `numpy`
- `h5py`
- `ase`
- `calorine`
- `hiphive`
- `trainstation`
- `phonopy`
- `phono3py`
- `seekpath`
- `matplotlib`

并确保终端可直接调用 `phono3py` 命令。

### 快速开始

1) 准备结构和势函数

- 结构：`--poscar`（例如 `Structure/POSCAR`）
- 势函数：`--nep_model`（例如 `NEP/Si_2025_Xuke.txt`）

如需自动构建结构：

```bash
python Structure/build_and_relax_prim.py --model film --thick 8 --vac 15
```

2) 运行主流程

```bash
python main.py input_1.txt   # Finite Displacement
python main.py input_2.txt   # HiPhive
```

或命令行参数方式：

```bash
python main.py \
  --poscar Structure/POSCAR \
  --nep_model NEP/Si_2025_Xuke.txt \
  --do_relax false \
  --dim 4 4 1 \
  --mesh 21 21 1 \
  --temps 100 1000 50 \
  --fc2fc3 true \
  --use_hiphive false \
  --method rta \
  --wigner true
```

### 两种力常数路径

- `--use_hiphive false`：Finite Displacement
- `--use_hiphive true`：HiPhive

### 关键参数

- `--fc2fc3`：`true` 重新计算 FC2/FC3；`false` 使用已有力常数
- `--temps`：只允许 1 个值（单温）或 3 个值（`tmin tmax tstep`）
- `--method`：`lbte` 或 `rta`
- `--wigner`：是否开启 `phono3py --wigner`

### 输出文件

常见输出：

- `POSCAR_relaxed`（当 `--do_relax true`）
- `phono3py_disp.yaml`
- `fc2.hdf5` / `fc3.hdf5`
- `FORCE_CONSTANTS` / `FORCE_CONSTANTS_3RD`
- `kappa-m{mesh}.hdf5`

### 绘图

```bash
python plot.py \
  --poscar Structure/POSCAR \
  --model NEP/Si_2025_Xuke.txt \
  --h5 kappa-m21211.hdf5 \
  --supercell 4 4 1 \
  --output NEP_results.pdf
```

可选参数：

- `--dos-mesh`
- `--freq-max`
- `--temp`
- `--truncate-label`

查看完整参数：

```bash
python plot.py --help
```

不传参数时默认使用：

- `NEP/POSCAR`
- `NEP/Si_2025_Xuke.txt`
- `NEP/kappa-m212121.hdf5`

### 常见问题

- `Error: Could not import 'workflow.py'`
  - 通常是依赖未安装（常见：`phono3py`、`hiphive`、`trainstation`）

- `phono3py failed with return code ...`
  - 检查 `fc2.hdf5` / `fc3.hdf5` 是否生成
  - 检查 `--dim` / `--mesh` 是否过大

- `--temps` 报错
  - 仅支持 1 个或 3 个数值输入
