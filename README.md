# NEP-kappa

[![Docs](https://readthedocs.org/projects/nep-kappa/badge/?version=latest)](https://nep-kappa.readthedocs.io/en/latest/)

[English](#english-version) | [中文](#中文版)


## English Version

NEP-kappa is a workflow for lattice thermal conductivity calculations:

1. Build/relax a structure with NEP
2. Generate force constants (Finite Displacement or HiPhive)
3. Run `phono3py` for thermal conductivity (`kappa`)
4. Plot example results for bulk and film

### Publication
If you use **NEP-kappa** in your research, please cite the following references.

[1] F. Yin, et al.,
*Accelerated phonon transport calculations for nanostructures: Combining neuroevolution potentials and compressed sensing*,
Journal of Applied Physics **139**, 135103 (2026). https://doi.org/10.1063/5.0324012


### Repository Layout

- `nepkappa.py`: workflow entry point (CLI args or input file)
- `workflow.py`: core workflow implementation
- `input_1.txt`: example input (Finite Displacement)
- `input_2.txt`: example input (HiPhive)
- `Structure/build_and_relax_prim.py`: example script to build and relax Si structure
- `Structure/POSCAR`: output structure file (always overwritten by build script)
- `NEP/Si_2025_Xuke.txt`: example Si NEP model
- `Example-Results/`: example results and plotting scripts
- `Example-Results/Bulk/plot_bulk.py`: example bulk plotting script (4 panels)
- `Example-Results/Film-1nm/plot_film.py`: example film plotting script (4 panels)

### Requirements

Recommended conda env (e.g. `nepkappa-env`) with:

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
- `phono3py`

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
- Default NEP model is `NEP/Si_2025_Xuke.txt`.
- This build script is an example implementation for Si models.
- If you have other structures, place them directly in `Structure/` and pass them via `--poscar`.

### 2) Run Workflow

Run with input files:

```bash
python nepkappa.py input_1.txt   # Finite Displacement
python nepkappa.py input_2.txt   # HiPhive
```

Or run with CLI arguments:

```bash
python nepkappa.py \
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

### Key Arguments

- `--fc2fc3`: `true` recomputes FC2/FC3; `false` reuses existing FCs
- `--use_hiphive`: `false` for finite displacement, `true` for HiPhive
- `--temps`: must be either 1 value (`Ts`) or 3 values (`tmin tmax tstep`)
- `--method`: `lbte` or `rta`
- `--wigner`: pass `--wigner` to `phono3py`

### Typical Outputs

- `POSCAR_relaxed` (if `--do_relax true`)
- `phono3py_disp.yaml`
- `fc2.hdf5` / `fc3.hdf5`
- `kappa-m{mesh}.hdf5`

### 3) Example Results and Plot

Bulk (4 panels):

```bash
python Example-Results/Bulk/plot_bulk.py
```

Output:

- `Example-Results/Bulk/Comparison_NEP_vs_DFT.pdf`

Film (4 panels):

```bash
python Example-Results/Film-1nm/plot_film.py
```

Output:

- `Example-Results/Film-1nm/Si_film_1nm_4panel.pdf`

Notes:

- These plotting scripts are examples.
- You can copy/modify them for your own structures and result files.

### Contact info

If you have any questions, please email sxliu98@gmail.com or yinfei0426@outlook.com

---

## 中文版

NEP-kappa 是一个晶格热导率计算流程：

1. 用 NEP 构建/弛豫结构
2. 生成力常数（Finite Displacement 或 HiPhive）
3. 用 `phono3py` 计算热导率 `kappa`
4. 绘制 bulk 和 film 示例结果

### 论文引用

如果你在科研工作中使用了 **NEP-kappa**，请引用以下文献：

[1] F. Yin, et al.,
*Accelerated phonon transport calculations for nanostructures: Combining neuroevolution potentials and compressed sensing*,
Journal of Applied Physics **139**, 135103 (2026). https://doi.org/10.1063/5.0324012

### 仓库结构

- `nepkappa.py`：流程入口（支持命令行参数或输入文件）
- `workflow.py`：核心计算逻辑
- `input_1.txt`：Finite Displacement 示例输入
- `input_2.txt`：HiPhive 示例输入
- `Structure/build_and_relax_prim.py`：构建并弛豫 Si 结构的示例脚本
- `Structure/POSCAR`：结构输出文件（build 脚本会固定覆盖）
- `NEP/Si_2025_Xuke.txt`：示例 Si 势函数
- `Example-Results/`：示例结果和绘图代码
- `Example-Results/Bulk/plot_bulk.py`：bulk 绘图示例脚本（4 图）
- `Example-Results/Film-1nm/plot_film.py`：film 绘图示例脚本（4 图）

### 环境依赖

建议使用 conda 环境（如 `nepkappa-env`），并安装：

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
- `phono3py`


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
- 默认势函数为 `NEP/Si_2025_Xuke.txt`
- 这个 build 脚本只是 Si 结构示例实现
- 如果你有其他结构，直接放到 `Structure/` 目录，并通过 `--poscar` 指定

### 2) 运行主流程

使用输入文件：

```bash
python nepkappa.py input_1.txt   # Finite Displacement
python nepkappa.py input_2.txt   # HiPhive
```

或直接命令行：

```bash
python nepkappa.py \
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

### 关键参数

- `--fc2fc3`：`true` 重新计算 FC2/FC3；`false` 复用已有力常数
- `--use_hiphive`：`false` 为有限位移；`true` 为 HiPhive
- `--temps`：只能传 1 个值（单温）或 3 个值（`tmin tmax tstep`）
- `--method`：`lbte` 或 `rta`
- `--wigner`：向 `phono3py` 传递 `--wigner`

### 常见输出

- `POSCAR_relaxed`（当 `--do_relax true`）
- `phono3py_disp.yaml`
- `fc2.hdf5` / `fc3.hdf5`
- `kappa-m{mesh}.hdf5`

### 3) 示例结果与绘图

Bulk（4 图）：

```bash
python Example-Results/Bulk/plot_bulk.py
```

输出：

- `Example-Results/Bulk/Comparison_NEP_vs_DFT.pdf`

Film（4 图）：

```bash
python Example-Results/Film-1nm/plot_film.py
```

输出：

- `Example-Results/Film-1nm/Si_film_1nm_4panel.pdf`

说明：

- 这些绘图脚本是示例脚本
- 你可以按自己的结构和结果文件复制修改后使用

### 联系我们

如果您有任何问题，请联系 sxliu98@gmail.com 或者 yinfei0426@outlook.com

也可加入 QQ 群交流：

![QQ Group](assets/qrcode.png)
