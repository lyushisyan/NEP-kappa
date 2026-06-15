# NEP-kappa

[![Docs](https://readthedocs.org/projects/nep-kappa/badge/?version=latest)](https://nep-kappa.readthedocs.io/en/latest/)

[English](#english-version) | [中文](#中文版)

## English Version

NEP-kappa is an installable workflow package for lattice thermal conductivity
calculations. It can:

1. Relax an input structure.
2. Generate `fc2.hdf5`, `fc3.hdf5`, and `phono3py_disp.yaml`.
3. Compute thermal conductivity with `phono3py`.
4. Plot phonon and thermal-transport results from existing HDF5 outputs.

Available commands:

```bash
nepkappa relax input.yaml
nepkappa fc input.yaml
nepkappa kappa input.yaml
nepkappa plot input.yaml
nepkappa run input.yaml
nepkappa info input.yaml
```

- `nepkappa relax`: relax the input structure
- `nepkappa fc`: generate `fc2.hdf5`, `fc3.hdf5`, and `phono3py_disp.yaml`
- `nepkappa kappa`: compute thermal conductivity from existing force constants
- `nepkappa plot`: plot dispersion, DOS, volume heat capacity, group velocity, relaxation time, and thermal conductivity
- `nepkappa run`: run `relax`, `fc`, and `kappa` in sequence
- `nepkappa info`: print the parsed configuration without running

### Publication

If you use **NEP-kappa** in your research, please cite:

[1] F. Yin, et al.,
*Accelerated phonon transport calculations for nanostructures: Combining neuroevolution potentials and compressed sensing*,
Journal of Applied Physics **139**, 135103 (2026). https://doi.org/10.1063/5.0324012

### Repository Layout

- `pyproject.toml`: package metadata and install configuration
- `src/nepkappa/`: installable Python package
- `examples/`: packaged POSCAR and YAML examples
- `potentials/`: packaged Si NEP model files
- `tests/`: parser and CLI smoke tests
- `docs/`: Sphinx documentation

Generated outputs are written to the configured `result_dir`, commonly
`results/...`, and are not tracked by the repository.

### Requirements

Python dependencies:

- `python>=3.9`
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

If an older `phono3py` is already installed:

```bash
python -m pip install --upgrade -e .
```

### Example Inputs

The packaged examples are:

- `examples/1-bulk-nep-rta.yaml`: bulk, NEP forces, finite displacement, RTA
- `examples/2-bulk-nep-hiphive-rta.yaml`: bulk, NEP forces, HiPhive, RTA
- `examples/3-bulk-nep-lbte.yaml`: bulk, NEP forces, finite displacement, LBTE
- `examples/4-bulk-nep-rta-wigner.yaml`: bulk, NEP forces, finite displacement, Wigner transport
- `examples/5-bulk-vasp-rta.yaml`: bulk, VASP relaxation and VASP forces, finite displacement, RTA
- `examples/6-bulk-vasp-hiphive-rta.yaml`: bulk, VASP relaxation and VASP forces, HiPhive, RTA
- `examples/7-film-nep-hiphive-rta.yaml`: film, NEP forces, HiPhive, RTA

Use `info` before running an expensive calculation:

```bash
nepkappa info examples/1-bulk-nep-rta.yaml
```

Run the complete workflow:

```bash
nepkappa run examples/1-bulk-nep-rta.yaml
```

Or run the stages separately:

```bash
nepkappa relax examples/1-bulk-nep-rta.yaml
nepkappa fc examples/1-bulk-nep-rta.yaml
nepkappa kappa examples/1-bulk-nep-rta.yaml
nepkappa plot examples/1-bulk-nep-rta.yaml
```

When `relaxation.enabled: true`, `nepkappa fc` expects
`POSCAR_relaxed` to already exist in `result_dir`. Use `nepkappa relax` first,
or use `nepkappa run`.

### YAML Structure

The recommended YAML structure is:

```yaml
structure:
  poscar: examples/POSCAR_bulk

calculator:
  name: nep
  nep_model: potentials/Si_Bulk_Fan.txt

relaxation:
  enabled: true

force-constant:
  dim: [3, 3, 3]
  fc_calculator: traditional
  use_hiphive: false

kappa:
  mesh: [21, 21, 21]
  temps: [100, 1000, 50]
  method: rta
  wigner: false

plot:
  layout: separate
  path: seekpath
  tau: total
  kappa: all
  temperature: 300
  dpi: 300

output:
  progress: true
  result_dir: results/1-bulk-nep-rta
```

For VASP calculations, use `calculator.name: vasp` and set
`vasp_command` or `vasp_path`, plus `potcar_path`. If `potcar_path` is a
directory, NEP-kappa assembles a combined `POTCAR` in POSCAR element order.

### Key Options

- `calculator.name`: `nep` or `vasp`
- `calculator.nep_model`: NEP model path; required for `nep`
- `calculator.vasp_command`: full VASP command, e.g. with `mpirun`
- `calculator.vasp_path`: VASP executable path
- `calculator.potcar_path`: POTCAR file or potential-library directory
- `relaxation.enabled`: whether to run structure relaxation
- `force-constant.dim`: supercell dimension for force constants
- `force-constant.use_hiphive`: `false` for finite displacement, `true` for HiPhive
- `force-constant.fc_calculator`: finite-displacement solver, e.g. `symfc`
- `kappa.method`: `rta` or `lbte`
- `kappa.wigner`: use `phono3py-wte` via `--tt wte`
- `plot.layout`: `separate`, `combined`, or `both`
- `plot.path`: high-symmetry path source, `seekpath` or `custom`
- `plot.tau`: relaxation-time channel, `total`, `normal`, `umklapp`, or `all`
- `plot.kappa`: thermal-conductivity component, `x`, `y`, `z`, or `all`
- `plot.temperature`: target temperature for relaxation-time plots
- `plot.dpi`: output figure resolution
- `output.result_dir`: directory for generated files and `run.log`

For a custom dispersion path, define fractional reciprocal coordinates:

```yaml
plot:
  path: custom
  path_points:
    G: [0.0, 0.0, 0.0]
    X: [0.5, 0.0, 0.5]
    U: [0.625, 0.25, 0.625]
    K: [0.375, 0.375, 0.75]
    L: [0.5, 0.5, 0.5]
    W: [0.5, 0.25, 0.75]
  path_segments:
    - [G, X]
    - [X, U]
    - [K, G]
    - [G, L]
    - [L, W]
    - [W, X]
```

Disconnected neighboring segments are shown with a combined label, e.g.
`[X, U]` followed by `[K, G]` is plotted as `U|K`.

### Plotting

`nepkappa plot input.yaml` reads `phono3py_disp.yaml`, `fc2.hdf5`, and
`kappa-m*.hdf5` from `output.result_dir`. It always generates six standard
figures: phonon dispersion, DOS, volume heat capacity, group velocity,
relaxation time, and thermal conductivity.

The `plot` section controls figure layout and selected data channels. Use
`layout: separate` for individual PNG files, `layout: combined` for one 2-by-3
summary figure, or `layout: both` for both outputs. Group velocity is plotted
in km/s. The figures use large axis and tick labels for publication use and do
not add subplot titles.

### Typical Outputs

- `run.log`
- `POSCAR_relaxed`
- `phono3py_disp.yaml`
- `fc2.hdf5`
- `fc3.hdf5`
- `hiphive_model.fcp` when HiPhive is used
- `vasp-relax/` and `vasp-runs/` when VASP is used
- `kappa-m{mesh}.hdf5`
- `phono3py.yaml`
- `plots/dispersion.png`
- `plots/dos.png`
- `plots/heat_capacity.png`
- `plots/group_velocity.png` with group velocity in km/s
- `plots/relaxation_time.png`
- `plots/kappa.png`
- `plots/combined.png` when `plot.layout` is `combined` or `both`

### Contact

For questions, please email sxliu98@gmail.com or yinfei0426@outlook.com.

---

## 中文版

NEP-kappa 是一个可安装的软件包，用于晶格热导率计算。它可以：

1. 弛豫输入结构。
2. 生成 `fc2.hdf5`、`fc3.hdf5` 和 `phono3py_disp.yaml`。
3. 调用 `phono3py` 计算热导率。
4. 基于已有 HDF5 结果绘制声子和热输运图像。

可用命令为：

```bash
nepkappa relax input.yaml
nepkappa fc input.yaml
nepkappa kappa input.yaml
nepkappa plot input.yaml
nepkappa run input.yaml
nepkappa info input.yaml
```

- `nepkappa relax`：弛豫输入结构
- `nepkappa fc`：生成 `fc2.hdf5`、`fc3.hdf5` 和 `phono3py_disp.yaml`
- `nepkappa kappa`：使用已有力常数计算热导率
- `nepkappa plot`：绘制色散关系、态密度、体积热容、群速度、弛豫时间和热导率
- `nepkappa run`：连续执行 `relax`、`fc` 和 `kappa`
- `nepkappa info`：只打印解析后的配置，不运行计算

### 论文引用

如果你在科研工作中使用了 **NEP-kappa**，请引用：

[1] F. Yin, et al.,
*Accelerated phonon transport calculations for nanostructures: Combining neuroevolution potentials and compressed sensing*,
Journal of Applied Physics **139**, 135103 (2026). https://doi.org/10.1063/5.0324012

### 仓库结构

- `pyproject.toml`：软件包元数据和安装配置
- `src/nepkappa/`：可安装的 Python 软件包
- `examples/`：打包的 POSCAR 和 YAML 示例
- `potentials/`：打包的 Si NEP 势函数
- `tests/`：解析器和 CLI 基础测试
- `docs/`：Sphinx 文档

计算输出会写入 YAML 中的 `result_dir`，通常为 `results/...`，
仓库不会跟踪这些本地运行产物。

Python 依赖：

- `python>=3.9`
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

如果环境中已有旧版 `phono3py`：

```bash
python -m pip install --upgrade -e .
```

### 示例输入

当前提供的 YAML 示例为：

- `examples/1-bulk-nep-rta.yaml`：bulk，NEP 力，有限位移，RTA
- `examples/2-bulk-nep-hiphive-rta.yaml`：bulk，NEP 力，HiPhive，RTA
- `examples/3-bulk-nep-lbte.yaml`：bulk，NEP 力，有限位移，LBTE
- `examples/4-bulk-nep-rta-wigner.yaml`：bulk，NEP 力，有限位移，Wigner 输运
- `examples/5-bulk-vasp-rta.yaml`：bulk，VASP 弛豫和 VASP 力，有限位移，RTA
- `examples/6-bulk-vasp-hiphive-rta.yaml`：bulk，VASP 弛豫和 VASP 力，HiPhive，RTA
- `examples/7-film-nep-hiphive-rta.yaml`：film，NEP 力，HiPhive，RTA

正式运行前建议先检查配置：

```bash
nepkappa info examples/1-bulk-nep-rta.yaml
```

运行完整流程：

```bash
nepkappa run examples/1-bulk-nep-rta.yaml
```

也可以分步运行：

```bash
nepkappa relax examples/1-bulk-nep-rta.yaml
nepkappa fc examples/1-bulk-nep-rta.yaml
nepkappa kappa examples/1-bulk-nep-rta.yaml
nepkappa plot examples/1-bulk-nep-rta.yaml
```

当 `relaxation.enabled: true` 时，单独运行 `nepkappa fc` 会要求
`result_dir` 中已经存在 `POSCAR_relaxed`。这种情况下请先运行
`nepkappa relax`，或者直接使用 `nepkappa run`。

### YAML 结构

推荐 YAML 结构为：

```yaml
structure:
  poscar: examples/POSCAR_bulk

calculator:
  name: nep
  nep_model: potentials/Si_Bulk_Fan.txt

relaxation:
  enabled: true

force-constant:
  dim: [3, 3, 3]
  fc_calculator: traditional
  use_hiphive: false

kappa:
  mesh: [21, 21, 21]
  temps: [100, 1000, 50]
  method: rta
  wigner: false

plot:
  layout: separate
  path: seekpath
  tau: total
  kappa: all
  temperature: 300
  dpi: 300

output:
  progress: true
  result_dir: results/1-bulk-nep-rta
```

VASP 计算使用 `calculator.name: vasp`，并设置 `vasp_command` 或
`vasp_path`，以及 `potcar_path`。如果 `potcar_path` 是目录，NEP-kappa
会按照 POSCAR 元素顺序自动拼接 `POTCAR`。

### 关键选项

- `calculator.name`：`nep` 或 `vasp`
- `calculator.nep_model`：NEP 势函数路径；使用 `nep` 时必需
- `calculator.vasp_command`：完整 VASP 命令，例如带 `mpirun`
- `calculator.vasp_path`：VASP 可执行文件路径
- `calculator.potcar_path`：POTCAR 文件或势函数库目录
- `relaxation.enabled`：是否进行结构弛豫
- `force-constant.dim`：力常数超胞尺寸
- `force-constant.use_hiphive`：`false` 为有限位移，`true` 为 HiPhive
- `force-constant.fc_calculator`：有限位移求解器，例如 `symfc`
- `kappa.method`：`rta` 或 `lbte`
- `kappa.wigner`：通过 `phono3py-wte` 使用 `--tt wte`
- `plot.layout`：`separate`、`combined` 或 `both`
- `plot.path`：高对称路径来源，`seekpath` 或 `custom`
- `plot.tau`：弛豫时间通道，`total`、`normal`、`umklapp` 或 `all`
- `plot.kappa`：热导率方向，`x`、`y`、`z` 或 `all`
- `plot.temperature`：弛豫时间图使用的目标温度
- `plot.dpi`：输出图片分辨率
- `output.result_dir`：结果文件和 `run.log` 的输出目录

如果需要自定义色散路径，可以写分数倒空间坐标：

```yaml
plot:
  path: custom
  path_points:
    G: [0.0, 0.0, 0.0]
    X: [0.5, 0.0, 0.5]
    U: [0.625, 0.25, 0.625]
    K: [0.375, 0.375, 0.75]
    L: [0.5, 0.5, 0.5]
    W: [0.5, 0.25, 0.75]
  path_segments:
    - [G, X]
    - [X, U]
    - [K, G]
    - [G, L]
    - [L, W]
    - [W, X]
```

相邻两段如果不连续，会在横坐标断点处合并显示，例如 `[X, U]`
后接 `[K, G]` 会显示为 `U|K`。

### 绘图功能

`nepkappa plot input.yaml` 会从 `output.result_dir` 读取
`phono3py_disp.yaml`、`fc2.hdf5` 和 `kappa-m*.hdf5`。它固定生成六类标准图：
声子色散、态密度、体积热容、群速度、弛豫时间和热导率。

`plot` 部分用于控制图像布局和数据通道。`layout: separate` 会分别输出 PNG，
`layout: combined` 会输出一个 2 行 3 列的总图，`layout: both` 会同时输出两种。
群速度单位为 km/s。图像采用较大的坐标轴和刻度字号，适合论文排版，并且不添加子图小标题。

### 常见输出

- `run.log`
- `POSCAR_relaxed`
- `phono3py_disp.yaml`
- `fc2.hdf5`
- `fc3.hdf5`
- 使用 HiPhive 时的 `hiphive_model.fcp`
- 使用 VASP 时的 `vasp-relax/` 和 `vasp-runs/`
- `kappa-m{mesh}.hdf5`
- `phono3py.yaml`
- `plots/dispersion.png`
- `plots/dos.png`
- `plots/heat_capacity.png`
- `plots/group_velocity.png`，群速度单位为 km/s
- `plots/relaxation_time.png`
- `plots/kappa.png`
- `plots/combined.png`，当 `plot.layout` 为 `combined` 或 `both` 时生成

### 联系我们

如果您有任何问题，请联系 sxliu98@gmail.com 或者 yinfei0426@outlook.com。

也可加入 QQ 群交流：

![QQ Group](assets/qrcode.png)
