
输入文件与参数
========================



``structure``：晶体结构
----------------------------------------

``structure`` 用于指定晶体结构、结构弛豫及低维体系的几何归一化参数。

.. list-table:: 结构参数
   :header-rows: 1
   :widths: 38 62

   * - 参数
     - 说明
   * - ``poscar``
     - 晶体结构文件，采用 VASP POSCAR 格式。
   * - ``relaxation``
     - 是否进行结构弛豫，支持布尔值或详细配置。
   * - ``dimensionality``
     - 体系维度：``1``、``2`` 或 ``3``，默认 ``3``。
   * - ``effective-thickness``
     - 二维体系的有效厚度，单位 Å。
   * - ``effective-area``
     - 一维体系的有效横截面积，单位 Å²。
   * - ``vacuum-axis``
     - 二维体系的真空方向，可选。
   * - ``periodic-axis``
     - 一维体系的周期方向，可选。


结构弛豫
~~~~~~~~~~~~~~~~~~~~~~~~

设置 ``relaxation: true`` 启用结构弛豫。

对于 VASP，还可以设置多阶段弛豫：

.. code-block:: yaml

   structure:
     poscar: POSCAR
     relaxation:
       enabled: true
       workdir: vasp-relax
       stages:
         coarse: {ediffg: -0.05, prec: Normal}
         fine: {ediffg: -0.01, prec: Accurate}

.. note::

  ``nepkappa run`` 会自动执行结构弛豫及后续计算。
  使用 ``fc2``、``fc2fc3`` 或 ``fc2fc3fc4`` 独立生成力常数时，
  如果启用了结构弛豫，需要先运行 ``nepkappa relax``，
  生成 ``POSCAR_relaxed``。
  QHA 使用指定的输入结构，要求 ``relaxation: false``。


低维体系的几何归一化
~~~~~~~~~~~~~~~~~~~~~~~~

对于包含真空层的二维薄膜和一维纳米线，
需要使用有效厚度或横截面积进行物理体积归一化。

.. list-table:: 低维体系设置
   :header-rows: 1
   :widths: 20 72 38

   * - 体系
     - 必需参数
     - 可选参数
   * - 二维薄膜
     - ``dimensionality: 2`` 和
       ``effective-thickness``
     - ``vacuum-axis``
   * - 一维纳米线
     - ``dimensionality: 1`` 和
       ``effective-area``
     - ``periodic-axis``

上述参数用于绘图时的体积热容和热导率归一化，
不修改原始 phono3py HDF5 数据。

.. note::

   QHA、SSCHA 和 FourPhonon 目前仅支持三维体系
   （``dimensionality: 3``）。
   几何归一化参数不会改变这些方法的适用范围。



``calculator``：原子间相互作用
----------------------------------------

``calculator`` 指定能量和原子力的计算方法。
NEP-kappa 支持 NEP、VASP、MACE 和其他 ASE 兼容计算器。

.. list-table:: 计算器类型
   :header-rows: 1
   :widths: 16 16 64

   * - 计算器
     - ``name``
     - 主要配置
   * - NEP
     - ``nep``
     - 通过 ``model`` 指定 NEP 势文件，默认使用 CPU 计算。
   * - MACE
     - ``mace``
     - 通过 ``model`` 或 ``foundation`` 指定模型，可使用 ``device`` 和 ``dtype`` 控制计算设备与精度。
   * - VASP
     - ``vasp``
     - 配置 ``vasp-command``、``vasp-path`` 和 ``potcar-path``。
   * - ASE
     - ``ase``
     - 使用 ``factory``、``kwargs`` 和 ``model-files`` 配置其他计算器。


NEP GPU 加速
~~~~~~~~~~~~~~~~~~~~~~~~

对于大量位移结构，可以使用 GPU 批量计算原子力：

.. code-block:: yaml

   calculator:
     name: nep
     model: /path/to/nep.txt
     nep-force-backend: gpu-batch
     nep-batch-size: 128
     nep-command: /path/to/nep

主要参数：

- ``nep-force-backend``：选择 ``cpu`` 或 ``gpu-batch``。
- ``nep-batch-size``：每批处理的结构数量上限。
- ``nep-command``：GPUMD ``nep`` 可执行文件的路径。

GPU 批处理适用于原生 phono3py FC2/FC3 和 HiPhive 位移力计算，目前不支持 QHA、SSCHA、Thirdorder 和 Fourthorder。

.. note::

   GPU 批处理需要 Calorine 3.5 及以上版本、
   GPUMD 的 ``nep`` 可执行文件和兼容的 GPU 环境。


``force-constant``：力常数
----------------------------------------

``force-constant`` 控制二阶、三阶和四阶力常数的生成。
可使用有限位移方法或 HiPhive 拟合，并支持不同的力常数计算后端。

.. list-table:: 力常数参数
   :header-rows: 1
   :widths: 35 65

   * - 参数
     - 说明
   * - ``dim-fc2`` / ``dim-fc3`` / ``dim-fc4``
     - 分别指定 FC2、FC3 和 FC4 的超胞维度。
   * - ``fc3-backend``
     - 默认使用 phono3py，也可选择 ``thirdorder``。
   * - ``cutoff-fc3``
     - 三阶力常数的位移截断设置，具体定义取决于计算后端。
   * - ``cutoff-fc4``
     - Fourthorder 的四阶相互作用截断设置。
   * - ``use-hiphive``
     - 启用 HiPhive 力常数拟合。
   * - ``format``
     - 输出格式：``phono3py``、``shengbte`` 或 ``both``。
   * - ``thirdorder-command`` / ``fourthorder-command``
     - Thirdorder 和 Fourthorder 可执行程序的路径或启动命令。
   * - ``qha`` / ``sscha`` / ``four-phonon``
     - 启用相应的计算模块，并设置其参数。

力常数超胞和截断
~~~~~~~~~~~~~~~~~~~~~~~~

``dim-fc2``、``dim-fc3`` 和 ``dim-fc4`` 分别控制不同阶次的超胞。
超胞大小会影响计算成本及力常数的收敛性，应针对不同材料进行测试。

``cutoff-fc3`` 的定义取决于 FC3 后端：

- **phono3py**：位移原子对的距离截断，单位 Å。
- **Thirdorder**：负整数表示邻居壳层，正值表示距离截断，单位 nm。

``cutoff-fc4`` 传递给 Fourthorder，负整数表示邻居壳层。

HiPhive 使用 ``n-structures``、``rattle-std``、``cutoffs`` 和 ``min-dist``
等参数控制训练结构与力常数拟合，其中 ``cutoffs`` 的单位为 Å。


QHA 与 SSCHA
~~~~~~~~~~~~~~~~~~~~~~~~

NEP-kappa 支持准谐近似（QHA）、固定体积随机自洽谐近似（SSCHA），
以及 QHA 平衡体积上的 SSCHA 计算。

三种路线的主要区别如下：

.. list-table:: 温度相关声子计算
   :header-rows: 1
   :widths: 25 75

   * - 方法
     - 计算内容
   * - QHA
     - 通过体积扫描和状态方程拟合获得温度相关平衡体积，描述热膨胀对声子性质的影响。
   * - SSCHA
     - 在固定体积下进行有限温度采样，自洽更新辅助谐近似力常数。
   * - QHA + SSCHA
     - 先由 QHA 确定各温度的平衡体积，再在对应体积上进行固定体积 SSCHA。

主要参数：

- ``qha.volume-ratios``：体积扫描范围。
- ``qha.temps``：QHA 的温度范围。
- ``sscha.temps``：SSCHA 采样温度。
- ``sscha.snapshots``：每次迭代的采样结构数量。
- ``sscha.iterations``：自洽迭代次数。
- ``run-transport``：是否在相应温度计算热导率。

详细的计算方法与完整示例见 :doc:`temperature`。



``kappa``：声子散射与热导率
----------------------------------------

``kappa`` 用于配置声子散射、玻尔兹曼输运方程（BTE）求解方法及晶格热导率计算。
NEP-kappa 支持 phono3py 和 FourPhonon 两种计算引擎，可选择不同的散射机制与求解方法。


常用参数
~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: 热输运参数
   :header-rows: 1
   :widths: 36 64

   * - 参数
     - 说明
   * - ``engine``
     - 计算引擎：``phono3py`` 或 ``fourphonon``。
   * - ``mesh``
     - 布里渊区 q 网格，例如 ``[31, 31, 31]``。
   * - ``temps``
     - 计算温度，单位 K。
   * - ``method``
     - phono3py 求解方法：``rta`` 或 ``lbte``。
   * - ``method-3ph`` / ``method-4ph``
     - FourPhonon 中三声子与四声子的求解方法。
   * - ``isotope``
     - 是否考虑同位素散射。
   * - ``bfmp``
     - 边界散射的平均自由程参数，单位 μm。
   * - ``wigner``
     - 是否启用 Wigner 相干输运计算。
   * - ``force-constants``
     - 已有力常数的来源目录，仅适用于普通 phono3py 单阶段输运。
   * - ``command``
     - 自定义 phono3py 启动命令。


计算引擎与求解方法
~~~~~~~~~~~~~~~~~~~~~~~~

phono3py 主要用于三声子热输运，FourPhonon 则同时考虑三声子和四声子散射。
两种引擎均支持 RTA 和迭代 LBTE，但可用的求解组合不同。

.. list-table:: 支持的热输运方法
   :header-rows: 1
   :widths: 25 30 45

   * - 计算引擎
     - 求解方法
     - 计算内容
   * - phono3py
     - ``rta``
     - 三声子弛豫时间近似。
   * - phono3py
     - ``lbte``
     - 三声子线性化 BTE 迭代求解。
   * - FourPhonon
     - ``rta / rta``
     - 三声子与四声子均采用 RTA。
   * - FourPhonon
     - ``lbte / rta``
     - 三声子迭代求解，四声子采用 RTA。
   * - FourPhonon
     - ``lbte / lbte``
     - 三声子与四声子均采用迭代求解。

RTA 使用各声子模式的弛豫时间近似处理碰撞项；
迭代 LBTE 则进一步考虑不同声子模式之间的线性耦合。

FourPhonon 使用 ``method-3ph`` 和 ``method-4ph`` 分别指定求解方法。
目前不支持 ``rta / lbte`` 组合。

Wigner 相干输运
~~~~~~~~~~~~~~~~~~~~~~~~

设置 ``wigner: true`` 可启用 Wigner 输运计算，
将晶格热导率分为粒子贡献（particle）和相干贡献（coherence）。

总热导率为：

.. math::

   \kappa_{\mathrm{total}}
   = \kappa_{\mathrm{P}} + \kappa_{\mathrm{C}}.

其中，粒子项描述声子布居的输运贡献，
相干项描述不同振动模式之间的非对角相干贡献。

Wigner 并非独立的 BTE 求解方法，其支持范围取决于计算引擎：

- **phono3py**：支持三声子 RTA 或 LBTE，并通过 SMM19 计算相干贡献。
- **FourPhonon**：需要 Wigner_Park，目前仅支持三声子和四声子均采用 RTA。

具体输入和计算示例见 :doc:`thermal_transport`。


计算温度
~~~~~~~~~~~~~~~~~~~~~~~~

``temps`` 支持单温度和等间隔温度范围两种格式：

.. list-table:: 温度设置
   :header-rows: 1
   :widths: 35 65

   * - 格式
     - 示例
   * - 单温度 ``[T]``
     - ``temps: [300]``
   * - 温度范围 ``[Tmin, Tmax, step]``
     - ``temps: [100, 1000, 50]``

温度单位为 K，必须大于零。范围格式包含起点和终点，
步长必须为正，且能够恰好到达终点。

目前不支持直接输入任意不等间隔的温度列表。

对于 QHA 或 SSCHA 耦合输运，热导率温度还必须满足对应的采样范围要求，
详见 :doc:`temperature`。


复用已有力常数
~~~~~~~~~~~~~~~~~~~~~~~~

如果已经生成 FC2 和 FC3，可以跳过原子力计算，
直接使用 phono3py 计算不同温度或 q 网格下的热导率。

例如：

.. code-block:: yaml

   kappa:
     force-constants: calculations/3C-SiC/fc
     engine: phono3py
     mesh: [41, 41, 41]
     temps: [300]
     method: rta

   output:
     result-dir: calculations/3C-SiC/kappa-m41

首先检查输入，再执行计算：

.. code-block:: bash

   nepkappa info kappa.yaml --for kappa
   nepkappa kappa kappa.yaml

``force-constants`` 指定已有力常数目录，
``output.result-dir`` 指定新的计算结果目录。

来源目录必须包含以下配套文件：

- ``fc2.hdf5``：二阶力常数。
- ``fc3.hdf5``：三阶力常数。
- ``phono3py_disp.yaml``：晶体结构、位移及相关元数据。

当来源目录与结果目录不同时，
程序会将所需输入文件复制到结果目录，不修改原始文件。


自定义 phono3py 命令
~~~~~~~~~~~~~~~~~~~~~~~~

高级用户可以通过 ``kappa.command`` 指定完整的 phono3py 启动命令，
代替 NEP-kappa 自动生成的命令行。

自定义命令中的网格和温度必须与 YAML 配置一致。启用后，NEP-kappa 不再自动附加其他计算参数。

命令直接在结果目录中执行，不经过 shell，因此不支持管道或重定向操作。

更多热输运计算方法与完整示例见 :doc:`thermal_transport`。


``tdbte``：非平衡声子动力学
----------------------------------------

``tdbte`` 用于求解含时玻尔兹曼输运方程，
模拟非平衡声子占据数随时间的演化和模式间能量交换。

与热导率计算不同，TD-BTE 不需要预先求解稳态输运方程，
但需要匹配的 FC2、FC3 和声子结构信息。

计算路线
~~~~~~~~~~~~~~~~~~~~~~~~

目前支持两种非平衡动力学路线：

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - 路线
     - 主要功能
   * - 空间均匀 TD-BTE
     - 求解非平衡声子占据数的时间演化，分析声子模式间的能量交换与热化。
   * - 一维空间 TD-BTE
     - 在一维周期实空间中同时考虑声子传播和非线性三声子碰撞。

两种路线共享三声子碰撞核。
空间输运通过可选的 ``tdbte.spatial`` 配置启用，
声子 q 网格仍为三维。

常用参数
~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: TD-BTE 参数
   :header-rows: 1
   :widths: 35 65

   * - 参数
     - 说明
   * - ``force-constants``
     - 已有 FC2/FC3 及匹配元数据的来源目录。
   * - ``kernel``
     - 已生成的碰撞核，用于避免重新构建。
   * - ``mesh``
     - 碰撞积分使用的声子 q 网格。
   * - ``temperature``
     - 初始 Bose 分布的参考温度。
   * - ``branches``
     - 指定激发的声子分支，索引从零开始。
   * - ``excitation``
     - 选定模式占据数的相对增加幅度。
   * - ``duration-ps``
     - 总演化时间，单位 ps。
   * - ``max-step-ps``
     - 最大时间积分步长，单位 ps。
   * - ``samples``
     - 保存结果的时间点数量。
   * - ``spatial``
     - 一维空间输运的长度、网格、边界及初始条件设置。

均匀 TD-BTE 中，``excitation`` 将指定分支的初始占据数乘以
``1 + excitation``，不是温升或激光能流密度。

.. note::

   当前 TD-BTE 不包含四声子、Wigner 相干项或电子耦合。

空间均匀动力学见 :doc:`tdbte`，
一维空间输运见 :doc:`tdbte_spatial`。


``plot``：计算结果绘图
----------------------------------------

``plot`` 用于控制声子性质、散射率、热导率和非平衡动力学结果的可视化。
绘图直接读取已有数据，不会重新执行原子力或热输运计算。

.. list-table:: 绘图参数
   :header-rows: 1
   :widths: 35 65

   * - 参数
     - 说明
   * - ``layout``
     - ``separate``、``combined`` 或 ``both``，控制单图及组合图。
   * - ``path``
     - 声子色散的高对称路径，可使用 ``seekpath`` 自动生成。
   * - ``path-points`` / ``path-segments``
     - 自定义色散路径。
   * - ``temperature``
     - 散射率绘图所用温度，默认从已有结果中选择匹配数据。
   * - ``kappa``
     - ``x``、``y``、``z`` 或 ``all``，控制热导率分量。
   * - ``tau``
     - 选择绘制的散射机制。
   * - ``scattering-rate-min``
     - 散射率显示下限，单位 ps⁻¹。
   * - ``dpi``
     - PNG 图像分辨率。
   * - ``tdbte-times``
     - 均匀 TD-BTE 色散快照的三个时间点。

绘图内容与数据要求
~~~~~~~~~~~~~~~~~~~~~~~~

不同类型的图像需要不同的计算结果：

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - 图像
     - 所需数据
   * - 声子色散与 DOS
     - FC2 及匹配的结构元数据。
   * - 体积热容与群速度
     - FC2 及对应声子信息。
   * - 声子散射率
     - 声子线宽或散射率数据。
   * - 晶格热导率
     - 已完成的热输运计算结果。
   * - 非平衡动力学
     - TD-BTE 保存的时间演化结果。

``layout: both`` 同时生成单独图像和组合图。

散射率绘图
~~~~~~~~~~~~~~~~~~~~~~~~

``tau`` 支持以下选项：

- ``phph``：声子–声子散射。
- ``isotope``：同位素散射。
- ``total``：声子–声子、同位素及二者之和。
- ``normal`` / ``umklapp``：N / U 散射。
- ``nu``：同时展示 N 和 U 散射。
- ``all``：绘制所有可用散射信息。

程序使用 ``4*pi*gamma`` 将 phono3py 线宽转换为散射率。
``total`` 不包含边界散射。
使用 ``nu`` 时需要同时保存 N 和 U 散射数据。

对于 FourPhonon，三声子与四声子散射率分别绘制，
不会在图中自动求和。

热导率绘图
~~~~~~~~~~~~~~~~~~~~~~~~

普通 LBTE 结果若包含 RTA 数据，可在同一图中比较 LBTE 和 RTA。
Wigner 输运分别展示粒子贡献、相干贡献及总热导率。

``kappa: all`` 表示三个对角热导率分量的平均值，
而不是热导率张量的全部分量。


``output``：结果输出
----------------------------------------

``output`` 用于设置计算结果目录与运行进度显示。

.. list-table:: 输出参数
   :header-rows: 1
   :widths: 35 65

   * - 参数
     - 说明
   * - ``result-dir``
     - 计算结果保存目录。
   * - ``progress``
     - 是否显示动态计算进度。

建议为不同材料、温度或 q 网格设置独立的结果目录，
方便保存、比较和复用数据。

运行日志保存在结果目录的 ``run.log`` 中。
计算完成后，可通过 ``nepkappa plot`` 和 ``nepkappa report``
读取已有数据，生成图像和计算报告。


``parallel``：并行计算
----------------------------------------

``parallel`` 为可选模块，用于配置力常数与热输运阶段的并行资源。
NEP-kappa 支持 Slurm 任务调度，并可分别设置不同阶段的计算资源。

以下示例使用 Slurm 分配力常数任务和 phono3py 输运任务：

.. code-block:: yaml

   parallel:
     force-constant:
       backend: slurm
       jobs: 16
       submit: false

     kappa:
       backend: slurm
       jobs: 4
       max-concurrent: 4
       cpus-per-task: 8
       collect-cpus-per-task: 8
       submit: false

.. list-table:: 常用 Slurm 参数
   :header-rows: 1
   :widths: 35 65

   * - 参数
     - 说明
   * - ``backend``
     - 选择计算执行方式，例如 ``slurm``。
   * - ``jobs``
     - 分配的工作任务数量。
   * - ``max-concurrent``
     - 同时运行的最大任务数。
   * - ``cpus-per-task``
     - 每个任务使用的 CPU 核心数。
   * - ``collect-cpus-per-task``
     - 结果收集阶段使用的 CPU 核心数。
   * - ``submit``
     - 是否自动提交 Slurm 任务。
   * - ``extra-sbatch``
     - 传递给 Slurm 的额外资源选项。

``submit: false`` 只生成 Slurm 脚本而不调用 ``sbatch``，
但准备阶段仍可能执行部分本地工作。

phono3py 的 RTA 与 LBTE 均支持将不可约 q 点分配到多个 Slurm 任务。
