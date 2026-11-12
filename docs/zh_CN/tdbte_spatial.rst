空间非均匀非线性 TD-BTE
====================================

适用范围与输入
--------------------------

``tdbte.spatial`` 启用一维周期实空间网格，每个单元使用相同的非线性三声子碰撞积分。
声子 q 网格仍是三维。省略该映射即使用 :doc:`tdbte` 的均匀路线。
这里不采用 RTA、Callaway 或线性化时间演化，也没有四声子或电子耦合。

从仓库根目录可检查完整示例，随后按需要计算：

.. code-block:: bash

   nepkappa info examples/tdbte-spatial.yaml --for run
   nepkappa run examples/tdbte-spatial.yaml
   nepkappa plot examples/tdbte-spatial.yaml
   nepkappa report calculations/3C-SiC/tdbte-spatial

``run`` 先生成 FC2/FC3。已有力常数使用下面的精简输入和 ``nepkappa tdbte``，不重新生成原子力。

.. code-block:: yaml

   tdbte:
     force-constants: calculations/3C-SiC/fc
     mesh: [3, 3, 3]
     temperature: 300
     duration-ps: 0.1
     max-step-ps: 0.01
     samples: 6
     spatial:
       length-nm: 64
       cells: 16
       direction: [1, 0, 0]
       boundary: periodic
       cfl: 0.8
       batch-size: 8
       initial:
         kind: ttg
         delta-temperature-K: 10
   output:
     result-dir: calculations/3C-SiC/spatial-ttg

来源目录需要匹配的 ``phono3py_disp.yaml``、``fc2.hdf5`` 和 ``fc3.hdf5``。
也可用已有 ``kernel`` 并省略来源和网格，但计算群速度仍需其记录的匹配力常数。
``spatial.phonons`` 可指向迁移后的来源目录，三份文件的哈希都要一致。
缺少原生网格 / 来源身份的旧核不足以准备空间速度；NAC、偏移 / 广义网格和一 / 二维 q 网格均不支持。
示例设置不代表收敛，已有 ``result-dir/tdbte`` 会被拒绝覆盖。

初态
--------------------------

.. list-table:: ``spatial.initial`` 支持的模式
   :header-rows: 1
   :widths: 24 76

   * - ``kind``
     - 定义与参数
   * - ``ttg``
     - 瞬态热光栅：局部温度为 T₀ + ΔT cos(2πx/L)，各单元初始化为该温度的 Bose 分布。
       用 ``delta-temperature-K`` 指定 ΔT；至少四个单元，且 :math:`|\Delta T| < T_0`。
   * - ``gaussian``
     - 热高斯脉冲。指定温升、``sigma-nm`` 和可选 ``center-nm``，中心默认 L/2。
       使用周期镜像求和，使宽脉冲在周期接缝处仍平滑。
   * - ``selective``
     - 仅将指定 ``branches`` 的参考占据数乘以 1 + δg(x)。``amplitude`` 是相对占据数增加，
       不是开尔文或激光能流密度。默认高斯包络，也可 ``profile: uniform``。
   * - ``uniform``
     - 所有单元采用 ``temperature + delta-temperature-K`` 的 Bose 分布，用作静止数值对照。

例如热脉冲的初态：

.. code-block:: yaml

   initial:
     kind: gaussian
     delta-temperature-K: 30
     sigma-nm: 8
     center-nm: 32

两原子 3C-SiC 的选择性光学支脉冲示例：

.. code-block:: yaml

   initial:
     kind: selective
     branches: [3, 4, 5]
     amplitude: 0.2
     sigma-nm: 8
     center-nm: 32

分支编号仍按原生频率排序，不跨支交叉跟踪极化；激发简并组中单独一个成员可能依赖本征矢基。

速度、时间推进与单位
--------------------------

``direction`` 是笛卡尔三维向量，内部归一化，与原胞晶轴不同。
群速度由 FC2 动力学矩阵的解析导数给出，投影到与碰撞后端相同的原生本征矢基中。
倒易格矢平移会改变本征矢相位，因此不能把原生本征矢与仅取模的网格坐标任意混用。
输出同时保存索引地址和原生分数 q 坐标，模式排序为 q 在前、支在后。

原生 THz·Å 导数转换到 nm/ps 的因子是 0.1。固定 Gamma 平移的占据数和漂移速度均为零。
原胞体积使用 nm³；此体相三维 q 方程不引入真空厚度或横截面积修正。

每步先执行半步非线性碰撞，再执行整步迎风漂移，最后再执行半步碰撞。
实际步长同时受 ``max-step-ps`` 和 CFL 限制：

.. math::

   \max |v_x|\,\Delta t/\Delta x \leq \mathrm{cfl} \leq 1.

迎风漂移采用守恒凸组合，在 CFL 条件下保持非负。
``spatial.collisions: false`` 仅关闭碰撞，可作为弹道对照。
碰撞采用 DOP853，默认 ``rtol=1e-11``、``atol=1e-14``。
若产生非正的活动占据数，从原状态以更小内部步长重试，最多八次减半；失败报错，不裁剪占据数或重标度速率。
刚性碰撞可能使显式积分昂贵，当前没有隐式刚性求解器。

.. note::

   半碰撞—漂移—半碰撞的 Strang 排列，并不会使当前完整方法成为二阶。
   漂移采用前向 Euler 一阶迎风步，总体时间和空间误差仍为一阶，应检查网格和 CFL 引起的数值扩散。

``samples`` 决定保存时刻；``snapshot-times-ps`` 可显式指定从 0 到 ``duration-ps`` 的递增列表。
推进会落在保存时刻，密集输出可能缩短实际分裂步；时间步收敛比较应采用相同保存时刻。

所有单元共用不可变事件核。每次 RHS 计算将各单元分批处理，每个事件块只打开一次。
``batch-size`` 约束临时数组，不改变模型；状态和 DOP853 中间数组仍随单元数 × 模式数增长。

HDF5 输出与能量连续性
--------------------------

``tdbte/spatial.hdf5`` 只写请求的快照：

.. list-table:: 主要数据集
   :header-rows: 1
   :widths: 49 51

   * - 名称
     - 形状 / 物理意义
   * - ``occupation``
     - (时间, 单元, 模式)，无量纲占据数。
   * - ``temperature_K``
     - (时间, 单元)，由局部热能反演的等效温度。
   * - ``heat_flux_eV_nm2_ps``
     - 单元中心物理声子热流，单位 eV/(nm²·ps)。
   * - ``branch_energy_eV_nm3``
     - (时间, 单元, 分支)，热能密度，不含零点能。
   * - ``branch_heat_flux_eV_nm2_ps``
     - 各分支的单元中心物理热流。
   * - ``branch_collision_power_eV_nm3_ps``
     - 非线性碰撞导致的瞬时分支能量交换。
   * - ``branch_face_flux_integral_eV_nm2``
     - 自前一快照起累计的数值右界面通量。
   * - ``branch_collision_energy_integral_eV_nm3``
     - 两次碰撞半步产生的实际累计能量变化。
   * - ``branch_continuity_residual_eV_nm3``
     - 保存区间的离散分支能量收支残差。
   * - ``departure_from_local_BE``
     - 相对局部等能 Bose 态的能量加权 L1 偏差，除以局部能量。
   * - ``temperature_harmonics_real_K`` / ``temperature_harmonics_imag_K``
     - 前四个复 Fourier 系数，掩码排除 Nyquist 或混叠谐波。
   * - ``ttg_cosine_amplitude_K``
     - 基频系数实部的两倍，对应 TTG 余弦温度幅值。

物理常数、单位转换、坐标、模式身份、速度、本征矢和来源均存入文件。
热流乘 1.602176634×10¹¹ 得 W/m²，能量密度乘 1.602176634×10⁸ 得 J/m³。

分支预算使用漂移实际采用的界面通量：

.. math::

   e_s(t_j)-e_s(t_{j-1})
   +\nabla_x\cdot\int_{t_{j-1}}^{t_j}F_s^{\mathrm{upwind}}\,dt
   =\int_{t_{j-1}}^{t_j}P_s^{3\mathrm{ph}}\,dt.

因此局部散射与空间输运由同一离散收支分开，而不是从稀疏输出的两个差分导数相减得到。
单元中心物理热流和数值迎风界面热流是不同数据集。
``audit.json`` 记录能量、连续性、碰撞能量、熵、正性、CFL、RHS 次数和重试次数。
``plot`` 读取场量生成温度 / 热流时空图、分支能量和温度谐波，不加载全部占据数快照。

验证与已记录结果
--------------------------

可选本地脚本 ``tools/validate_spatial_tdbte.py`` 检查六模式非线性代数模型的弹道、网格 / 时间步收敛、
弱激发长波 Fourier 极限、脉冲和不同 TTG 幅值。``tools/`` 被 Git 忽略，新克隆中不包含该脚本。
模型系数任意，不是 SiC；其 Fourier 参考由 Bose 多项式和热容独立推导，不是拟合时间演化得到。

.. code-block:: bash

   OMP_NUM_THREADS=1 RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
     python tools/validate_spatial_tdbte.py \
     --output calculations/spatial-tdbte-validation

已有三维 SiC 核时可添加 ``--kernel`` 做材料数值演示；``--material-only`` 跳过代数模型。

:download:`2026-11-07 空间验证 JSON <../source/_static/tdbte-spatial-validation-20261107.json>` 记录：

* 弹道传播误差在 32 / 64 / 128 个单元上近似一阶下降。
* 固定空间网格、相对 0.00125 ps 参考，0.04 / 0.02 / 0.01 ps 的占据数 L2 误差为
  0.01089 / 0.005514 / 0.002871，显示整体一阶时间误差。
* 代数模型 300 K 的独立 Fourier 热扩散率为 0.507074 nm²/ps。
  64 nm、0.01 K 的 TTG 演化 20 ps，128 / 256 / 512 单元的衰减率误差为 41.86% / 17.26% / 4.95%。
* 热高斯、选择性高斯与 1 / 30 / 90 K TTG 均配置无碰撞对照，通过守恒、熵与正性检查。
* 八个 SiC 短算例使用 3×3×3 q 网格、162 模式、833,976 个能量壳事件、8 单元和 0.02 ps。
  最大相对总能漂移约 1.36e-16，连续性残差约 1.70e-16；这不代表材料输运已经收敛。

原生 SiC 沿输运方向的简并速度非对角诊断值为 0.527 nm/ps，输出记录该数值，标量占据数模型不包含对应相干。
所有算例仍标记 ``validated_for_physical_dynamics: false``。
绝对速率、q 网格、力常数超胞、空间 / 时间收敛及简并基近似仍需独立验证。
不能把无关的稳态热导率直接当作这个离散算子的 Fourier 极限参考。
