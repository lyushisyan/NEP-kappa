时变声子 BTE 与非平衡弛豫
====================================

TD-BTE 在固定声子频率下，利用非线性三声子散射演化占据数。
本章关注空间均匀的选择性光学支激发：没有空间漂移、外部热浴或持续泵浦。
:doc:`tdbte_spatial` 增加周期一维实空间输运，保留三维 q 网格和非线性碰撞。

.. important::

   当前实现通过离散守恒等数值检查，但生产碰撞核的绝对弛豫速率尚未完成独立物理验证。
   2026-11-07 审计记录了小 q 声学端点的离壳差异；该诊断本身不等于在壳散射率发散，
   也不阻止空间求解器的软件开发。详细约定、数据和限制见 :doc:`tdbte_validation`。

.. toctree::
   :maxdepth: 1

   tdbte_spatial
   tdbte_validation

准备力常数并执行
--------------------------

复用快速开始生成的配套 ``fc2.hdf5``、``fc3.hdf5`` 和 ``phono3py_disp.yaml``。
元数据提供结构、超胞和原胞映射，热导率结果不是前提。

.. code-block:: bash

   mkdir -p ~/work/3C-SiC/tdbte
   cd ~/work/3C-SiC/tdbte
   ls -lh ../results/fc2.hdf5 ../results/fc3.hdf5 ../results/phono3py_disp.yaml

下载 :download:`tdbte.yaml <../source/inputs/tdbte.yaml>`：

.. literalinclude:: ../source/inputs/tdbte.yaml
   :language: yaml

.. code-block:: bash

   nepkappa info tdbte.yaml --for tdbte
   nepkappa tdbte tdbte.yaml

示例的 6×6×6 网格包含 216 个 q 点。对于两原子 3C-SiC，0–2 为声学支、3–5 为光学支；
编号从零开始，在各 q 点按频率排序，并不是适用于任意材料的光学支编号。
``excitation: 0.01`` 将所选支的初始 Bose 占据数增加 1%，不等于温度增加 1 K。

``samples: 201`` 决定 20 ps 内的保存时刻数，不决定内部积分步长。
求解器采用自适应 DOP853，``max-step-ps`` 是最大步长。
长工作站计算可使用 :ref:`tmux <transport-tmux>`，或在外层 Slurm 分配中执行相同 CLI。

从头生成力常数
--------------------------

仓库的 ``examples/tdbte.yaml`` 是完整六部分输入，用 ``nepkappa run`` 先生成配套 FC2/FC3，再构建碰撞核并积分。
同一文件交给 ``nepkappa tdbte`` 时，直接读取结果目录中已经存在的力常数，不补算原子力。
当前 TD-BTE 只使用 phono3py FC2/FC3；QHA、SSCHA、FC4 和 Wigner 设置不会自动转入该求解器。

关键参数
--------------------------

.. list-table:: 均匀动力学参数
   :header-rows: 1
   :widths: 32 68

   * - 参数
     - 约定
   * - ``force-constants``
     - 已有力常数目录，不是生成原子力的请求。
   * - ``kernel``
     - 已完成的碰撞核；指定它时省略力常数来源和网格。
   * - ``mesh``
     - 碰撞积分的完整、无偏移 q 网格，与 ``kappa.mesh`` 无关。
   * - ``temperature``
     - 初始 Bose 参考温度，默认 300 K。
   * - ``branches``
     - 明确指定零起始分支编号；跨交叉不跟踪极化身份。
   * - ``excitation``
     - 所选初始占据数乘以 ``1 + excitation``，默认 0.01。
   * - ``duration-ps`` / ``max-step-ps`` / ``samples``
     - 默认 20 ps / 0.5 ps / 201 个保存时刻。

三维网格三个分量均至少为 2；实验一维网格仅一个分量大于 1，另外两个为 1。
二维网格不支持。均匀积分容差为 ``rtol=1e-11``、``atol=1e-14``，
时间步长与 q 网格收敛应分开检查。该初态不规定激光能流密度，也不为各模式注入相同能量。

碰撞核、归一化与资源
--------------------------

构建器逐个母 q 点计算 phono3py 原生最短矢量三声子相互作用强度。
三维网格单元使用六个 Freudenthal 四面体，对线性插值能量壳采用正权重、对称二次曲面求积。
有序子模式事件系数为 ``4*pi*U*strength*surface_weight``，其中 ``U`` 为原生虚部自能单位转换。
原生强度已包含 ``1/Nq``，碰撞积分不能再次乘相同网格归一化。

频率与熵变量使用相同重心权重插值：

.. math::

   \eta_\lambda=\ln\frac{1+n_\lambda}{n_\lambda}.

子模式占据数更新用对应的伴随沉积，保证离散能量守恒与 Bose 细致平衡。
这个弱离散算子不等同于 phono3py 原生四面体碰撞算子，也不是 RTA 指数衰减。

事件以最多 20,000 个一块写入 ``tdbte-kernel/``。每次导数计算只映射、归约和关闭一个块，
不把全体事件或稠密的模式碰撞矩阵放进内存。可用时自动使用 Numba，并在积分前与 NumPy 参考值比较。
审计记录实际后端，编译加速不改变物理模型。

分块不能消除细网格代价。原生顶角计算约随 q 点数平方和支数立方增长，
单个母 q 的原生强度、声子和力常数仍需内存。
未压缩事件数组约 112 字节 / 事件，数千万事件可能占数 GB 磁盘，反复读取也可能受 I/O 限制。
先用较小网格估计规模，必要时使用计算节点本地 scratch。

稳定性与声学端点
--------------------------

三维构建要求 Gamma 仅有三个近零刚性平移（阈值 1e-4 THz），并固定为零。
其他虚频或额外近零模式被拒绝，不通过裁剪或静默对称化修复。
这些刚性腿的熵变量和顶角固定为零，其声学极限近似仍需独立检验。
非零耦合下的平坦共振四面体会报错。

当前不支持 NAC、偏移网格和广义网格；不要删除 Born 电荷数据以绕过检查。

实验性一维链 / 纳米带
--------------------------

``mesh: [51,1,1]`` 沿第一倒易坐标进行一维碰撞积分，其他方向必须在物理上由足够真空隔离。
位移分量和所有声子支仍保留。这仍是空间均匀占据数模型，不会沿纳米带传播热脉冲或施加端部热库。

.. code-block:: yaml

   tdbte:
     force-constants: calculations/graphene-ribbon/force-constants
     mesh: [51, 1, 1]
     temperature: 300
     branches: [4]
     excitation: 0.01
     duration-ps: 20
     max-step-ps: 0.5
     samples: 201
   output:
     result-dir: calculations/graphene-ribbon/tdbte-1d

分支 4 只是示例，应根据实际谱选择。程序不会从 q 网格自动判断真空是否充分。
一维路线用区间代替四面体，在分段线性能量差的简单根上使用斜率倒数与共享端点半权重。
占据数方程不再引入横截面积因子。

Gamma 允许三个刚性平移，或再加一个轴向刚性转动。质量加权近零本征子空间需在 ``1e-5`` 容差内
匹配相应刚性运动。转动测试要求原子位置在真空两侧形成连续纳米带表示。
额外软模、不稳定模和非刚性零空间会被拒绝。

零线性顶角插值区间记入 ``zero_strength_intervals``，不代表区间内部真实散射一定为零。
近切根、支交叉、转动和声学端点都需进一步验证；石墨烯纳米带的绝对弛豫速率尚无完整基准。

复用与失败处理
--------------------------

完成构建后发布 ``manifest.json``，包含源哈希、网格、模式排序、后端版本、归一化、边界假设与文件校验和。
丢失、重复、修改或未完成的数据块会在传播前被拒绝。运行时将核文件视为不可变数据。

.. code-block:: yaml

   tdbte:
     kernel: calculations/SiC/time-bte/tdbte-kernel/manifest.json
     temperature: 300
     branches: [3, 4, 5]
   output:
     result-dir: calculations/SiC/new-excitation

网格属于核，不能另行指定。自动缓存复用检查来源、网格和构建 / 后端版本。
失败构建不自动续算；保留 ``build-failure.json`` 和旧目录，修正原因后使用新目录。
已有 ``result-dir/tdbte`` 不会被覆盖。

旧的单文件 NPZ 加同名 JSON sidecar 仍可读。研究原型的 chunk manifest 不能仅改格式名冒充版本化产物。
单文件核包含 ``frequency_THz``、``zero_modes``、``parent``、``daughter1``、``daughter2``、
``weights``、``coefficients``；子模式数组为 (事件数, V)，三维 V=4，一维 V=2。
权重非负且逐行和为 1；频率、系数非负，索引为整数。
sidecar 需有 ``model``、``coefficient_time_unit: ps``、明确的 ``mode_branch_indices`` 与 ``uniform_mode_weight``。
不支持不等权重不可约网格；读取器不能替任意外部核证明归一化声明真实可靠。

输出、绘图与解释
--------------------------

.. code-block:: bash

   tail -n 20 results/optical-m666/run.log
   nepkappa info tdbte.yaml --for plot
   nepkappa plot tdbte.yaml

``tdbte/`` 中：

* ``trajectories.npz`` 保存占据数、频率、时刻、支编号、熵和独立积分的未激发对照。
* ``branch-energy.csv`` 保存相对初始 Bose 态的支超额能量，单位 eV / 原胞，含 ``1/Nq`` 权重。
* ``audit.json`` 记录来源、核身份、积分设置、注入能量、能量漂移、平衡漂移和采样熵检查。

均匀 TD-BTE 六面板上排是三个时间的超额能量着色色散，下排是分支超额能量、声学支与光学支能量等效温度。
着色量为 ``h*f*(n(t)-n_Bose(T0))``，单位 eV。三张色散共用色标，频率保持固定。
周期三线性插值只将网格上的能量映射到路径，不会额外计算路径上的动力学。
分支温度通过各支的热能匹配 Bose 分布得到，不表示该支内部已经达到平衡。

默认取起点、中点和终点的最近保存时刻，``plot.tdbte-times`` 可指定递增且位于保存区间内的三个时间。
输出包括 ``tdbte_combined.png``、六个单图、``tdbte_branch_temperature.csv`` 和绘图元数据。
完整绘图会核对来源哈希、频率和审计；CSV-only 旧结果只支持能量图。

占据数不会被静默裁剪。失败审计保留诊断文件并报错，不能作为接受的轨迹。
守恒检查即使在错误的全局散射前因子下也可能成立，因此不能替代绝对速率和网格收敛验证。
孤立系统不要求回到初始 300 K；有限时窗内各支温度不同也不能称为完全平衡。
当前未实现频率 / q 窗口激发、规定能量注入、持续泵浦、checkpoint 重启或专用 TD-BTE Slurm 调度。
