
TD-BTE：非平衡声子输运
========================================

上一章通过空间均匀 TD-BTE 研究声子模式之间的能量交换。
本章进一步引入一维实空间坐标，
同时考虑声子传播与非线性三声子碰撞。

NEP-kappa 使用三维声子 q 网格和一维周期实空间网格，
模拟温度、热流和声子分布的时空演化。

本章以 3C-SiC 的瞬态热光栅（TTG）为例，
介绍空间 TD-BTE 的输入、计算方法和结果分析。


计算方法
------------------------

空间 TD-BTE 同时包含声子漂移和碰撞：

.. math::

   \frac{\partial n_\lambda(x,t)}{\partial t}
   +
   v_{\lambda,x}
   \frac{\partial n_\lambda(x,t)}{\partial x}
   =
   \mathcal C_\lambda^{3\mathrm{ph}}[n(x,t)].

其中：

- :math:`n_\lambda(x,t)` 为声子占据数的时空分布。
- :math:`v_{\lambda,x}` 为群速度在输运方向上的分量。
- :math:`\mathcal C_\lambda^{3\mathrm{ph}}` 为非线性三声子碰撞积分。

与空间均匀 TD-BTE 相比，
这里增加了声子漂移项，
用于描述声子携带能量在不同空间位置之间传播。

.. note::

   本章的“一维”指实空间输运方向只有一个，
   声子仍然使用完整三维 q 网格。

   这与一维纳米线所用的 ``mesh: [N, 1, 1]``
   碰撞积分模型不同。


瞬态热光栅
------------------------

瞬态热光栅（Transient Thermal Grating, TTG）
使用空间周期性的初始温度分布：

.. math::

   T(x,0)
   =
   T_0+\Delta T
   \cos\left(\frac{2\pi x}{L}\right).

其中，:math:`T_0` 为背景温度，
:math:`\Delta T` 为初始温度调制幅度，
:math:`L` 为温度调制周期。

初始时，各空间单元中的声子服从对应局部温度的 Bose-Einstein 分布。
随后系统在没有持续外部加热的条件下自由演化。

声子传播和碰撞共同作用，
使初始温度分布随时间改变。


准备计算
------------------------

延续前一章的工作目录：

.. code-block:: bash

   conda activate nepkappa
   cd ~/work/3C-SiC

确认已有 FC2、FC3 和声子元数据：

.. code-block:: bash

   ls results/fc2.hdf5
   ls results/fc3.hdf5
   ls results/phono3py_disp.yaml

空间 TD-BTE 使用这些数据计算声子频率、
群速度与三声子碰撞核。

与空间均匀模型一样，
本章不需要预先计算晶格热导率。


配置输入
------------------------

创建 ``tdbte-spatial.yaml``：

.. code-block:: yaml

   # 空间非平衡声子动力学
   tdbte:
     force-constants: results        # 复用已有 FC2/FC3
     mesh: [3, 3, 3]                # 三维声子 q 网格
     temperature: 300               # 背景温度（K）
     duration-ps: 0.1               # 总演化时间（ps）
     max-step-ps: 0.01              # 最大时间步长（ps）
     samples: 6                     # 保存时刻数量

     # 一维实空间输运
     spatial:
       length-nm: 64                # 空间周期长度（nm）
       cells: 16                    # 实空间网格单元数
       direction: [1, 0, 0]        # 笛卡尔输运方向
       boundary: periodic           # 周期边界条件
       cfl: 0.8                     # CFL 数
       batch-size: 8                # 碰撞计算的单元批大小

       # 初始温度分布：瞬态热光栅
       initial:
         profile: ttg
         modes: all
         excitation: temperature
         delta: 10   # 初始温度调制幅度

   # 绘图设置
   plot:
     layout: both
     dpi: 150

   # 结果目录
   output:
     result-dir: results/tdbte-spatial

本例使用 3×3×3 的三维声子 q 网格，
同时沿 x 方向建立长度为 64 nm 的一维周期实空间。

实空间包含 16 个单元，网格间距为：

.. math::

   \Delta x=\frac{L}{N_x}=4\ \mathrm{nm}.

初始温度在 300 K 附近呈余弦分布，
调制幅度为 10 K。

该输入使用较小的 q 网格和较短的演化时间，
用于展示空间 TD-BTE 的计算流程。


主要空间参数
~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: 空间 TD-BTE 参数
   :header-rows: 1
   :widths: 36 64

   * - 参数
     - 说明
   * - ``length-nm``
     - 一维周期空间的长度，单位 nm。
   * - ``cells``
     - 实空间网格单元数。
   * - ``direction``
     - 输运方向的笛卡尔向量。
   * - ``boundary``
     - 空间边界条件，目前使用 ``periodic``。
   * - ``cfl``
     - 控制漂移积分稳定性的 CFL 数。
   * - ``threads``
     - 空间单元碰撞并行线程数，默认 1；大于 1 时需要 Numba。
   * - ``batch-size``
     - 碰撞积分的空间单元处理批大小。
   * - ``initial``
     - 初始温度或声子占据数的空间分布。

``length-nm`` 与 ``cells`` 决定实空间网格，
而 ``mesh`` 决定声子碰撞积分所用的 q 网格。
两者分别控制空间分辨率与声子模式分辨率。


运行计算
------------------------

检查输入并执行空间 TD-BTE：

.. code-block:: bash

   nepkappa info tdbte-spatial.yaml --for tdbte
   nepkappa tdbte tdbte-spatial.yaml

程序读取力常数和声子群速度，
构建三声子碰撞核并设置初始温度分布。

随后在各空间单元之间进行声子输运，
同时计算每个单元内部的非线性三声子碰撞。

计算结果保存至 ``results/tdbte-spatial``。


多核碰撞计算
------------------------

空间 TD-BTE 可以在多个 CPU 线程上计算不同空间单元的碰撞积分。
在 ``tdbte.spatial`` 下设置：

.. code-block:: yaml

   threads: 8
   batch-size: 16

各线程处理不同空间单元，共享当前事件分块。
单元内部仍按原有顺序累加散射事件，避免多个线程同时写同一模式的导数。
默认 ``threads: 1`` 使用串行实现；多线程需要 Numba，
后端不可用或线程数超过 Numba 上限时会报错。
申请的线程数记录在 ``audit.json`` 的 ``spatial_options`` 中。

单节点 Slurm 作业应申请 ``--ntasks=1 --cpus-per-task=8``，
并在启动 Python 前设置 ``NUMBA_NUM_THREADS=8``。
BLAS 线程保持为 1，避免嵌套并行；仅修改 ``OMP_NUM_THREADS`` 不会开启此功能。
``batch-size`` 限制一批可以同时处理的单元数。
对于 16 个空间单元，可先测试 ``threads: 8``、``batch-size: 16``；
线程数多于每批单元数时无法增加单元并行度。

时间积分、漂移及事件分块读取仍为串行。
加速比取决于碰撞计算量、内存带宽和存储速度，不能保证 8 线程加速 8 倍。
无空间维度的均匀 TD-BTE 仍使用串行的单状态碰撞计算。
此选项不包含 MPI 或跨节点并行。


时间推进方法
------------------------

空间 TD-BTE 使用算子分裂方法，
分别处理碰撞和漂移两个过程。

一个时间步按照以下顺序执行：

1. 半个时间步的非线性三声子碰撞。
2. 一个完整时间步的空间漂移。
3. 半个时间步的非线性三声子碰撞。

碰撞部分使用 DOP853 时间积分，
空间漂移使用守恒的一阶迎风格式。

时间步长同时受到 ``max-step-ps`` 和 CFL 条件限制：

.. math::

   \max_\lambda |v_{\lambda,x}|
   \frac{\Delta t}{\Delta x}
   \leq \mathrm{cfl}\leq1.

当前漂移格式为一阶迎风格式，
因此整体方法具有一阶时间与空间精度。
计算时应检查网格和 CFL 引起的数值扩散。


结果与分析
------------------------

完成计算后生成图像和报告：

.. code-block:: bash

   nepkappa plot tdbte-spatial.yaml
   nepkappa report results/tdbte-spatial

图像保存在 ``results/tdbte-spatial/plots`` 目录。

.. .. figure:: ../source/_static/tdbte_spatial.png
..    :width: 95%
..    :align: center
..    :alt: 3C-SiC spatial TD-BTE transient thermal grating

..    3C-SiC 的空间 TD-BTE 计算结果，
..    展示瞬态热光栅中温度、热流和声子分布的时空演化。


温度的时空演化
~~~~~~~~~~~~~~~~~~~~~~~~

空间 TD-BTE 输出每个位置和时刻的声子占据数。
程序根据各单元中的总声子热能反演能量等效温度：

.. math::

   E(x,t)=E_{\mathrm{BE}}[T(x,t)].

其中，:math:`T(x,t)` 为局部能量等效温度，
并不要求每个空间单元内部的声子分布已经达到局部平衡。

对于 TTG 初始条件，可以观察周期性温度扰动如何随时间演化。

温度调制的衰减受到声子传播、模式碰撞
以及空间尺度的共同影响。


热流与能量输运
~~~~~~~~~~~~~~~~~~~~~~~~

空间 TD-BTE 同时计算声子热流的时空分布。

不同声子模式根据其频率、群速度与占据数
对热流产生贡献。

通过分析温度和热流的空间分布，
可以研究非平衡声子输运与普通扩散行为的差异。

程序还保存各声子分支的能量密度和热流，
用于区分不同分支的输运贡献。


TTG 温度谐波
~~~~~~~~~~~~~~~~~~~~~~~~

初始 TTG 为单一空间余弦扰动。
随着时间演化，程序对温度分布进行 Fourier 分解，
记录各空间谐波的时间变化。

对于基频：

.. math::

   A_1(t)
   =
   \frac{2}{L}
   \int_0^L
   [T(x,t)-\overline T(t)]
   \cos\left(\frac{2\pi x}{L}\right)\,dx.

通过分析 :math:`A_1(t)`，
可以研究不同空间周期下温度调制的衰减过程。

在扩散极限且热扩散率恒定时，
TTG 基频满足：

.. math::

   A_1(t)
   =
   A_1(0)
   \exp\left[
   -\alpha\left(\frac{2\pi}{L}\right)^2t
   \right].

其中，:math:`\alpha` 为热扩散率。

这一表达式可作为扩散行为的参考，
而空间 TD-BTE 可用于研究偏离简单扩散的输运过程。


独立设置初始条件
------------------------

``spatial.initial`` 将空间分布、模式选择和激发方式分开设置：

.. list-table:: 初始条件参数
   :header-rows: 1
   :widths: 25 40 35

   * - 参数
     - 可选值
     - 含义
   * - ``profile``
     - ``ttg``、``gaussian``、``uniform``
     - 空间包络，默认 ``ttg``。
   * - ``modes``
     - ``all``、``branches``、``frequency-window``
     - 模式选择，默认 ``all``。
   * - ``excitation``
     - ``temperature``、``relative``
     - 占据初始化方式，默认 ``temperature``。

三组参数共 18 种组合均受支持。``uniform`` 专指空间均匀，
``all`` 表示所有有效声子模式。分支和频率筛选作用于碰撞核所含的全部 q 点，
零频刚体模式保持固定。

以 :math:`p(x)` 表示空间包络，:math:`S` 表示选中的模式。
TTG 为 :math:`p(x)=\cos(2\pi x/L)`，空间均匀为 :math:`p(x)=1`。
高斯包络采用周期高斯，必须给出宽度 ``sigma-nm``，
中心 ``center-nm`` 默认为 :math:`L/2`，并在中心归一化为 1。
选中模式 :math:`\lambda\in S` 的两种初始化方式为：

.. math::

   n_\lambda(x,0)=
   \begin{cases}
   n_{\mathrm{BE}}(\omega_\lambda,T_0+\Delta T p(x)),
      &\text{temperature},\\
   n_{\mathrm{BE}}(\omega_\lambda,T_0)[1+A p(x)],
      &\text{relative}.
   \end{cases}

未选中模式保持背景温度 :math:`T_0` 的 Bose 占据。
选择部分模式时，``temperature`` 指定这些模式初始 Bose 占据对应的温度，
并不等于整个空间单元的能量等效温度。
即使选择所有模式，相对占据增幅通常也会产生非热频谱分布。

两种方式统一使用 ``delta``：
``excitation: temperature`` 时单位为 K，默认 1 K；
``excitation: relative`` 时为无量纲占据变化比例，默认取 ``tdbte.excitation`` （0.01）。
例如 ``delta: 0.2`` 分别表示 0.2 K 或 20%。
TTG 要求 :math:`|\Delta T|<T_0` 或 :math:`|A|<1`；
余弦为负的半周期会降低相对于背景的占据。
高斯或均匀空间分布要求 :math:`T_0+\min(\Delta T,0)>0` 或 :math:`A>-1`。

指定分支的 TTG 温度激发
~~~~~~~~~~~~~~~~~~~~~~~~

将 ``spatial.initial`` 改为：

.. code-block:: yaml

   initial:
     profile: ttg
     modes: branches
     branches: [3, 4, 5]  # 从 0 编号，需核对材料和碰撞核。
     excitation: temperature
     delta: 10

改为 ``modes: all`` 并删除 ``branches``，即为全模式热 TTG。

指定频率窗口的高斯激发
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   initial:
     profile: gaussian
     sigma-nm: 8
     center-nm: 32
     modes: frequency-window
     frequency-range-THz: [10, 15]
     excitation: relative
     delta: 0.2

频率窗口包含两个端点，要求 :math:`0\le f_{min}<f_{max}`。
只有窗口内的有效模式改变占据，在高斯中心增加 20%。
程序拒绝空窗口选择和碰撞核中不存在的分支。
改为 ``profile: uniform`` 并删除高斯参数，即为空间均匀激发。
不同频段采用同一相对增幅时，注入能量一般不同。

旧输入兼容
~~~~~~~~~~~~~~~~~~~~~~~~

原有 ``kind: ttg``、``gaussian``、``uniform`` 保留全模式温度初始化行为。
``kind: selective`` 保留指定分支的相对占据激发，空间 ``profile`` 可选高斯或均匀。
旧热初态的 ``delta-temperature-K`` 和选择性初态的 ``amplitude`` 会转换为 ``delta``；
同时填写 ``delta`` 和对应的旧参数会报错。
新组合统一使用 ``delta`` 及上述独立参数，并删除 ``kind``。
旧 ``kind`` 与新 ``modes`` 或 ``excitation`` 混用将报错。
旧热初态下曾被忽略的 ``branches`` 或 ``profile`` 现在会报错，
避免误以为只激发了部分模式。


主要输出
------------------------

空间 TD-BTE 将时空分辨结果保存为：

.. code-block:: text

   results/tdbte-spatial/tdbte/spatial.hdf5

.. list-table:: 主要输出数据
   :header-rows: 1
   :widths: 43 57

   * - 数据集
     - 内容
   * - ``occupation``
     - 各位置、时刻和声子模式的占据数。
   * - ``temperature_K``
     - 局部能量等效温度。
   * - ``heat_flux_eV_nm2_ps``
     - 空间分辨声子热流。
   * - ``branch_energy_eV_nm3``
     - 各声子分支的能量密度。
   * - ``branch_heat_flux_eV_nm2_ps``
     - 各声子分支的热流。
   * - ``branch_collision_power_eV_nm3_ps``
     - 三声子碰撞引起的分支间能量交换率。
   * - ``ttg_cosine_amplitude_K``
     - TTG 基频余弦温度幅值。
   * - ``departure_from_local_BE``
     - 与局部等能 Bose 分布的偏离程度。

此外，``audit.json`` 记录能量守恒、分支能量收支、
熵、占据数正性和 CFL 等数值检查。

程序生成温度、热流、分支能量及温度谐波的图像，
用于分析空间非平衡输运过程。
