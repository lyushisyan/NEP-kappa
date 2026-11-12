谐近似声子与温度重整化
================================

本章介绍谐近似、准谐近似（QHA）、固定体积随机自洽谐近似（SSCHA），以及 QHA 平衡体积上的 SSCHA。
热导率计算见 :doc:`thermal_transport`。

准备与路线选择
------------------------

延续快速开始的工作目录：

.. code-block:: bash

   mkdir -p ~/work/3C-SiC/phonons
   cd ~/work/3C-SiC/phonons
   cp ../POSCAR .
   cp ../nep_3C-SiC.txt .

下列示例使用相同模型、独立的结果目录。下载文件中的相对路径以当前目录为准。

.. list-table:: 各路线的区别
   :header-rows: 1

   * - 路线
     - 新增计算与含义
     - 命令
   * - 谐近似
     - 先弛豫，再生成 FC2。
     - ``relax`` → ``fc2``
   * - QHA
     - 扫描体积，拟合自由能，再计算平衡体积声子。
     - ``fc2``
   * - 固定体积 SSCHA
     - 从初始 FC2 出发，在各温度热采样并自洽更新辅助 FC2。
     - ``run``
   * - QHA+SSCHA
     - QHA 给出体积，各温度再执行固定体积 SSCHA。
     - ``run``

谐近似声子
------------------------

下载 :download:`harmonic.yaml <../source/inputs/harmonic.yaml>`：

.. literalinclude:: ../source/inputs/harmonic.yaml
   :language: yaml

.. code-block:: bash

   nepkappa info harmonic.yaml --for relax
   nepkappa relax harmonic.yaml
   nepkappa info harmonic.yaml --for fc2
   nepkappa fc2 harmonic.yaml
   nepkappa info harmonic.yaml --for plot
   nepkappa plot harmonic.yaml

启用弛豫后，单独的 ``fc2`` 读取 ``results/harmonic/POSCAR_relaxed``。
这里只生成 FC2，不需要 FC3。结果包含匹配的晶胞与位移元数据。
``layout: both`` 生成色散、DOS、体积定容热容、群速度四幅单图和 2×2 ``combined.png``。
q 网格用于布里渊区积分，色散路径控制曲线展示，两者用途不同。

QHA：热膨胀及体积引起的频移
------------------------------------

下载 :download:`qha.yaml <../source/inputs/qha.yaml>`。示例从已弛豫的谐近似结构开始，扫描五个体积比。

.. literalinclude:: ../source/inputs/qha.yaml
   :language: yaml

.. code-block:: bash

   nepkappa info qha.yaml --for fc2
   nepkappa fc2 qha.yaml
   nepkappa info qha.yaml --for plot
   nepkappa plot qha.yaml

``equilibrium-phonons: true`` 在拟合的各温度平衡体积上额外计算 FC2。
``run-transport: false`` 表示本例到声子和热力学性质为止。
体积扫描必须包围整个温度区间中的平衡体积。

总览包含自由能—体积曲线、平衡体积—温度、体积热膨胀和温度相关色散 / DOS。
输出为 PNG 和 PDF，包括 ``qha_combined`` 与四个单图。
扫描中体积比为 1 的点提供固定体积参考。这里的频移来自体积改变。

SSCHA：固定体积的声子重整化
------------------------------------

下载 :download:`sscha.yaml <../source/inputs/sscha.yaml>`：

.. literalinclude:: ../source/inputs/sscha.yaml
   :language: yaml

.. code-block:: bash

   nepkappa info sscha.yaml --for run
   nepkappa run sscha.yaml
   nepkappa info sscha.yaml --for plot
   nepkappa plot sscha.yaml

当前路线使用 Phonopy 的随机 SSCHA 迭代，由 NEP 提供总能和原子力，在固定体积下更新 FC2。
示例采样 300、600、900 K；快照数、迭代数、超胞与采样网格都需要收敛检查。

禁用输运时，总览叠加所有请求温度的色散、DOS 和群速度。
每个热容点使用相应温度的 FC2 和体积，并配合固定体积谐近似参考。
``plot.temperature`` 不在此模式下挑选单一温度。

.. note::

   导出的 FC2 是辅助谐近似矩阵，不是自由能 Hessian；据此计算的热容属于辅助谐模型估计。
   可选 FC3 bubble 频移单独保存，不修改该 FC2，也不自动修改热导率。

QHA+SSCHA：体积膨胀与重整化
------------------------------------

下载 :download:`qha-sscha.yaml <../source/inputs/qha-sscha.yaml>`：

.. literalinclude:: ../source/inputs/qha-sscha.yaml
   :language: yaml

.. code-block:: bash

   nepkappa info qha-sscha.yaml --for run
   nepkappa run qha-sscha.yaml
   nepkappa info qha-sscha.yaml --for plot
   nepkappa plot qha-sscha.yaml

QHA 先确定平衡体积，再在各温度和对应体积上运行 SSCHA。SSCHA 迭代过程中该体积保持固定，
采样温度必须位于 QHA 范围内。这个顺序耦合过程不使用 SSCHA 自由能重新最小化体积。

谱图使用相同 QHA 体积上的谐近似作为对照。热容的灰色参考曲线另使用体积比为 1 的固定体积 QHA 点，
SSCHA 热容显示为独立彩色点，不用连线暗示连续自由能导数。
需要温度相关热导率时，查看 :doc:`input_files` 中的 ``run-transport`` 与温度网格关系。
