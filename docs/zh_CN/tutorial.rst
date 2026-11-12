进阶用法与结果复用
================================

分阶段执行
------------------------

普通三声子完整流程可拆为：

.. code-block:: bash

   nepkappa relax input.yaml
   nepkappa fc2fc3 input.yaml
   nepkappa kappa input.yaml
   nepkappa plot input.yaml
   nepkappa report input.yaml

只研究谐近似性质时用 ``fc2`` 替代 ``fc2fc3``。
若启用弛豫，单独 FC 阶段要求结果目录已有 ``POSCAR_relaxed``；已弛豫输入应显式关闭弛豫。
更改网格或温度时，可保留匹配的 FC2 / FC3 / 元数据，只运行 ``kappa``，不重新计算原子力。

需要同时使用两种格式时选 ``force-constant.format: both``。
``shengbte`` 文本本身不足以进行普通 phono3py HDF5 输运。
力缓存复用与外部求解器 checkpoint 重启是不同机制，缓存格式与校验见 :doc:`architecture`。

截断与拟合
------------------------

原生 phono3py 的 ``cutoff-fc3`` 是位移对距离，单位 Å；
Thirdorder 的负整数表示邻居壳层，正值单位 nm；
HiPhive ``cutoffs`` 是各阶拟合半径，单位 Å；
Fourthorder 的负 ``cutoff-fc4`` 表示邻居壳层。
这些截断不能互相照搬，要与超胞共同检查收敛。

替换计算器
------------------------

VASP 可从 ``examples/vasp-rta-3ph.yaml`` 修改启动命令与 POTCAR 来源。
``potcar-path`` 可指向合并文件或库目录，程序按结构元素顺序组装。
生成位移力前应检查分阶段电子参数并验证收敛。

MACE 后端随标准安装提供，使用时需配置匹配的 checkpoint。
只有环境支持时才选 ``device: cuda``。其他 ASE 计算器用 ``factory`` 与 ``kwargs`` 接入。
模型、工厂和缓存依赖不明的后端，只在同一计算 / Slurm 任务范围内复用力，见 :doc:`architecture`。

薄膜需要物理有效厚度，纳米线需要物理横截面积；不能把真空尺寸直接当作这些数值。
相应网格和超胞也需遵循周期方向，绘图归一化不会改写原始 HDF5。

Slurm 执行
------------------------

可把完整工作流放入外层 batch，也可按已支持路线设置 ``parallel.force-constant`` 或 ``parallel.kappa``。
先用 ``submit: false`` 检查生成脚本，再按需要提交。
``info`` 不执行准备工作；实际执行阶段即使不提交，也可能先发现 q 点或执行此前的本地阶段。
QHA / SSCHA 温度循环不能嵌套异步输运。

GPU 位移力还需 ``nep-force-backend: gpu-batch``、Calorine >= 3.5 与 ``nep`` 程序，并为 worker 申请 GPU。
并发太多短任务不一定更快，调度和 I/O 开销需要在目标设备实测。

只绘制已有 FC2
------------------------

结果目录包含 ``fc2.hdf5`` 及其匹配的 ``phono3py_disp.yaml``、``phonopy.yaml`` 或 ``phonopy_disp.yaml``。
精简输入例子：

.. code-block:: yaml

   output:
     result-dir: calculations/my-material/results
   kappa:
     mesh: [21, 21, 21]
     temps: [0, 1000, 100]
   plot:
     layout: both
     path: seekpath

.. code-block:: bash

   nepkappa info plot.yaml --for plot
   nepkappa plot plot.yaml

这会生成色散、DOS、体积热容与群速度的单图和 2×2 组合图，不需要势函数或 FC3。
这里的温度网格服务于谐近似绘图，允许 0 K，与热输运温度必须正的要求不同。
同名图像可能被替换。仅在 YAML 写 ``plot`` 不会自动生成图像。

存在热导率 HDF5 时，网格必须与保存数据匹配；线宽支持散射率，``mode_kappa`` 支持累积热导率。
``report`` 只总结已有产物，不补算缺失物理量。论文定制绘图脚本与内置标准绘图分开维护。

.. _change-the-transport-approximation:
.. _three-and-four-phonon-conductivity:
.. _qha-and-qha-volume-sscha:
.. _time-dependent-bte:

按物理问题选择章节
------------------------

* :doc:`harmonic_temperature`：谐近似、QHA、SSCHA 与 QHA+SSCHA。
* :doc:`thermal_transport`：LBTE、同后端 3ph / 4ph 比较、tmux 与 Wigner Slurm。
* :doc:`tdbte`：非平衡激发、数值审计和弛豫图。
* :doc:`tdbte_spatial`：TTG / 脉冲、空间漂移与能量连续性。

温度相关热导率还需相应 ``run-transport`` 开关和相容的 QHA / SSCHA / 输运温度网格，见 :doc:`input_files`。
