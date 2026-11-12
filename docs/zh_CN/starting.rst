快速开始：三声子 RTA 与 Wigner 输运
============================================

本例使用随仓库提供的 3C-SiC NEP 模型，依次完成结构弛豫、FC2/FC3 生成和
三声子 RTA / SMM19 Wigner 输运。输入中的数值是流程演示设置，需要针对研究对象检查收敛。

准备工作目录
------------------------

在已安装项目的仓库根目录运行：

.. code-block:: bash

   conda activate nepkappa
   mkdir -p ~/work/3C-SiC
   cp examples/structures/3C-SiC/POSCAR_primitive ~/work/3C-SiC/POSCAR
   cp examples/potentials/3C-SiC/nep_3C-SiC.txt ~/work/3C-SiC/
   cd ~/work/3C-SiC

普通相对路径以启动命令时的目录为准。后续章节在这个目录下使用 ``phonons/``、``transport/``
和 ``tdbte/`` 分开保存结果。

准备并检查输入
------------------------

下载 :download:`quickstart.yaml <../source/inputs/quickstart.yaml>`，在工作目录保存为 ``input.yaml``：

.. literalinclude:: ../source/inputs/quickstart.yaml
   :language: yaml

本例采用 3×3×3 力常数超胞、31×31×31 q 网格，在 100–1000 K 范围内每隔 50 K 计算一次。
``wigner: true`` 启用粒子与相干贡献，``isotope: true`` 包括同位素散射。

.. code-block:: bash

   nepkappa info input.yaml --for run
   nepkappa run input.yaml

第一条命令仅检查输入，第二条才开始计算。加 ``-v`` 可查看详细终端输出：
``nepkappa run input.yaml -v``。两种输出模式都将完整诊断写入 ``results/run.log``。

长时间工作站计算可使用 :ref:`tmux <transport-tmux>`；集群计算应在合适的
:ref:`Slurm 资源分配 <transport-slurm>` 内运行，并匹配 MPI / 线程数。

检查结果
------------------------

.. code-block:: bash

   tail -n 20 results/run.log
   ls -lh results/fc2.hdf5 results/fc3.hdf5 results/phono3py_disp.yaml
   ls -lh results/kappa-m313131.hdf5

``POSCAR_relaxed`` 是弛豫后的结构。HDF5 力常数和 ``phono3py_disp.yaml`` 必须配套使用。
``kappa-m313131.hdf5`` 保存指定网格上的热导率及模式数据；网格变化时文件名也变化。

绘图与报告
------------------------

.. code-block:: bash

   nepkappa info input.yaml --for plot
   nepkappa plot input.yaml
   nepkappa report input.yaml
   ls results/plots

``layout: both`` 生成各单图和 ``combined.png``。标准六面板包括色散、DOS、体积定容热容、群速度、
散射率和热导率。Wigner 热导率区分粒子项、相干项与总和。
``plot.temperature`` 选择散射率所用温度，热导率图保留全部计算温度。

报告保存为 ``results/report.md`` 和 ``results/report.yaml``。
绘图与报告读取已有数据，不重新计算力或热导率。完整输入目录见 :doc:`examples`。
