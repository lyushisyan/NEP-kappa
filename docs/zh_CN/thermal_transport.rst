声子散射与晶格热导率
================================

本章比较三声子 LBTE、三加四声子输运，以及 Wigner 的粒子 / 相干贡献。
复用数据时，力常数、晶胞坐标和晶格平移标记必须遵循同一约定。

准备已有力常数
------------------------

延续 :doc:`starting` 的结果：

.. code-block:: bash

   mkdir -p ~/work/3C-SiC/transport
   cd ~/work/3C-SiC/transport
   ls -lh ../results/fc2.hdf5 ../results/fc3.hdf5 ../results/phono3py_disp.yaml
   cp ../results/POSCAR_relaxed POSCAR
   cp ../nep_3C-SiC.txt .

LBTE 示例复用已有 FC2/FC3；四声子示例重新生成一组匹配的 FC2/FC3/FC4，供后续 Wigner 计算复用。

三声子 LBTE
------------------------

RTA 分别处理各模式对平衡分布的偏离；LBTE 保留碰撞矩阵中模式之间的线性耦合。
下载 :download:`lbte.yaml <../source/inputs/lbte.yaml>`：

.. literalinclude:: ../source/inputs/lbte.yaml
   :language: yaml

.. code-block:: bash

   nepkappa info lbte.yaml --for kappa
   nepkappa kappa lbte.yaml
   tail -n 20 results/lbte/run.log
   ls -lh results/lbte/kappa-m151515.hdf5
   nepkappa info lbte.yaml --for plot
   nepkappa plot lbte.yaml

命令会将匹配的来源输入复制到结果目录供 phono3py 与绘图使用，不重新计算原子力。
LBTE 的碰撞矩阵通常比 RTA 消耗更多内存，示例网格不是收敛结论。

普通 LBTE 结果中的 ``kappa`` 与 ``kappa_RTA`` 对应同次计算的 LBTE 和 RTA 热导率，
绘图分别用实线、虚线。``plot.kappa: x`` 选择 xx 分量，``all`` 选择三个对角分量的平均。
标准六面板显示四种谐近似性质、指定温度的散射率和热导率温度曲线。

``tau: total`` 的曲线来自保存的 ph-ph 与同位素通道及其和，这个显示总和不含边界散射。
本例显示下限为 1e-8 ps⁻¹。N/U 图需要单独保存的 N、U 线宽；累积热导率需要模式分辨数据。

四声子输运与同后端比较
------------------------------

先准备 Thirdorder、Fourthorder 和 FourPhonon CPU 程序：

.. code-block:: bash

   command -v thirdorder_vasp.py
   command -v Fourthorder_vasp.py
   command -v ShengBTE_cpu

下载 :download:`fourphonon.yaml <../source/inputs/fourphonon.yaml>`：

.. literalinclude:: ../source/inputs/fourphonon.yaml
   :language: yaml

Thirdorder 生成 FC3，Fourthorder 生成 FC4，原子力由 NEP 提供；脚本文件名中含 ``vasp``
不意味着必须使用 VASP 力。负截断值选择邻居壳层。
``format: both`` 保留 ShengBTE 文本并导出用于分析的 HDF5。
示例使用四个 MPI 进程、每进程一个 OpenMP 线程，需按实际资源调整。

``compare-3ph: true`` 顺序执行同后端的纯三声子参考与三加四声子计算，保持 FC2/FC3、结构、网格、
温度和同位素等设置一致。这里三声子用 LBTE，四声子部分用 RTA。

.. _transport-tmux:

使用 tmux 保持远程会话
------------------------------

.. code-block:: bash

   nepkappa info fourphonon.yaml --for run
   tmux new -s fourphonon
   conda activate nepkappa
   cd ~/work/3C-SiC/transport
   nepkappa run fourphonon.yaml

依次按 ``Ctrl-b``、松开，再按 ``d`` 可离开会话，SSH 断开后服务器中的进程仍可运行。

.. code-block:: bash

   tmux ls
   tmux attach -t fourphonon

``Ctrl-b`` 后按 ``[`` 可滚动，``q`` 退出滚动。
tmux 只维持终端会话，不替代调度器，也不会延长已有的资源分配时限。
详细输出使用 ``-v``；每个 FourPhonon 工作目录另有 ``fourphonon.log``。

结果与比较图
------------------------

.. code-block:: bash

   tail -n 20 results/fourphonon/run.log
   ls results/fourphonon/fourphonon/threephonon
   ls results/fourphonon/fourphonon
   nepkappa plot fourphonon.yaml

参考目录保存 ``kappa3-rta.dat``、``kappa3-iterative.dat``；组合目录保存 ``kappa4-rta.dat``、
``kappa4-iterative.dat``。本例后者是混合的 3ph LBTE + 4ph RTA 解。
两种 method 都为 ``lbte`` 时，它才是三加四声子的迭代解。

热导率图比较 3ph RTA、3ph LBTE、3ph+4ph RTA 与混合解。
``scattering_rate_3ph_4ph.png`` 单独标记两种散射率，面板中的曲线不额外加入同位素，
尽管热导率计算可以包含同位素。加四个谐近似面板和热导率面板，组合图共六面板。

四声子 Wigner 输运
------------------------

Wigner 将热导率分为粒子项和相干项，二者之和为总热导率。
此例需要 Wigner_Park，三、四声子都使用 RTA，只计算组合结果。
在独立工作目录的 ``inputs/`` 中放置同一来源的：

.. code-block:: text

   POSCAR_relaxed
   FORCE_CONSTANTS_2ND
   FORCE_CONSTANTS_3RD
   FORCE_CONSTANTS_4TH
   fc2.hdf5
   fc3.hdf5
   phono3py_disp.yaml

下载 :download:`wigner-4ph.yaml <../source/inputs/wigner-4ph.yaml>`，替换其中的可执行文件占位路径：

.. literalinclude:: ../source/inputs/wigner-4ph.yaml
   :language: yaml

.. code-block:: bash

   nepkappa info wigner-4ph.yaml --for kappa

这是复用输入，不需要计算器。自动 CONTROL 关闭 NAC；需要长程非解析修正时，应提供带有经过验证
Born 电荷和介电数据的自洽 CONTROL。

.. _transport-slurm:

使用 Slurm 分配资源
------------------------

保存为 ``submit-wigner-4ph.sh``，将分区、环境路径和资源替换为所在集群设置。
下面的四个任务需与 YAML 中 ``parallel.kappa.mpi-processes: 4`` 匹配。

.. code-block:: bash

   #!/bin/bash -l
   #SBATCH --job-name=SiC-Wigner-4ph
   #SBATCH --partition=YOUR_PARTITION
   #SBATCH --nodes=1
   #SBATCH --ntasks=4
   #SBATCH --cpus-per-task=1
   #SBATCH --time=48:00:00
   #SBATCH --output=logs/wigner-4ph-%j.out
   #SBATCH --error=logs/wigner-4ph-%j.err

   set -eo pipefail
   cd "$SLURM_SUBMIT_DIR"
   source /path/to/miniconda3/etc/profile.d/conda.sh
   conda activate nepkappa
   export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
   nepkappa kappa wigner-4ph.yaml

按站点要求加载编译器与 MPI 模块。外层只运行一次 NEP-kappa，MPI 后端由 YAML 设置启动。

.. code-block:: bash

   mkdir -p logs
   bash -n submit-wigner-4ph.sh
   sbatch submit-wigner-4ph.sh
   squeue -u "$USER"

``bash -n`` 仅检查脚本语法；``sbatch`` 才会提交任务。
结束后用真实作业 ID 替换 ``JOB_ID``：

.. code-block:: bash

   sacct -j JOB_ID --format=JobID,State,Elapsed,ExitCode
   tail -n 20 results/wigner-4ph/run.log

需同时核对 ``COMPLETED``、退出码 ``0:0`` 与计算日志；从队列消失本身不代表成功。

Wigner 绘图
------------------------

标准绘图还需要同一输入集的谐近似 HDF5 和元数据：

.. code-block:: bash

   cp inputs/fc2.hdf5 inputs/phono3py_disp.yaml results/wigner-4ph/
   nepkappa info wigner-4ph.yaml --for plot
   nepkappa plot wigner-4ph.yaml

粒子、相干和总热导率分别导出为 ``kappa4-rta.dat``、``kappa4-rta-coherence.dat``、
``kappa4-rta-total.dat``。接近零的相干曲线应先检查数值，再判断是否可以描述为零。
