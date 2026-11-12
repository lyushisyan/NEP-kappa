输入文件与参数
========================

完整输入与单阶段输入
--------------------------

``nepkappa run`` 要求完整路线的六个部分。声子 / 输运路线使用 ``structure``、``calculator``、
``force-constant``、``kappa``、``plot``、``output``；非平衡动力学用 ``tdbte`` 替代 ``kappa``。
未配置选项的部分可写成 ``{}``。两种路线不能在同一输入中同时出现 ``kappa`` 和 ``tdbte``。
可选的 ``parallel`` 是额外的资源设置部分。

.. code-block:: yaml

   structure:
     poscar: examples/structures/3C-SiC/POSCAR_primitive
     relaxation: true
   calculator:
     name: nep
     model: examples/potentials/3C-SiC/nep_3C-SiC.txt
   force-constant:
     dim-fc2: [3, 3, 3]
     dim-fc3: [3, 3, 3]
   kappa:
     engine: phono3py
     mesh: [31, 31, 31]
     temps: [100, 1000, 50]
     method: rta
   plot:
     layout: both
   output:
     result-dir: calculations/example-runs/3c-sic-nep-rta

直接命令只执行对应的计算阶段，允许使用精简输入，也接受完整输入。

.. list-table:: 单阶段所需内容
   :header-rows: 1
   :widths: 24 76

   * - 命令
     - 必需内容和说明
   * - ``relax``
     - ``structure``、``calculator``。
   * - ``fc2`` / ``fc2fc3``
     - ``structure``、``calculator``、``force-constant``。
   * - ``fc2fc3fc4``
     - 同上；依次生成 FC2、FC3、FC4，不执行结构弛豫或热导率计算。
   * - 普通 phono3py ``kappa``
     - ``kappa``；可用 ``kappa.force-constants`` 选择已有力常数。
   * - ``tdbte``
     - ``tdbte``；指定力常数目录与网格，或指定已有 ``kernel``。
   * - ``plot``
     - 可只写 ``output`` 使用默认绘图，或加入 ``plot``；按需增加网格、温度及几何参数。
   * - ``report``
     - ``output.result-dir``，也可直接把结果目录作为命令参数。

除报告 YAML 外，``output`` 可省略，默认结果目录为启动目录下的 ``result``。
建议显式指定输出路径。``info --for <命令>`` 检查对应的阶段；不带 ``--for`` 时检查完整 ``run`` 输入。
精简的复用输入不能交给 ``run``。内部 FC4、SSCHA、bubble 和耦合步骤由输入开关选择，没有独立公开命令。

``fc2fc3fc4`` 组合执行已有的 FC2/FC3 和 FC4 阶段，FC4 通过 Fourthorder
生成 ``FORCE_CONSTANTS_4TH``，FC2/FC3 的后端和输出格式遵循输入设置。
可使用精简力常数输入或普通四声子完整输入；此命令不要求额外启用
``force-constant.four-phonon``。QHA、SSCHA 和 TD-BTE 配置应使用 ``run``。
FC2/FC3 的原子力任务交给 Slurm 时，命令先暂停，生成的收集脚本会在收集完成后继续生成 FC4。

.. code-block:: bash

   nepkappa info input.yaml --for fc2fc3fc4
   nepkappa fc2fc3fc4 input.yaml

.. important::

   普通相对路径以启动目录为准，不以 YAML 文件所在目录为准。
   ``info`` 检查输入结构与参数关系，不会验证模型精度、外部程序可运行性或数值收敛。

structure：结构与归一化
--------------------------

``poscar`` 指向结构。``relaxation: true`` 启用弛豫，也可以配置 VASP 的分阶段弛豫：

.. code-block:: yaml

   structure:
     poscar: POSCAR
     relaxation:
       enabled: true
       workdir: vasp-relax
       stages:
         coarse: {ediffg: -0.05, prec: Normal}
         fine: {ediffg: -0.01, prec: Accurate}

启用弛豫后，单独运行 ``fc2`` / ``fc2fc3`` / ``fc2fc3fc4`` 之前需先运行 ``relax``，因为力常数阶段读取 ``POSCAR_relaxed``。
完整 ``run`` 会安排顺序。QHA 直接读取指定结构，因此它的输入要求 ``relaxation: false``。

``dimensionality`` 可为 1、2 或 3。薄膜需指定物理有效厚度 ``effective-thickness``（Å），
可选 ``vacuum-axis``；纳米线需指定 ``effective-area``（Å²），可选 ``periodic-axis``。
这些因子修正绘图中的热容和热导率，不改写原始 phono3py HDF5。
当前 QHA、SSCHA 和 FourPhonon 只接受 ``dimensionality: 3``，单独设置几何归一化不能消除这一限制。

calculator：能量和原子力
--------------------------

.. list-table:: 计算器选项
   :header-rows: 1
   :widths: 30 70

   * - 选项
     - 含义
   * - ``name``
     - 默认 ``nep``；可选 ``vasp``、``mace``、``ase`` 或已安装的插件名。
   * - ``model``
     - NEP 势文件，或 MACE checkpoint / 基础模型工厂认可的标识。
   * - ``nep-model``
     - 保留的 NEP 旧键，不能与 ``model`` 同时出现，即使值相同。
   * - ``vasp-command`` / ``vasp-path``
     - VASP 启动命令 / 可执行文件，前者优先。
   * - ``potcar-path``
     - 合并的 POTCAR 文件或势库目录，必须在执行主机可访问。
   * - ``foundation``、``device``、``dtype``
     - MACE 工厂、设备和精度设置。
   * - ``factory``、``kwargs``、``model-files``
     - 通用 ASE 工厂路径、参数和显式纳入身份记录的模型文件。

GPU 位移力使用：

.. code-block:: yaml

   calculator:
     name: nep
     model: /path/to/nep.txt
     nep-force-backend: gpu-batch
     nep-batch-size: 128
     nep-command: /path/to/nep

``nep-force-backend`` 默认 ``cpu``。``nep-command`` 必须是一份可执行文件，不是 shell 命令或 ``gpumd``。
批处理保留结构和原子顺序，拒绝形状错误或非有限力；失败时不切换 CPU。
CPU / GPU 缓存身份不同，GPU 记录模型、程序及批大小等信息。
GPU 目前只用于原生 FC2/FC3 与 HiPhive 位移力，弛豫仍为 CPU；不适用于 QHA、SSCHA 和 Thirdorder / Fourthorder 路线。
安装与节点检查见 :doc:`installation`。

force-constant：力常数与温度修正
------------------------------------

.. list-table:: 力常数选项与单位
   :header-rows: 1
   :widths: 30 70

   * - 选项
     - 含义
   * - ``dim-fc2`` / ``dim-fc3`` / ``dim-fc4``
     - 分别指定二、三、四阶力常数超胞。
   * - ``fc3-backend``
     - 默认原生 phono3py；``thirdorder`` 使用 Thirdorder 生成 FC3。
   * - ``cutoff-fc3``
     - 原生 phono3py：位移原子对距离截断，单位 Å；Thirdorder：负整数为邻居壳层，正值单位 nm。
   * - ``cutoff-fc4``
     - 传给 Fourthorder；负值选择邻居壳层，需独立检查收敛。
   * - ``use-hiphive``
     - 启用 HiPhive 拟合，配合 ``n-structures``、``rattle-std``、``cutoffs`` 和 ``min-dist``。
   * - ``format``
     - ``phono3py``、``shengbte`` 或 ``both``；phono3py 需要 HDF5，FourPhonon 需要文本 FC2/FC3/FC4。
   * - ``thirdorder-command`` / ``fourthorder-command``
     - 外部位移生成 / 收集程序的启动命令。
   * - ``qha`` / ``sscha`` / ``four-phonon``
     - 可为布尔值，也可为含 ``enabled`` 和具体参数的映射。

原生 FC3 的截断限制位移原子对，不等于把张量在同一距离处直接截断。
HiPhive 的 ``cutoffs`` 是按阶设置的拟合截断，单位 Å。原生 FC2 没有独立的同类实空间截断，需检查超胞。

QHA 与 SSCHA 的配置关系
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   force-constant:
     dim-fc2: [3, 3, 3]
     dim-fc3: [3, 3, 3]
     qha:
       enabled: true
       volume-ratios: [0.94, 0.97, 1.00, 1.03, 1.06]
       temps: [0, 1000, 10]
     sscha:
       enabled: true
       run-transport: true
       temps: [300, 900, 300]
       snapshots: 200
       iterations: 4
   kappa:
     engine: phono3py
     temps: [300, 900, 300]
     method: rta

QHA 各向同性缩放体积并拟合状态方程。常用参数还有 ``dim-fc2``、``mesh``、``eos``、``pressure``、
``relax-internal`` 以及 VASP 弛豫 / 静态计算参数。
仅 QHA 的 ``run-transport: true`` 路线在各温度的 QHA 平衡体积上重新生成 FC2/FC3，
目前只支持 phono3py RTA；它不包含显式非谐频率重整化。

``qha.equilibrium-phonons: true`` 会在体积拟合后，为各输出温度的平衡体积额外计算 FC2，
写入 ``qha-phonons/T.../`` 和 ``qha-phonons.yaml``。平衡体积超出扫描范围时拒绝计算。
随后 ``plot`` 读取这些结果生成自由能、平衡体积、体积热膨胀及色散 / DOS 四面板。
缺少这些平衡体积 FC2 的旧结果，需要重新执行对应计算，绘图本身不生成力。

SSCHA 在每个温度重新拟合辅助谐近似 FC2。常用参数为 ``snapshots``、``iterations``、``transient``、
``sscha-mesh``、``initial-fc2``。QHA+SSCHA 在 QHA 给定的平衡体积上进行固定体积 SSCHA，
不会用 SSCHA 自由能重新优化体积。该路线中的 FC3/FC4 在对应体积重新生成，但不作显式温度重整化。

QHA / SSCHA 默认 ``run-transport: false``。同时启用时，任一 ``run-transport: true`` 可选择组合输运，
显式冲突的开关会报错。耦合温度循环必须同步完成各温度的输运，拒绝嵌套 Slurm 输运，
即使 ``submit: false`` 也不例外。可把整个流程放在外层 Slurm 分配中，内部使用 ``backend: none``。

``sscha.bubble`` 保存由输入 FC3 得到的独立对角 on-shell 频移，不更新导出的 FC2 或热导率。
详细输入见 :doc:`harmonic_temperature`。

kappa：散射与热导率
--------------------------

``engine: phono3py`` 使用 ``method: rta`` 或 ``lbte``；``isotope`` 启用同位素散射，
``bfmp`` 的单位为微米，``wigner`` 启用 SMM19 相干贡献。
FourPhonon 使用 ``engine: fourphonon``，并分别设置 ``method-3ph`` 与 ``method-4ph``：
支持 ``rta/rta``、``lbte/rta``、``lbte/lbte``，不支持 ``rta/lbte``。
四声子 Wigner 需要 ``rta/rta`` 和 Wigner_Park。

``four-phonon.compare-3ph: true`` 用同一 FourPhonon 程序顺序运行纯三声子与三加四声子两组结果，
保持结构、FC2/FC3、网格、温度、展宽和同位素设置一致。
参考结果在 ``fourphonon/threephonon/``，组合结果在 ``fourphonon/``。
该比较要求自动生成 CONTROL，不能与自定义 CONTROL、Wigner 或 only-harmonic 组合。

FC2、Thirdorder FC3、Fourthorder FC4 和自动 CONTROL 共用规范化晶胞坐标：
周期位置折叠到 [0, 1)，距边界 1e-10 内的浮点噪声归零，原子顺序与晶格保持不变。
只修改 CONTROL 不能转换旧力常数中的晶格平移标记；不一致的旧数据应重新生成或按原始元数据显式重映射。
自动 CONTROL 默认关闭 NAC；极性材料需要相应的 Born 电荷和介电数据时，应使用经过验证的一致输入。

温度网格
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``qha.temps``：``[Tmin, Tmax, step]``，要求 ``0 <= Tmin < Tmax``。
* ``sscha.temps``：同样是范围写法，允许 ``Tmin == Tmax`` 表示单点，步长仍须为正。
* ``kappa.temps``：``[T]`` 或 ``[Tmin, Tmax, step]``，温度必须大于零，不接受任意不等间隔列表。

步长须为正并恰好覆盖终点。SSCHA 范围必须位于 QHA 范围内，允许对 QHA 平衡体积插值。
启用 SSCHA 输运时，每个热导率温度必须是实际 SSCHA 采样点；它们可以是采样点的子集。
例如 QHA ``[0,1000,10]``、SSCHA ``[300,900,100]`` 与输运 ``[300,900,300]`` 相容。
耦合输运必须显式写出 ``kappa.temps``；禁用输运时，它不会限制 QHA / SSCHA 温度范围。

复用已有力常数
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   kappa:
     force-constants: calculations/3C-SiC/fc
     method: rta
     mesh: [5, 5, 5]
     temps: [300]
   output:
     result-dir: calculations/3C-SiC/kappa-reuse

使用 ``nepkappa kappa kappa.yaml``。来源目录必须包含配套的 ``phono3py_disp.yaml``、``fc2.hdf5``、``fc3.hdf5``。
来源与输出不同时，普通 phono3py 路线复制三份文件及可选 BORN；冲突输入会被拒绝，源文件不改动。
该来源选项用于普通直接 phono3py 输运，不用于 FourPhonon、耦合路线或从头执行的 ``run``。

自定义命令
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``kappa.command`` 替代整个 phono3py 命令，不会再自动附加其他计算选项。
程序可读取命令里的 ``--mesh``、``--ts`` / ``--temperature``，或完整的 ``--tmin/--tmax/--tstep``。
命令行 ``--ts`` 是实际温度点列表，与 YAML 的三数范围约定不同；多点须递增且等间隔。
若同时写 YAML 网格 / 温度，二者必须一致。隐藏在外部脚本或配置文件里的采样设置不会自动推断。
命令在结果目录中执行，不经过 shell，因此管道和重定向不会被解释。

tdbte：非平衡动力学
--------------------------

从头计算用 :download:`完整示例 <../../examples/tdbte.yaml>` 和 ``run``；
已有力常数用 ``tdbte.force-constants`` 与 ``tdbte.mesh``，再执行 ``nepkappa tdbte``。
已有碰撞核则只指定 ``kernel``，不能同时指定网格或力常数来源。
三维碰撞网格三个分量均须至少为 2；实验一维路线只允许一个分量大于 1，二维网格不支持。

``branches`` 是从零开始、按各 q 点频率排序的分支编号。
``excitation`` 表示所选占据数乘以 ``1 + excitation``，不是温升或激光能流密度。
``temperature`` 为初始 Bose 参考温度；``duration-ps``、``max-step-ps``、``samples`` 分别控制时长、
最大步长和保存时刻数量。默认值及数值模型见 :doc:`tdbte`。
可选 ``spatial`` 打开周期一维实空间输运，此时 q 网格仍为三维，详见 :doc:`tdbte_spatial`。
QHA、SSCHA、FC4 和 Wigner 设置不会自动传入 TD-BTE。

plot：绘图
--------------------------

.. list-table:: 常用绘图选项
   :header-rows: 1

   * - 选项
     - 含义
   * - ``layout``
     - ``separate``、``combined`` 或 ``both``。
   * - ``path`` / ``path-points`` / ``path-segments``
     - 自动或自定义色散路径。
   * - ``temperature``
     - 标准输运图中散射率所用温度；缺少精确匹配时取最近保存点并提示。
   * - ``kappa``
     - ``x``、``y``、``z`` 或 ``all``；``all`` 表示三个对角分量的平均。
   * - ``tau``
     - ``phph``、``isotope``、``total``、``normal``、``umklapp``、``nu`` 或 ``all``。
   * - ``scattering-rate-min``
     - 仅显示层面的下限，单位 ps⁻¹；默认 0，等于下限的点保留，不改写线宽或热导率。
   * - ``dpi`` / ``tdbte-times``
     - PNG 分辨率 / 均匀 TD-BTE 色散快照的三个时间。

FC2 和匹配结构元数据足以绘制色散、DOS、体积热容和群速度。
有热导率、线宽和模式热导率时，再分别增加对应面板，裸 FC2 数组不足以建立色散。
``nu`` 要求 ``gamma_N``、``gamma_U`` 均存在。
线宽转换为散射率使用 ``4*pi*gamma``；无温度轴的同位素线宽会广播到各温度。
``total`` 比较 ph-ph、同位素和二者之和，图中的这个和不含边界散射。
显式要求同位素但结果缺少 ``gamma_isotope`` 时，不会把不完整的和当作总散射率。
FourPhonon 的 3ph / 4ph 散射率分别标记，不在面板中求和或额外加入同位素。

普通非 Wigner LBTE 数据含 ``kappa_RTA`` 时，同图比较 LBTE 实线与 RTA 虚线，应用相同几何修正。
Wigner 读取器验证 ``kappa_intra_exact`` / ``kappa_intra`` 与 ``kappa_inter`` 的存储约定，
不会用 RTA 粒子解代替 LBTE 粒子解。

禁用输运时，SSCHA 总览叠加所有采样温度的谱和群速度；热容使用各温度 FC2 的独立彩色点。
灰色谐近似热容参考使用固定体积 FC2，QHA+SSCHA 则优先使用体积比为 1 的 QHA 扫描点。
该热容不是完整 SSCHA 自由能对温度的导数。绘图始终不会补算原子力。

output 与 parallel
--------------------------

``output.result-dir`` 指定结果目录，``progress`` 控制活动进度显示。
为不同材料和参数组使用独立目录。可选资源配置例如：

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

``submit: false`` 写出脚本而不调用 ``sbatch``，但执行阶段仍可能进行准备工作；例如输运需要本地发现 q 点。
phono3py 的 RTA / LBTE 都支持 Slurm，worker 以轮转方式分配不可约 q 点。
``jobs`` 控制子集数，``max-concurrent`` 限制并行 worker 数；各任务共享结果文件系统。
脚本从 ``SLURM_CPUS_PER_TASK`` 设置 Rayon 与 OpenMP 线程；BLAS 线程需单独管理。

RTA 使用 ``--write-gamma/--read-gamma``，每次准备新建 ``rta-slurm/run-*`` 工作区，避免旧整网格结果
遮蔽 worker 输出；LBTE 使用 ``--write-pp/--read-pp``。
只有所有 worker 成功后才收集并发布 ``kappa-m*.hdf5``，随后标记完成。
准备或提交并不等于科学计算已完成；收集后再执行 ``plot``。
GPU 资源可通过 ``extra-sbatch: ["--gres=gpu:1"]`` 交给目标集群，但实际语法与环境需按站点配置。
