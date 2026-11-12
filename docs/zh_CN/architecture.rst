代码架构与计算阶段
================================

两条科学路线
------------------------

.. code-block:: text

   结构 + 计算器 → 位移结构 → 原子力 → FC2 / FC3 / FC4
                                           │
                  ┌────────────────────────┴─────────────────┐
                  │                                          │
            声子与热输运                                非线性 TD-BTE
         FC2 → 谐近似性质                            FC2 + FC3 → 碰撞核
         FC2 + FC3 → RTA / LBTE / Wigner              → n(q, branch, t)
         可选 FC4 → 四声子输运                        可选周期实空间网格
         可选 QHA / SSCHA → 温度相关算例               → n(x, q, branch, t)

两条路线共用力常数准备，但热导率不是 TD-BTE 输入。
QHA / SSCHA 拥有自己的体积、温度和迭代循环，其设置不会自动传入 TD-BTE。

模块职责
------------------------

.. list-table:: 从入口到结果
   :header-rows: 1
   :widths: 37 63

   * - 模块
     - 职责
   * - ``cli.py``、``command_registry.py``
     - 公开命令、终端 / 日志、配置选择与命令级溯源。
   * - ``config.py``、``config_schema.py``
     - 严格 YAML、键规范和 argparse 类型转换 / 默认值。
   * - ``config_plan.py``、``config_validation.py``
     - 科学开关推导阶段，并检查所选操作的前提。
   * - ``config_models.py``
     - 可变主配置及按需创建的类型化分区快照。
   * - ``application.py``、``stages/``
     - 执行阶段、计时、处理延期、结构弛豫和力常数生成。
   * - ``workflow.py``
     - 共享结构、计算器、力作业和导出上下文，不负责执行完整计划。
   * - ``calculators.py``、``adapters/vasp.py``
     - ASE 兼容计算器构造，以及 VASP 专用文件 / 执行逻辑。
   * - ``transport.py``、``fourphonon.py``、``qha*.py``、``sscha.py``
     - 输运接口和对应科学循环。
   * - ``scheduler.py``、``slurm.py``、``run_state.py``
     - 调度提交、脚本 / manifest 与延期作业状态。
   * - ``transport_data.py``、``plot*.py``、``report.py``
     - 数据解释、单位 / 几何归一化、绘图和报告。

耦合温度流程通过显式注入的 ``runner_factory`` 创建各算例阶段执行器，不反向导入应用层。
``identity.py`` 提供缓存和溯源共用的哈希工具；``temperatures.py`` 提供温度网格和目录名称；
``command_sampling.py`` 负责 phono3py 命令采样与常规结果名称。公共工具不依赖上层科学工作流。

配置和执行契约
------------------------

输入路径为 YAML → 科学开关展开 → argparse 参数 → ``WorkflowConfig`` → 操作校验 → 有序阶段。
Argparse 继续负责标量转换、默认值和 choices，避免维护第二套类型转换规则。
``WorkflowConfig`` 是主状态，``WorkflowSections`` 按需构造；算例复制配置后再更改体积、温度和路径。

``run`` 要求对应路线的六部分；直接命令只要求自身阶段内容。
``complete``、``prepared``、``submitted`` 分别表示完成、准备成功和提交成功。
准备或提交延期任务后，本地计划停止，不能继续声称科学结果已完成。
QHA / SSCHA 温度循环要求同步输运，并在计算前拒绝不支持的嵌套 Slurm 设置。

力缓存与结果身份
------------------------

``force_identity.py`` 定义 schema-3 力输入与显式依赖清单：结构、计算器参数、模型 / POTCAR 内容、
力求值代码和已安装提供者的 Python / 原生文件身份。
本仓库 Python 用规范化语法摘要排除注释和文档字符串；共享上下文只选择明确的力相关方法。
新增影响力的函数或后端时，应同步审查依赖清单。

``ForceArtifactStore`` 纳入精确 ASE 数组、晶胞和周期性。
Slurm 准备保存带哈希的 ``atoms.json``，worker 恢复原状态，避免 POSCAR 往返丢失附加数组或精度。
无法确定隐藏输入的任意 info / constraint 对象不允许持久复用。

当前支持跨独立运行复用的后端包括 NEP、ASE LJ / EMT / Morse 工厂和直接识别的原生 VASP 命令。
不透明插件、MACE 工厂及 VASP 包装 / launcher 命令限定于同一计算或 Slurm 管线内复用，
worker 与收集阶段共用 scope token，独立运行则产生新 token。

安装提供者的摘要在进程内缓存，所以运行期间必须保持代码和安装文件不可变；
模型输入在后续指纹调用中重新核对。这是保守的依赖追踪，不是数值等价性的数学证明。

读取缓存时检查输入 JSON、输入摘要、输出 SHA-256、形状与有限性；``forces.npy`` 和 Thirdorder / Fourthorder 的 XML 都有内容校验。
发布先撤销有效标记，再替换输出、写记录，最后发布标记。
数组在验证后才替换，XML 后端可能在求值中覆盖旧文件，因此调用前就撤销旧标记。
旧 schema 直接视为缓存未命中，不会把旧记录自动认证成新身份。
该流程保证单写者的进程中断恢复，不提供并发写锁或断电持久性协议。

``provenance.py`` 独立记录更广泛的软件版本（含 SciPy）、输入、执行历史和产物哈希。
调度作业结束不会回写原始命令的历史溯源记录。

TD-BTE 子系统
------------------------

``tdbte_builder.py`` 准备原生模式、顶角和速度；``tdbte_quadrature.py`` 积分能量壳；
``tdbte_storage.py`` 验证并流式读取事件块；``tdbte_accumulate.py`` 提供 NumPy / Numba 归约。
``tdbte.py`` 负责均匀演化和公共产物加载。

``tdbte_streaming.py`` 负责周期守恒漂移，``tdbte_observables.py`` 负责温度、能量、热流与熵，
``tdbte_spatial.py`` 组合非线性碰撞与漂移、审计并保存 HDF5。
事件分块不消除单元 × 模式状态内存，也不替代物理速率与收敛验收。

维护原则
------------------------

保持模块围绕具体职责，只有多个实际调用者需要时才提取公共工具。
保留科学数据的单位、模式顺序和拒绝策略；不要为减少代码行数混淆数值算法与调度逻辑。
测试与文档约定见 :doc:`development`。
