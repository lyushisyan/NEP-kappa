输入助手
================================

仓库的 ``.agents/skills/nepkappa-input/`` 保存输入助手技能，用于根据需求准备和检查 YAML。
技能以当前检出的 parser、示例和文档为准。只有能访问本地文件系统和已安装环境的助手，才能真正核对文件并执行校验。

使用示例
------------------------

在项目中向支持该技能的助手提出具体需求：

.. code-block:: text

   使用 $nepkappa-input，为体相 Si 准备 300 K 三声子 RTA 输入。
   使用 examples/structures/Si/POSCAR_bulk 和
   examples/potentials/Si/Si_Bulk_Fan.txt。
   保存为 si-rta.yaml 并验证，不启动计算。

复用 FC2 绘图：

.. code-block:: text

   使用 $nepkappa-input，检查 calculations/3C-SiC/results。
   只复用现有 FC2 和元数据，准备色散、DOS、体积热容与群速度的
   单图和组合图输入。先验证，不生成图像，也不重新计算力。

服务器环境排查：

.. code-block:: text

   使用 $nepkappa-input，检查服务器上的 VASP 与 POTCAR 路径，
   核对 input.yaml 的启动命令、MPI 环境和 Slurm 设置，不提交任务。

远程文件检查需要对目标主机的实际访问。提供“准备输入、绘图、还是执行计算”的明确范围，
可避免助手把已有结果复用扩展成重新计算。

路线与输出约定
------------------------

声子 / 输运采用六部分输入，QHA、SSCHA、四声子开关位于 ``force-constant``。
TD-BTE 从头计算使用六部分中的 ``tdbte``；已有核或 FC 使用精简 ``tdbte`` / ``output``。
旧的顶层 ``workflow``、``qha``、``scph`` 和 preset 格式不再作为公共输入。

助手应交付 YAML、启动目录、将要执行的命令、材料 / 方法假设，以及文件和配置的验证状态。
完整输入通常验证为：

.. code-block:: bash

   nepkappa info input.yaml --for run

单阶段换成 ``--for kappa``、``--for plot`` 等。
配置通过、文件存在、设备可用和物理模型有效是不同结论，应分别说明。
准备输入不会隐含授权启动力计算、绘图或提交调度作业。

技能位置与维护
------------------------

在当前仓库打开支持 Agent Skills 的助手时，使用仓库内的技能目录；
其他助手可按其规则安装或链接整个 ``nepkappa-input``，不要只复制一个文件而丢失引用资源。
例如从仓库根目录创建可选项目链接，且仅在目标不存在时执行：

.. code-block:: bash

   mkdir -p .claude/skills
   ln -s ../../.agents/skills/nepkappa-input .claude/skills/nepkappa-input

   mkdir -p .cursor/skills
   ln -s ../../.agents/skills/nepkappa-input .cursor/skills/nepkappa-input

复制后的技能仍需访问匹配版本的文档与示例。没有文件访问能力的助手只能提供文本，不能声称已在本地验证。
参数和模板变化时同步维护技能，并使用同一源码安装执行校验。
