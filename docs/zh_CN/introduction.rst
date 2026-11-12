
软件概述
====================

NEP-kappa 通过单一 YAML 配置文件，连接原子力计算、力常数提取、声子性质分析及热输运模拟。

- **多种计算后端**：支持 NEP、VASP、MACE 及其他 ASE 兼容计算器。
- **集成声子计算工具**：连接 Phonopy、phono3py 和 FourPhonon，实现声子性质、非简谐散射及晶格热导率计算。
- **非平衡声子动力学**：基于含时玻尔兹曼输运方程（TD-BTE），模拟非平衡声子分布的时间演化、能量弛豫与声子热化过程。

计算流程
--------------------

.. graphviz::
   :align: center
   :caption: NEP-kappa 计算工作流

   digraph NEPkappa {
       graph [rankdir=TB, bgcolor="transparent", pad="0.20",
              nodesep="0.35", ranksep="0.48", splines=polyline];

       node [shape=box, style="rounded,filled",
             color="#A7B8CF", fillcolor="#F3F6FA",
             fontname="sans-serif", fontsize=12,
             margin="0.18,0.12"];

       edge [color="#64748B", penwidth=1.3, arrowsize=0.65];

       structure [label="晶体结构弛豫"];
       forces [label="能量与原子力计算\nNEP / VASP / MACE / ASE"];
       ifcs [label="力常数生成\nFC2 / FC3 / FC4"];

       structure -> forces -> ifcs;

       thermal [label="声子性质与热输运\n调用 Phonopy / Phono3py / ShengBTE / FourPhonon", fillcolor="#DFEAF8"];
       td [label="非平衡声子动力学\n自研 TD-BTE", fillcolor="#DFEAF8"];

       { rank=same; thermal; td; }

       ifcs -> thermal;
       ifcs -> td;

       qha [label="QHA / SSCHA\n（可选）",
            style="rounded,dashed,filled", fillcolor="#FFFFFF"];
       phonon [label="色散 / DOS / 群速度 / 热容"];
       kappa [label="三声子 / 四声子散射\n晶格热导率（RTA / LBTE / Wigner）"];
       thermal_result [label="结果保存与绘图", fillcolor="#E9F4EF"];

       thermal -> qha [style=dashed];
       qha -> phonon [style=dashed];
       thermal -> phonon;
       phonon -> kappa -> thermal_result;

       excitation [label="初始声子分布与激发条件"];
       kernel [label="声子碰撞核与积分"];
       solver [label="求解含时 BT \n占据数演化 / 能量热化"];
       td_result [label="结果保存与绘图",
                  fillcolor="#E9F4EF"];

       td -> excitation -> kernel -> solver -> td_result;

       { rank=same; qha; excitation; }
       { rank=same; phonon; kernel; }
       { rank=same; kappa; solver; }
       { rank=same; thermal_result; td_result; }
   }

常用命令
--------------------

以下命令中的 ``input.yaml`` 为输入配置文件，可替换为实际文件路径。
通常先用 ``info`` 检查配置，再用 ``run`` 执行配置中启用的工作流步骤。

.. list-table::
   :header-rows: 1
   :widths: 14 36 46
   :width: 100%

   * - 类别
     - 命令
     - 用途
   * - 输入检查
     - ``nepkappa info input.yaml``
     - 校验并显示配置，不启动计算。
   * - 工作流
     - ``nepkappa run input.yaml``
     - 按顺序执行配置中启用的步骤。
   * - 结构优化
     - ``nepkappa relax input.yaml``
     - 弛豫晶体结构。
   * - 力常数
     - ``nepkappa fc2 input.yaml``
     - 生成二阶力常数。
   * - 力常数
     - ``nepkappa fc2fc3 input.yaml``
     - 生成二阶和三阶力常数。
   * - 力常数
     - ``nepkappa fc2fc3fc4 input.yaml``
     - 依次生成二阶、三阶和四阶力常数。
   * - 热输运
     - ``nepkappa kappa input.yaml``
     - 按配置选择的输运方案计算晶格热导率。
   * - 非平衡动力学
     - ``nepkappa tdbte input.yaml``
     - 求解含时 BTE，模拟声子分布演化。
   * - 后处理
     - ``nepkappa plot input.yaml``
     - 读取已有结果并绘图。
   * - 后处理
     - ``nepkappa report input.yaml``
     - 汇总已有结果，生成 Markdown 和 YAML 报告。
