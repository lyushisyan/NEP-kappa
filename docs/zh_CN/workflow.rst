
计算工作流与命令
========================

NEP-kappa 使用 YAML 文件统一配置晶体结构、原子力计算、
力常数生成、晶格热输运及非平衡声子动力学。

用户可以通过 ``nepkappa run`` 执行完整工作流，
也可以使用独立命令只运行指定的计算阶段。

本页介绍输入文件的组织方式、各模块的常用参数及其适用范围。
首次使用的读者建议先阅读 :doc:`starting`。


完整工作流
------------------------

完整的 YAML 输入文件包含六个主要部分：

.. list-table:: YAML 配置模块
   :header-rows: 1
   :widths: 30 70

   * - 模块
     - 功能
   * - ``structure``
     - 晶体结构与结构弛豫。
   * - ``calculator``
     - 原子间相互作用模型与能量、原子力计算。
   * - ``force-constant``
     - 二阶、三阶和四阶力常数的生成。
   * - ``kappa`` / ``tdbte``
     - 晶格热输运或非平衡声子动力学。
   * - ``plot``
     - 绘图内容、格式及输出设置。
   * - ``output``
     - 计算结果的保存路径。

``kappa`` 和 ``tdbte`` 对应两条独立的计算路线，
不能在同一个完整输入文件中同时使用。

未设置具体参数的模块可以写为 ``{}``。

保存输入文件为 ``input.yaml`` 后，可以使用以下命令：

.. code-block:: bash

   nepkappa info input.yaml
   nepkappa run input.yaml

``info`` 检查输入配置，``run`` 按照 YAML 设置自动组织并执行各计算阶段。

例如，三声子热输运工作流依次完成结构弛豫、FC2/FC3 生成和热导率计算，并保存各阶段的结果以供后续使用。


单阶段计算
------------------------

如果只需要执行部分计算，可以使用独立命令，
无需重新运行整个工作流。

例如，已生成 FC2/FC3 后，可以使用 ``nepkappa kappa``
直接计算不同温度或 q 网格下的热导率，而不重新计算原子力。

此时，输入文件只需要包含 ``kappa`` 和 ``output`` 等相关模块。
使用 ``kappa.force-constants`` 指定已有力常数目录，
使用 ``output.result-dir`` 指定新的结果目录。

.. code-block:: bash

   nepkappa info kappa.yaml --for kappa
   nepkappa kappa kappa.yaml

第一条命令检查单阶段输入，第二条执行热导率计算。

其他单阶段命令如下：

.. list-table:: 单阶段命令
   :header-rows: 1
   :widths: 28 72

   * - 命令
     - 功能与所需配置
   * - ``relax``
     - 结构弛豫，需要 ``structure`` 和 ``calculator``。
   * - ``fc2``
     - 生成 FC2，需要 ``structure``、``calculator`` 和 ``force-constant``。
   * - ``fc2fc3``
     - 生成 FC2/FC3，使用相同的配置模块。
   * - ``fc2fc3fc4``
     - 生成 FC2/FC3/FC4，使用相同的配置模块。
   * - ``kappa``
     - 使用已有力常数计算热导率，需要 ``kappa``。
   * - ``tdbte``
     - 使用已有力常数或碰撞核求解含时 BTE，需要 ``tdbte``。
   * - ``plot``
     - 读取已有结果并绘图，使用 ``plot`` 和 ``output``。
   * - ``report``
     - 生成计算报告，指定 ``output.result-dir`` 或直接提供结果目录。

独立命令既可以使用精简输入，也可以读取完整 YAML 中的相关配置。

.. note::

   - 相对路径以运行命令时的工作目录为准，而不是 YAML 文件所在目录。
   - 建议为不同计算使用独立的 ``output.result-dir``，避免结果混淆。
   - 精简输入应通过 ``nepkappa info input.yaml --for <命令>`` 检查配置。
   - 如果输入启用了结构弛豫，单独运行力常数命令前应先执行 ``relax``。
