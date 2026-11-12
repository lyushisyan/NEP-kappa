科学回归基线
================================

NEP-kappa 2.0 的本地科学基线以 ``benchmarks/scientific-baselines.yaml`` 保存紧凑的可机读结果。
它冻结部分已记录工作流的观测量，避免把大规模位移与输运目录加入仓库。

基线范围
------------------------

* 3C-SiC 的 DFT、专用 NEP 与 MACE 在 300 K 的 RTA 热导率随 q 网格变化。
* 3C-SiC NEP 的代表温度 QHA 性质。
* 已完成的 3C-SiC 300 K Phonopy SSCHA 冒烟结果。
* BAs DFT 与 NEP 的 FC2、FC3、FC4 比较。
* BAs DFT 与 NEP 的三加四声子热导率。

历史记录中未完成的 BAs MACE 项只是归档状态，不是实时调度器状态，也不能作为参考结果。

如何使用
------------------------

``benchmarks/`` 和相应本地测试由 Git 忽略，需要在本机实际存在。
BAs 来源摘要在 ``benchmarks/data/BAs/sources/``，可检查哈希而不读取私有计算目录。

.. code-block:: bash

   python -m pytest tests/test_scientific_baselines.py -q

快照检查不重新进行 DFT 或生产计算，只验证归档记录。
严格快照容差用于发现数据或解析变化，较宽的科学容差用于评估有意的算法 / 数值库升级。
通过归档检查不能替代重新计算的生产收敛，也不能证明新硬件或环境与旧结果完全等价。
