开发与文档维护
================================

版本控制范围
------------------------

仓库维护源码、文档、公共示例、模型来源说明、输入助手和 ``tests/public/``。
``tests/gpu/`` 中公开两个可选 GPU 验证 / 基准脚本。
其余本地测试、``benchmarks/``、整个 ``calculations/``、``tmp/``、``tools/``、``constraints/``
及生成输出目录由 Git 忽略。不要把私有计算、站点脚本或历史任务状态混入公共示例。

依赖本地 ``tools/`` 的可选审计测试在脚本缺失时跳过；求解器自包含检查仍然运行。
软件依赖由 ``pyproject.toml`` 声明，CI 不依赖本地 ``constraints/`` 快照。

运行检查
------------------------

在仓库根目录执行：

包含 ``.. graphviz::`` 的页面还需要 Graphviz 的 ``dot`` 程序。
macOS 使用 ``brew install graphviz``，Ubuntu 使用 ``sudo apt-get install graphviz fonts-noto-cjk``。
可用 ``dot -V`` 检查安装；Sphinx 的图扩展已包含在共享配置中，HTML 输出采用 SVG。

.. code-block:: bash

   python -m pip install -e '.[dev,docs,tdbte]'
   OMP_NUM_THREADS=1 RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
     python -m pytest tests/public -q
   python -m sphinx -b html -W --keep-going docs/source docs/build/en
   python -m sphinx -b html -W --keep-going docs/zh_CN docs/build/zh_CN

明确验证两条碰撞路径：

.. code-block:: bash

   NEPKAPPA_TEST_COLLISION_BACKEND=numpy python -m pytest tests/public -q
   NEPKAPPA_TEST_COLLISION_BACKEND=numba python -m pytest tests/public -q

这是测试专用选择，不改变生产 YAML。HiPhive 依赖的 Numba 仍可用；编译路径要求 Numba 实际可用。
NumPy 检查保留参考一致性与小规模收敛；完整收敛基准只在相应后端且本地工具可用时运行。
CI 使用 Python 3.10 / 3.12 的两条路径。GPU 与生产规模材料收敛需要独立环境和资源。

本地补充测试可用 ``python -m pytest tests --ignore=tests/gpu -q``。
只有具备所需私有基线和数据时才应比较相同范围的本地测试总数。
根目录不限定路径的 pytest 还可能发现 ``calculations/`` 下的研究测试，范围并不相同。

中英文文档独立维护
--------------------------

英文源在 ``docs/source/``，中文源在 ``docs/zh_CN/``，同名章节分别撰写。
``docs/conf_common.py`` 共用主题和 Sphinx 扩展；两份 ``conf.py`` 分别指定语言和页面标题。
两种语言共用 ``docs/source/inputs/`` 的 YAML 和 ``docs/source/_static/`` 的图像 / 数据，不复制数值输入。
中文使用 ``literalinclude`` 和 ``download`` 引用这些文件。

维护时保持命令、参数键、单位、支持范围和科学结论一致，说明文字可以采用不同表述。
先更新相关页面，再对两种语言执行严格构建。不要把“仅把导航语言改成中文”当作完成正文翻译。
本地预览目录的 ``zh_CN/index.html`` 可直接打开，也可以启动只绑定本机的 HTTP 服务：

.. code-block:: bash

   python -m http.server 8765 --bind 127.0.0.1 --directory docs/build

然后访问 ``http://127.0.0.1:8765/zh_CN/``。
在线发布的 Read the Docs 配置仍是英文入口；本地中文预览不自动配置线上翻译项目或发布内容。

关键回归保护
------------------------

公开检查覆盖解析、阶段选择、延期状态、力缓存身份 / 损坏 / 中断、worker / 收集一致性、
小型能量守恒、能量壳矩、输运单位归一化和空间 TD-BTE。
GPU 校验脚本为 ``tests/gpu/validate_nep_batch.py`` 和 ``benchmark_nep_batch.py``，需在真实设备上按说明运行。
实际 Slurm 集群性能与生产材料收敛不能由本地模拟脚本替代。

RTA 集群检查会在本地执行生成的 prepare / worker / collect 脚本并模拟 Slurm 变量，
检查结果与串行值、旧整网格文件隔离和缺失 worker 数据失败处理；它不会提交真实任务。

公共接口与示例
------------------------

新增阶段通过应用层接入，后端细节留在计算器 / adapter 中。
参数变化时同时更新 parser、校验、公共 YAML、两种语言文档以及 ``.agents/skills/nepkappa-input/``。
公共模板包含九个声子 / 输运和两个 TD-BTE 输入；资源占位值需要由使用者在目标站点替换。

``info`` 验证不能被描述为完成计算或验证物理精度。
缓存变更需升级身份并明确失效旧记录，不能伪造历史摘要来复用旧力。
科研基线说明见 :doc:`release_baseline`，碰撞核独立诊断见 :doc:`tdbte_validation`。

.. toctree::
   :maxdepth: 1

   release_baseline
