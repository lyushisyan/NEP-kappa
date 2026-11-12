常见问题与排查
================================

先确认环境与阶段
------------------------

在与任务相同的环境和启动目录检查：

.. code-block:: bash

   nepkappa --version
   nepkappa info input.yaml --for kappa

把 ``kappa`` 换成实际阶段。检查结果目录的 ``run.log`` 和 Slurm 标准输出 / 错误文件。
``nepkappa --debug <命令> input.yaml`` 会真正执行该命令，并在失败时输出 traceback；
如果只想检查配置，应使用 ``info``。

日志安静或输出过少
------------------------

默认终端只显示简洁进度；``-v`` 显示详细设置与外部程序输出。
``info`` 始终显示解析后的配置。FourPhonon 各工作目录还有 ``fourphonon.log``。
可使用 ``tail -f`` 跟随日志；某些数值阶段本来就长时间不输出。
``output.progress: false`` 可关闭活动状态行。
部分 Phonopy 自动原胞 / Gamma 提示只写日志，其他警告仍保留完整日志诊断。

找不到文件或路径不符合预期
------------------------------

* 普通相对路径从启动目录解析，不是从 YAML 文件位置解析。
* 核对 ``output.result-dir``，并把力常数和匹配的结构元数据一起保存。
* 需要 ``POSCAR_relaxed`` 时，先弛豫；已有优化结构则关闭该输入的弛豫开关。
* 可编辑安装会直接跟随当前源码，即使版本号没有变化。检查实际导入的包位置。
* 新缓存格式会拒绝旧身份记录。不要手工复制或重写有效标记来强制复用旧力。

phono3py 参数报错
------------------------

项目要求 phono3py >= 4.0.1。遇到已移除选项或版本不匹配时，核对实际执行的程序与 Python 环境：

.. code-block:: bash

   command -v phono3py
   python -m pip install --upgrade phono3py
   python -m pip install --upgrade -e .

环境升级可能改变数值结果，应在保留旧结果的前提下重新验证。

VASP 路径与资源
------------------------

示例中的程序、MPI 和 POTCAR 路径需要替换。
``vasp-command`` 优先于 ``vasp-path``，只改后者不会覆盖已经配置的命令。
目标主机加载环境后用 ``command -v vasp_std`` 检查；PATH 未命中也可能只是未加载模块。
模块加载放在 batch / preamble 中，不放进可执行文件路径。不要为了找路径就在登录节点启动 VASP。

缺少图或绘图失败
------------------------

* FC2-only 需要匹配元数据，不需要势或 FC3。
* 多个 ``kappa-m*.hdf5`` 并存时明确选择网格，不要靠重命名伪造匹配。
* 没有线宽就没有散射率图，没有模式热导率就没有累积图。
* Wigner-only、分裂 q 点或额外展宽维度的文件不一定是完整标准输运结果。
* 散射率使用最近的已保存温度，核对程序给出的提示。
* 虚频要求稳定性与收敛排查；热容排除非正模式并不意味着该相已稳定。

已有数据复用例子见 :doc:`tutorial`。

TD-BTE 输入与审计失败
------------------------------

提供匹配的 ``force-constants`` 与 ``mesh``，或者提供已有 ``kernel``，二者不能同时指定。
旧 ``experimental`` 键不再使用。三维来源只允许三个近零 Gamma 平移，不支持不稳定模式或含 NAC 的来源。
不要裁剪频率或丢弃 Born 数据以绕过校验。

已有 ``tdbte`` 目录或未完成核需要保留诊断，选择新的输出目录再运行。
数值审计失败的轨迹不能作为已接受结果；审计通过也不等于绝对散射速率已验证。
空间路线还需检查 CFL、占据数正性、连续性残差和数值扩散，见 :doc:`tdbte_spatial`。

Slurm 中的准备、提交与完成
--------------------------------

没有 ``sbatch`` 时，可用 ``parallel.kappa.submit: false`` 生成脚本。
结果目录必须是共享文件系统。RTA 的记录在 ``rta-slurm/submission.yaml``，
LBTE 在 ``lbte-slurm/submission.yaml``，对应日志在各自 ``logs/`` 下。
脚本自动设置 Rayon / OpenMP 线程，BLAS 线程需按库单独设置。

.. code-block:: bash

   squeue -u "$USER"
   sacct -j JOB_ID --format=JobID,State,Elapsed,ExitCode

``prepared`` 只表示准备成功，``submitted`` 只表示提交成功，科学结果仍未完成。
必须在 worker 成功收集、日志和退出码正常后再进行结果解释。
任务从队列消失不足以判断成功。没有 ``squeue`` 时，保存的 YAML 仍可查看。

联系与问题报告
------------------------

联系邮箱：``sxliu98@gmail.com``、``yinfei0426@outlook.com``。
报告问题时附上最小输入、启动目录、版本、错误日志和所用后端；避免发送授权势库或不应公开的数据。
