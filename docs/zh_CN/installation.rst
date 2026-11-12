
安装与环境配置
========================

NEP-kappa 支持 Python 3.10 及以上版本，推荐使用 Python 3.11。
Python 依赖（含 MACE）由 ``pyproject.toml`` 统一管理，外部求解程序需按计算任务配置。

快速安装
------------------------

推荐使用 Conda 创建独立的 Python 环境：

.. code-block:: bash

   conda create -n nepkappa python=3.11
   conda activate nepkappa

下载源码并安装：

.. code-block:: bash

   git clone https://github.com/lyushisyan/NEP-kappa.git
   cd NEP-kappa
   python -m pip install -e .

以上命令以可编辑模式安装 NEP-kappa，并一并安装 MACE 后端，适合需要更新或开发源码的用户。

检查安装结果：

.. code-block:: bash

   python -m pip check
   nepkappa --version
   nepkappa info examples/nep-rta-wigner-3ph.yaml

如果命令正常执行，说明软件已安装，并能解析示例输入文件。
``nepkappa info`` 仅检查输入配置，不会启动实际计算。

.. note::

   也可以使用 Python 自带的 ``venv`` 创建虚拟环境：

   .. code-block:: bash

      python -m venv .venv
      source .venv/bin/activate

计算后端与附加依赖
------------------------

上述安装步骤包含 NEP、MACE、Phonopy、phono3py 和 HiPhive 等依赖。
以下外部程序需要根据研究任务配置。

.. list-table::
   :header-rows: 1
   :widths: 25 75
   :width: 100%

   * - 计算后端
     - 附加要求
   * - **VASP**
     - 安装具有有效授权的 VASP，并配置 MPI、POTCAR 和相关环境。
   * - **Thirdorder / Fourthorder**
     - 配置 ``thirdorder_vasp.py`` 和 ``Fourthorder_vasp.py``。
   * - **FourPhonon**
     - 单独安装 FourPhonon；四声子 Wigner 计算还需要 Wigner_Park。

.. note::

   - 外部程序的安装与编译需要适配当前操作系统、编译器及 MPI 环境。
   - 并非所有计算路线都需要上述外部程序，按需安装即可。

更新软件
------------------------

进入已有的 NEP-kappa 目录，获取最新代码并更新安装：

.. code-block:: bash

   conda activate nepkappa
   cd /path/to/NEP-kappa
   git status
   git pull --ff-only
   python -m pip install -e .

检查更新后的环境：

.. code-block:: bash

   python -m pip check
   nepkappa --version

更新前建议检查本地修改，避免与远程代码发生冲突。

完成安装后，参见 :doc:`starting` 运行第一个示例。
