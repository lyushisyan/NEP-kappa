依赖与运行要求
========================

Python 依赖
------------------------

标准安装由 ``pyproject.toml`` 声明依赖，包括 Python >= 3.10、NumPy、SciPy、ASE、h5py、
Matplotlib、Phonopy >= 4.1、phono3py >= 4.0.1、HiPhive、trainstation、seekpath、
Calorine >= 3.5、mace-torch、tqdm 和 PyYAML。安装步骤见 :doc:`installation`。

不同路线的附加要求
------------------------

* NEP 需要与结构元素和目标相匹配的模型。元素覆盖并不等于目标材料的精度保证。
* GPU 位移力需要 GPUMD 的 ``nep`` 程序和兼容 GPU；Python 安装不会安装 GPUMD 或驱动。
* VASP 需要目标主机上的可执行文件、MPI 环境与 POTCAR 文件 / 库，示例路径需要替换。
* MACE 后端随标准安装提供，计算时需指定匹配的模型文件或基础模型。
* 通用 ASE / 插件计算器需要另行安装提供该计算器的软件包。
* SSCHA 需要计算器返回总能和原子力。Calorine NEP 后端符合这一接口要求。
* Thirdorder、Fourthorder 与 FourPhonon 需独立安装；四声子 Wigner 使用 Wigner_Park。
* Slurm 需要 ``sbatch`` 及登录 / 计算节点共享的结果文件系统。
* TD-BTE 会自动使用 Numba；HiPhive 的依赖中已包括 Numba。不可用时使用 NumPy 参考实现。
* 文档构建使用 ``python -m pip install -e '.[docs]'`` 安装 Sphinx 和主题。

外部程序编译使用的数值库与 MPI 库必须在计算节点可用。Python 环境安装不会替你安装有授权要求的程序。
普通输入路径相对于命令启动目录解析；仓库公共示例默认从仓库根目录启动。
