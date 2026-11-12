示例目录
========================

仓库包含九个声子 / 热输运输入与两个 TD-BTE 输入，均使用完整六部分布局。
下面的公共模板默认从仓库根目录运行，数值设置用于展示流程，不代表材料收敛结果。

.. list-table:: 公共输入
   :header-rows: 1
   :widths: 52 48

   * - 文件
     - 用途
   * - :download:`vasp-rta-3ph.yaml <../../examples/vasp-rta-3ph.yaml>`
     - VASP 原子力与三声子 RTA。
   * - :download:`nep-rta-wigner-3ph.yaml <../../examples/nep-rta-wigner-3ph.yaml>`
     - NEP、三声子 RTA 与 SMM19 Wigner。
   * - :download:`nep-lbte-wigner-3ph.yaml <../../examples/nep-lbte-wigner-3ph.yaml>`
     - NEP、三声子 LBTE 与 SMM19 Wigner。
   * - :download:`nep-rta-3ph-4ph.yaml <../../examples/nep-rta-3ph-4ph.yaml>`
     - NEP、FourPhonon 三加四声子 RTA。
   * - :download:`nep-rta-wigner-3ph-4ph.yaml <../../examples/nep-rta-wigner-3ph-4ph.yaml>`
     - Wigner_Park 粒子项与相干项。
   * - :download:`nep-lbte-3ph-rta-4ph.yaml <../../examples/nep-lbte-3ph-rta-4ph.yaml>`
     - 三声子 LBTE 与四声子 RTA 混合解。
   * - :download:`nep-lbte-3ph-lbte-4ph.yaml <../../examples/nep-lbte-3ph-lbte-4ph.yaml>`
     - 三、四声子均使用迭代 LBTE。
   * - :download:`nep-qha-rta-3ph.yaml <../../examples/nep-qha-rta-3ph.yaml>`
     - 在各 QHA 平衡体积重新生成 FC2/FC3 并计算 RTA。
   * - :download:`nep-qha-sscha-rta-3ph.yaml <../../examples/nep-qha-sscha-rta-3ph.yaml>`
     - QHA 体积、SSCHA FC2 与三声子 RTA。
   * - :download:`tdbte.yaml <../../examples/tdbte.yaml>`
     - 先生成 FC2/FC3，再演化空间均匀占据数。
   * - :download:`tdbte-spatial.yaml <../../examples/tdbte-spatial.yaml>`
     - 三维 q 网格、周期一维实空间的非线性动力学。

检查输入而不运行：

.. code-block:: bash

   nepkappa info examples/nep-rta-wigner-3ph.yaml --for run

执行完整路线用 ``run``，绘图和报告另行运行。直接 ``tdbte`` 只复用已有力常数，
完整 ``run`` 会先生成它们。TD-BTE 不需要已有热导率结果。

模型与路径
------------------------

公共晶体结构在 ``examples/structures/``，模型及来源说明在 ``examples/potentials/``。
3C-SiC 的专用模型是 ``examples/potentials/3C-SiC/nep_3C-SiC.txt``；
NEP89 是不同的通用模型，不能只因支持 C 和 Si 就将两者互换。

VASP 与 FourPhonon 程序路径是占位设置，必须替换为目标机器的安装路径。
计算案例和结果放在被 Git 忽略的 ``calculations/``，不要把本地集群账号和资源要求写入公共模板。

输入变体
------------------------

固定体积 SSCHA、QHA 和 FC2 结果复用的可下载示例见 :doc:`harmonic_temperature` 与 :doc:`tutorial`。
普通三声子输入使用 ``kappa.method``；FourPhonon 使用两个独立 method，支持 RTA/RTA、LBTE/RTA、LBTE/LBTE。
Wigner_Park 要求 RTA/RTA。QHA 和 SSCHA 的 ``run-transport`` 与温度关系见 :doc:`input_files`。

极性材料与收敛
------------------------

SiC 是极性材料。需要定量色散或输运结论时，应评估长程非解析修正，并提供经过验证的 Born 电荷和介电数据。
示例默认设置不自动完成这项验证。力常数超胞、截断、q 网格、采样和所用势函数都需要各自的收敛或精度检查。
