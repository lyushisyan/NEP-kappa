NEP-kappa
=========

通过单一 YAML 配置文件，连接原子力计算、力常数提取、声子性质分析及热输运模拟。

- **多种计算后端**：支持 NEP、VASP、MACE 及其他 ASE 兼容计算器。
- **集成声子计算工具**：连接 Phonopy、phono3py 和 FourPhonon，实现声子性质、非简谐散射及晶格热导率计算。
- **非平衡声子动力学**：基于含时玻尔兹曼输运方程（TD-BTE），模拟非平衡声子分布的时间演化、能量弛豫与声子热化过程。

.. toctree::
   :maxdepth: 1
   :numbered: 2
   :caption: 用户指南

   introduction
   installation
   starting
   input_files
   harmonic_temperature
   thermal_transport
   tdbte

.. toctree::
   :maxdepth: 1
   :caption: 示例说明

   examples
   tutorial
   requirements
   input_assistant
   troubleshooting
   reference

.. toctree::
   :maxdepth: 1
   :caption: 后续开发

   development
   architecture

`English：本地英文预览 <../en/index.html>`_ ·
`GitHub 源码 <https://github.com/lyushisyan/NEP-kappa>`_ ·
`已发布的英文文档 <https://nep-kappa.readthedocs.io/en/latest/>`_
