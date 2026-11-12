NEP-kappa
=========

NEP-kappa 通过统一的 YAML 输入文件，连接原子力计算、力常数生成、声子性质分析、晶格热输运及非平衡声子动力学模拟。

- **计算后端**：NEP、VASP、MACE 及其他 ASE 兼容计算器。
- **声子与输运**：集成 Phonopy、phono3py 和 FourPhonon。
- **非平衡动力学**：支持空间均匀及一维周期空间的含时 BTE。

建议初次使用的读者先完成第三、四章的两个示例，再根据需要阅读参数参考和专题章节。

.. toctree::
   :maxdepth: 1
   :numbered: 2
   :caption: 用户指南

   introduction
   installation
   starting_phonons
   starting
   workflow
   input_files
   temperature
   thermal_transport
   tdbte
   tdbte_spatial

.. toctree::
   :maxdepth: 1
   :caption: 项目与支持

   reference
   contact
   版本更新 <https://github.com/lyushisyan/NEP-kappa/releases>

`GitHub 源码 <https://github.com/lyushisyan/NEP-kappa>`_ ·
`英文在线文档 <https://nep-kappa.readthedocs.io/en/latest/>`_
