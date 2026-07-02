NEP-kappa Documentation
=======================

Welcome to the documentation of **NEP-kappa**.

Source code repository:
`GitHub <https://github.com/lyushisyan/NEP-kappa>`_

**NEP-kappa v1.1** is an installable workflow package for lattice thermal
conductivity calculations based on **NEP**, **VASP**, **HiPhive**, and
**phono3py**.

It supports:

- staged commands: ``info``, ``relax``, ``fc``, ``kappa``, ``plot``, ``compare``, and ``run``
- optional NEP or VASP structure relaxation
- force-constant generation by finite displacement or HiPhive
- thermal conductivity calculations with phono3py
- plotting for dispersion, DOS, heat capacity, group velocity, relaxation time, and kappa
- DFT-vs-NEP comparison plots
- repository-provided YAML examples for bulk and film systems

.. toctree::
   :maxdepth: 1
   :caption: Table of contents
   :numbered:

   introduction
   requirements
   installation
   starting
   input_files
   tutorial
   reference
   troubleshooting
