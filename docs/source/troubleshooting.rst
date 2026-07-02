Troubleshooting
===============

``phono3py`` option errors
--------------------------

NEP-kappa targets ``phono3py>=4.0.1``. If ``phono3py`` reports removed options
such as ``--fc2``, ``--fc3``, or ``--wigner``, upgrade the environment and
reinstall NEP-kappa:

.. code-block:: bash

   python -m pip install --upgrade phono3py
   python -m pip install --upgrade -e .

VASP path errors
----------------

The repository-provided VASP examples contain paths for a specific server. Edit
``calculator.vasp_command``, ``calculator.vasp_path``, and
``calculator.potcar_path`` before running on a different machine.

Slow runs
---------

For long VASP or large ``phono3py`` calculations, use a scheduler such as Slurm
when available. For workstation testing, run inside ``tmux`` so the calculation
continues after disconnecting.

Questions
---------

For questions, contact ``sxliu98@gmail.com`` or ``yinfei0426@outlook.com``.
