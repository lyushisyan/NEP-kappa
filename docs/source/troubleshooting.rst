Troubleshooting
=================

``phono3py`` option errors
----------------------------

NEP-kappa targets ``phono3py>=4.0.1``. If ``phono3py`` reports removed options
such as ``--fc2``, ``--fc3``, or ``--wigner``, upgrade the environment and
reinstall NEP-kappa:

.. code-block:: bash

   python -m pip install --upgrade phono3py
   python -m pip install --upgrade -e .

VASP path errors
------------------

The repository-provided VASP examples contain placeholder paths. Edit
``calculator.vasp_command``, ``calculator.vasp_path``, and
``calculator.potcar_path`` before running on your machine.

Slow runs
-----------

For long VASP or large ``phono3py`` calculations, use a scheduler such as Slurm
when available. For workstation testing, run inside ``tmux`` so the calculation
continues after disconnecting.

Slurm LBTE submission
-----------------------

If ``sbatch`` is unavailable, set ``kappa.parallel.submit: false`` to generate
the scripts without submitting them. The configured ``output.result_dir`` must
be on a filesystem shared by the login and compute nodes. Job IDs and generated
script paths are recorded in ``output.result_dir/lbte-slurm/submission.yaml``;
Slurm stdout and stderr files are written under ``lbte-slurm/logs``.

To inspect every force-array, LBTE, and FourPhonon submission associated with
one input file, run:

.. code-block:: bash

   nepkappa status input.yaml

The command still reports the stored submission state when ``squeue`` is not
available, for example after copying a result directory to another computer.

Questions
-----------

For questions, contact ``sxliu98@gmail.com`` or ``yinfei0426@outlook.com``.
