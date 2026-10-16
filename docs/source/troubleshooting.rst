Troubleshooting
=================

Configuration errors
------------------------------------------------

Run from the same launch directory and Python environment used for the job:

.. code-block:: bash

   nepkappa --version
   nepkappa info input.yaml --for kappa

Replace ``kappa`` with the intended stage. Inspect ``run.log`` and the Slurm
job's stdout/stderr. ``nepkappa --debug <command> input.yaml`` executes the
command and prints a traceback on failure. ``info`` does not execute calculations.

Screen output and logs
----------------------

Calculation commands use compact output by default. Add ``-v`` or ``--verbose``
to show detailed settings and external-program output. ``info`` always displays
the full parsed configuration. Detailed diagnostics are saved in ``run.log``;
FourPhonon also writes ``fourphonon.log`` in each transport working directory.

Interactive transport runs show an updating elapsed-time line. Redirected output
has no dynamic updates; ``output.progress: false`` disables the status line.
Logs are flushed as external-program lines arrive, so ``tail -f`` can follow
them. A quiet interval can reflect a computational stage that produces no output.

The Phonopy warnings ``PrimitiveMatrixAutoDefaultWarning`` and
``GammaAcousticWarning`` are logged without appearing on screen, including in
verbose mode. Other warnings retain their full diagnostics in the log while
the compact screen output omits Python source locations.

Missing files or unexpected defaults
--------------------------------------

Ordinary relative paths resolve from the launch directory, not the YAML's
directory. Check ``output.result-dir`` and keep force constants together with
their original structure metadata. If a stage expects ``POSCAR_relaxed``,
complete relaxation first or explicitly disable relaxation for an already
relaxed input structure.

Use the source and executable from the same installation. An editable install
follows local changes even when the displayed version number is unchanged. Supported
keys are listed in :doc:`input_files`.

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
``calculator.vasp-command``, ``calculator.vasp-path``, and
``calculator.potcar-path`` before running on your machine.
``vasp-command`` takes precedence: changing only ``vasp-path`` cannot override
it. Check ``command -v vasp_std`` on the execution host after loading its
environment. No PATH match may mean a missing module, not an absent installation.
Put module setup in a batch script/preamble, not in the command string. Never
start VASP on a login node merely to discover its location.

Plotting errors or missing panels
----------------------------------

* FC2-only plots require matching ``phono3py_disp.yaml``, ``phonopy.yaml``, or
  ``phonopy_disp.yaml``. FC3 and a potential are not needed.
* A configured mesh must match the conductivity filename. If several files
  exist, select the intended mesh rather than renaming files.
* No linewidths means no lifetime/scattering panels; no ``mode_kappa`` means
  no cumulative-conductivity panel. Comparisons show only shared capabilities.
* Wigner-only HDF5, split grid-point output, and extra broadening axes are not
  supported by the standard transport reader. Do not treat partial results
  as a completed standard kappa file. Custom manuscript plots are separate.
* Transport scatter plots use the nearest stored temperature and warn when it
  differs from the request. Check the reported temperature.
* Imaginary frequencies require stability and convergence checks.
  Harmonic heat capacity excludes nonpositive modes
  and is not a stable-phase prediction when imaginary modes remain.

See :doc:`tutorial` for a minimal harmonic-only plotting input. ``report``
summarizes artifacts; it does not calculate missing quantities.

TD-BTE input and audit errors
-----------------------------------

Supply ``force-constants`` (a directory with matching phono3py YAML, FC2 and FC3)
and ``mesh``, or supply ``kernel`` (NPZ plus sidecar, or a completed packaged
chunk manifest), not both. Remove the old ``experimental`` key. The source must
have only three near-zero Gamma translations and no unstable modes; NAC-bearing
inputs are not yet supported by the builder. Do not bypass these checks by
clipping frequencies or discarding Born charges. Use a fresh result directory
when ``result-dir/tdbte`` exists or a kernel build is incomplete; preserve the
earlier audit and ``build-failure.json`` for diagnosis.
Failed numerical audits are not accepted trajectories. Passing audits still
does not validate absolute physical rates; see :doc:`tdbte`.

Slow runs
-----------

For long VASP or large ``phono3py`` calculations, use a scheduler such as Slurm
when available. For workstation testing, run inside ``tmux`` so the calculation
continues after disconnecting.

Slurm LBTE submission
-----------------------

If ``sbatch`` is unavailable, set ``parallel.kappa.submit: false`` to generate
the scripts without submitting them. The configured ``output.result-dir`` must
be on a filesystem shared by the login and compute nodes. Job IDs and generated
script paths are recorded in ``<result-dir>/lbte-slurm/submission.yaml``;
Slurm stdout and stderr files are written under ``lbte-slurm/logs``.

To inspect active Slurm jobs, run:

.. code-block:: bash

   squeue -u "$USER"

Saved submission details remain in each result directory's submission YAML
files, including when ``squeue`` is unavailable.

Questions
-----------

For questions, contact ``sxliu98@gmail.com`` or ``yinfei0426@outlook.com``.
