Troubleshooting
=================

Configuration errors
------------------------------------------------

Run from the same launch directory and Python environment used for the job:

.. code-block:: bash

   nepkappa --version
   nepkappa validate input.yaml --for kappa
   nepkappa info input.yaml --for kappa

Replace ``kappa`` with the intended stage. Inspect ``run.log`` and the Slurm
job's stdout/stderr. ``nepkappa --debug <command> input.yaml`` executes the
command and prints a traceback on failure. ``validate`` and ``info`` do not
execute calculations. Compare/convergence inputs have separate parsers;
see :doc:`input_assistant`.

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
``calculator.vasp_command``, ``calculator.vasp_path``, and
``calculator.potcar_path`` before running on your machine.
``vasp_command`` takes precedence: changing only ``vasp_path`` cannot override
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

Supply both kernel NPZ and JSON metadata and explicitly enable experimental
mode. FC2/FC3 alone are not accepted as a kernel. Use a fresh result directory
when ``result-dir/tdbte`` already exists; preserve earlier audit results.
Failed numerical audits are not accepted trajectories. Passing audits still
does not validate absolute physical rates; see :doc:`tdbte`.

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
