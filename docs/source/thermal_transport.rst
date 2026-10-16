Phonon scattering and lattice thermal conductivity
======================================================

This chapter compares three-phonon LBTE, three-plus-four-phonon transport, and
Wigner particle/coherence contributions. It continues from the force constants
prepared in :doc:`starting`.

Workflow and preparation
----------------------------

.. code-block:: bash

   mkdir -p ~/work/3C-SiC/transport
   cd ~/work/3C-SiC/transport
   ls -lh ../results/fc2.hdf5 ../results/fc3.hdf5 ../results/phono3py_disp.yaml
   cp ../results/POSCAR_relaxed POSCAR
   cp ../nep_3C-SiC.txt .

The LBTE example reuses the existing FC2/FC3. The four-phonon example generates
a matching FC2/FC3/FC4 set, which is then reused by the Wigner example.
All force constants, cell coordinates and lattice-translation labels must
follow the same convention; see :doc:`input_files` for reuse checks.

Three-phonon scattering and LBTE
------------------------------------

Input and execution
~~~~~~~~~~~~~~~~~~~~~~~

RTA treats the departure of each mode from equilibrium independently. LBTE
retains linear coupling between modes through the collision matrix.
Save :download:`lbte.yaml <inputs/lbte.yaml>`:

.. literalinclude:: inputs/lbte.yaml
   :language: yaml

.. code-block:: bash

   nepkappa info lbte.yaml --for kappa
   nepkappa kappa lbte.yaml
   tail -n 20 results/lbte/run.log
   ls -lh results/lbte/kappa-m151515.hdf5

``kappa`` reuses the specified force constants and copies the matching inputs
into the result directory for subsequent plotting. LBTE collision matrices can
require substantially more memory than RTA; the example mesh is a demonstration
setting, not a convergence result.

Plot and compare
~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   nepkappa info lbte.yaml --for plot
   nepkappa plot lbte.yaml
   ls results/lbte/plots

The saved ``kappa`` and ``kappa_RTA`` tensors provide LBTE and RTA conductivities
from the same calculation. The conductivity panel overlays them as solid and
dashed curves. ``plot.kappa: x`` selects the xx component; ``all`` selects the
mean of the three diagonal components.

The six-panel combination contains four harmonic properties, the scattering
rate at 300 K, and conductivity versus temperature. ``tau: total`` displays
the saved ph-ph and isotope channels and their sum; boundary scattering is not
included in that plotted sum. The display threshold is 1e-8 ps^-1. N/U plots
require separately saved N and U linewidths, which the default LBTE command
does not request. Cumulative conductivity requires mode-resolved conductivity.

Four-phonon scattering and conductivity
-------------------------------------------

Input and execution
~~~~~~~~~~~~~~~~~~~~~~~

Install Thirdorder, Fourthorder and the FourPhonon CPU executable as described
in :doc:`installation`. Check their discovery:

.. code-block:: bash

   command -v thirdorder_vasp.py
   command -v Fourthorder_vasp.py
   command -v ShengBTE_cpu

Save :download:`fourphonon.yaml <inputs/fourphonon.yaml>`:

.. literalinclude:: inputs/fourphonon.yaml
   :language: yaml

Thirdorder generates FC3 and Fourthorder generates FC4 from NEP displacement
forces. Negative cutoffs select neighbor shells. ``format: both`` preserves
the ShengBTE text force constants and exports the HDF5 forms used for analysis.
The example requests four MPI processes with one OpenMP thread each; adjust
these to the resources available on your machine.

``compare-3ph: true`` runs a three-phonon reference and a three-plus-four-phonon
calculation sequentially with the same FourPhonon backend, structure, FC2/FC3,
mesh, temperatures and isotope settings. Here the three-phonon part uses LBTE
and the four-phonon part uses RTA.

.. _transport-tmux:

Long calculations in tmux
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Check the input, then create a session on the calculation server:

.. code-block:: bash

   nepkappa info fourphonon.yaml --for run
   tmux new -s fourphonon

Inside the session, activate the environment and start the calculation:

.. code-block:: bash

   conda activate nepkappa
   cd ~/work/3C-SiC/transport
   nepkappa run fourphonon.yaml

Press ``Ctrl-b``, release both keys, then press ``d`` to detach. The calculation
continues on the server when the SSH connection closes. Reconnect with:

.. code-block:: bash

   tmux ls
   tmux attach -t fourphonon

Use ``Ctrl-b`` then ``[`` to scroll earlier output, and ``q`` to leave scroll
mode. tmux keeps the shell session alive; it does not replace a cluster scheduler
or extend a scheduler allocation.

Thirdorder, Fourthorder and FourPhonon use compact screen output by default.
Add ``-v`` to the calculation command for detailed screen output. Details are
retained in ``run.log``; each FourPhonon working directory also contains
``fourphonon.log``.

Results and comparison plots
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After completion, inspect the log before running plotting:

.. code-block:: bash

   tail -n 20 results/fourphonon/run.log
   ls results/fourphonon/fourphonon/threephonon
   ls results/fourphonon/fourphonon
   nepkappa info fourphonon.yaml --for plot
   nepkappa plot fourphonon.yaml
   ls results/fourphonon/plots

The reference directory contains ``kappa3-rta.dat`` and
``kappa3-iterative.dat``. The combined directory contains ``kappa4-rta.dat``
and ``kappa4-iterative.dat``. With this input, the latter stores the mixed
3ph LBTE + 4ph RTA solution. If both ``method-3ph`` and ``method-4ph`` are
``lbte``, it instead stores the iterative three-plus-four-phonon solution.

The conductivity panel overlays four curves: 3ph RTA, 3ph LBTE, 3ph+4ph RTA,
and 3ph LBTE + 4ph RTA. Three- and four-phonon scattering rates share one panel,
``scattering_rate_3ph_4ph.png``. Together with dispersion, DOS, heat capacity
and group velocity, this gives six panels and seven PNG files with
``layout: both``. The plotted 3ph and 4ph rate series do not include the isotope
channel, although isotope scattering is enabled in the conductivity calculation.

Wigner transport with three- and four-phonon scattering
-----------------------------------------------------------

Prepare existing force constants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Wigner transport resolves conductivity into particle and coherence
contributions, with their sum giving total conductivity. This example uses
RTA for both scattering orders and the FourPhonon Wigner_Park executable.
It runs only the combined 3ph+4ph case.

Create a separate working directory on the target cluster. Transfer the
completed force-constant set into ``inputs/``, keeping these files together:

.. code-block:: text

   inputs/
     POSCAR_relaxed
     FORCE_CONSTANTS_2ND
     FORCE_CONSTANTS_3RD
     FORCE_CONSTANTS_4TH
     fc2.hdf5
     fc3.hdf5
     phono3py_disp.yaml

Save :download:`wigner-4ph.yaml <inputs/wigner-4ph.yaml>` and replace the
executable placeholder with the installed Wigner_Park binary:

.. literalinclude:: inputs/wigner-4ph.yaml
   :language: yaml

.. code-block:: bash

   nepkappa info wigner-4ph.yaml --for kappa

No calculator is needed because ``kappa`` reuses the prepared force constants.
Automatically generated CONTROL files disable non-analytic corrections.
Calculations requiring them need a consistent custom CONTROL with validated
Born charges and dielectric data.

.. _transport-slurm:

Submit through Slurm
~~~~~~~~~~~~~~~~~~~~~~~~

Save the following as ``submit-wigner-4ph.sh``. Replace ``YOUR_PARTITION`` and
the environment path with your site's settings. Its four tasks match
``parallel.kappa.mpi-processes: 4`` in the downloadable YAML.

.. code-block:: bash

   #!/bin/bash -l
   #SBATCH --job-name=SiC-Wigner-4ph
   #SBATCH --partition=YOUR_PARTITION
   #SBATCH --nodes=1
   #SBATCH --ntasks=4
   #SBATCH --cpus-per-task=1
   #SBATCH --time=48:00:00
   #SBATCH --output=logs/wigner-4ph-%j.out
   #SBATCH --error=logs/wigner-4ph-%j.err

   set -eo pipefail
   cd "$SLURM_SUBMIT_DIR"
   source /path/to/miniconda3/etc/profile.d/conda.sh
   conda activate nepkappa
   export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
   nepkappa kappa wigner-4ph.yaml

Load any required compiler/MPI modules before the command, following your site's
setup. NEP-kappa launches the MPI backend from the YAML settings; this script
runs the NEP-kappa command once within the allocation.

.. code-block:: bash

   mkdir -p logs
   bash -n submit-wigner-4ph.sh
   sbatch submit-wigner-4ph.sh
   squeue -u "$USER"

``bash -n`` checks shell syntax. ``sbatch`` submits the calculation and prints
its job ID. After completion, substitute that ID for ``JOB_ID`` below:

.. code-block:: bash

   sacct -j JOB_ID --format=JobID,State,Elapsed,ExitCode
   tail -n 20 results/wigner-4ph/run.log
   ls results/wigner-4ph/fourphonon

Check ``COMPLETED``, exit code ``0:0`` and ``Completed: kappa`` in the log.
Disappearance from ``squeue`` alone does not establish success.

Plot and interpret
~~~~~~~~~~~~~~~~~~~~~~

The FourPhonon stage links its text force constants into the backend working
directory. Standard plotting also needs the matching harmonic HDF5 and cell
metadata in ``output.result-dir``. After the calculation, copy these two files
from the same input set before plotting:

.. code-block:: bash

   cp inputs/fc2.hdf5 inputs/phono3py_disp.yaml results/wigner-4ph/
   nepkappa info wigner-4ph.yaml --for plot
   nepkappa plot wigner-4ph.yaml
   ls results/wigner-4ph/plots

The conductivity panel shows 3ph+4ph particle, coherence and total terms.
They are exported as ``kappa4-rta.dat``, ``kappa4-rta-coherence.dat`` and
``kappa4-rta-total.dat``. The 3ph and 4ph scattering rates share one panel,
so ``combined.png`` again has six panels. A coherence curve close to zero on
this scale should be checked in the numeric data before being described as zero.
