Software overview
=====================

Capabilities
----------------

NEP-kappa organizes crystal structures, atomic forces, force constants and result
analysis into one workflow. It supports NEP, VASP, MACE and other ASE calculators.
Phono3py finite displacements or HiPhive fitting generate FC2/FC3; Thirdorder
and Fourthorder provide the force constants used by FourPhonon.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Calculation
     - Method and outputs
   * - Harmonic phonons
     - Dispersion, density of states, volumetric heat capacity and group velocity.
   * - Temperature-dependent phonons
     - QHA thermal expansion, fixed-volume Phonopy SSCHA, and SSCHA at QHA equilibrium volumes.
   * - Three-phonon transport
     - phono3py RTA, LBTE and SMM19 Wigner particle/coherence contributions.
   * - Three- and four-phonon transport
     - FourPhonon RTA or iterative LBTE; Wigner transport with the Wigner_Park executable.
   * - Nonequilibrium relaxation
     - Time-dependent BTE for spatially homogeneous phonon populations at fixed frequencies.

Workflow
------------

A typical transport calculation follows this sequence::

   Structure → relaxation → FC2/FC3 → thermal transport → plots and report

Harmonic calculations stop after FC2. Four-phonon calculations additionally
require FC4. QHA and SSCHA introduce volume and temperature-dependent stages.
TD-BTE uses FC2/FC3 to evolve phonon occupations without requiring a preceding
conductivity calculation. Each route is developed in its own manual chapter.

The full ``run`` command orders the required stages from the input settings.
Individual commands can reuse completed stages. Force caches are reused when
their recorded inputs match. Long workstation calculations can run in tmux;
cluster calculations can use Slurm.

Input files
---------------

A complete phonon/transport input has six sections: ``structure``, ``calculator``,
``force-constant``, ``kappa``, ``plot`` and ``output``. TD-BTE replaces ``kappa``
with ``tdbte``. Optional execution resources belong in ``parallel``.
A direct stage command accepts a reduced input containing only its required
settings. See :doc:`input_files` for keys, units and examples.

Common commands
-------------------

.. list-table::
   :header-rows: 1
   :widths: 38 62

   * - Command
     - Purpose
   * - ``nepkappa info input.yaml``
     - Check the full input; add ``--for <command>`` for a single stage.
   * - ``nepkappa run input.yaml``
     - Execute the input-selected calculation sequence.
   * - ``nepkappa relax input.yaml``
     - Optimize the structure.
   * - ``nepkappa fc2 input.yaml``
     - Generate harmonic force constants or perform an enabled QHA scan.
   * - ``nepkappa fc2fc3 input.yaml``
     - Generate second- and third-order force constants.
   * - ``nepkappa kappa input.yaml``
     - Calculate thermal conductivity from prepared inputs.
   * - ``nepkappa tdbte input.yaml``
     - Evolve nonequilibrium occupations using existing force constants or a kernel.
   * - ``nepkappa plot input.yaml``
     - Plot saved results.
   * - ``nepkappa report input.yaml``
     - Write Markdown and YAML summaries.

Calculation commands show compact progress by default. Add ``-v`` for detailed
screen output; detailed diagnostics are saved in the result directory's
``run.log`` in either mode.

Related publications
-----------------------

S. Liu et al.,
*NEP-kappa 2.0: A unified framework for phonon transport and nonequilibrium thermalization with machine-learning interatomic potentials*,
arXiv preprint (2026),
`arXiv:2610.25477 <https://arxiv.org/abs/2610.25477>`_,
`doi:10.48550/arXiv.2610.25477 <https://doi.org/10.48550/arXiv.2610.25477>`_.

F. Yin, S. Liu, Y. Dong, A. A. Barinov, K. Xu and V. I. Khvesyuk,
*Accelerated phonon transport calculations for nanostructures: Combining
neuroevolution potentials and compressed sensing*, Journal of Applied Physics
**139**, 135103 (2026), `doi:10.1063/5.0324012 <https://doi.org/10.1063/5.0324012>`_.
See :doc:`reference` for BibTeX and model-specific references.

Related links
-----------------

* `NEP-kappa source code <https://github.com/lyushisyan/NEP-kappa>`_
* `NEP-kappa online documentation <https://nep-kappa.readthedocs.io/en/latest/>`_
* `Phonopy documentation <https://phonopy.github.io/phonopy/>`_
* `phono3py documentation <https://phonopy.github.io/phono3py/>`_
* `FourPhonon source code <https://github.com/FourPhonon/FourPhonon>`_
