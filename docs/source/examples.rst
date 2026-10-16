Example catalog
===================

The public catalog contains nine six-section phonon/transport 3C-SiC inputs and one
six-section time-dependent BTE input. Run from the repository root. The YAML
files are workflow demonstrations, not converged material results. See
``examples/README.md`` for prerequisites and resource guidance.

.. list-table::
   :header-rows: 1
   :widths: 38 62

   * - YAML file
     - Use
   * - :download:`vasp-rta-3ph.yaml <../../examples/vasp-rta-3ph.yaml>`
     - VASP force calculations and phono3py three-phonon RTA.
   * - :download:`nep-rta-wigner-3ph.yaml <../../examples/nep-rta-wigner-3ph.yaml>`
     - NEP, phono3py three-phonon RTA, and SMM19 Wigner transport.
   * - :download:`nep-lbte-wigner-3ph.yaml <../../examples/nep-lbte-wigner-3ph.yaml>`
     - NEP, phono3py three-phonon LBTE, and SMM19 Wigner transport.
   * - :download:`nep-rta-3ph-4ph.yaml <../../examples/nep-rta-3ph-4ph.yaml>`
     - NEP and FourPhonon combined three-plus-four-phonon RTA.
   * - :download:`nep-rta-wigner-3ph-4ph.yaml <../../examples/nep-rta-wigner-3ph-4ph.yaml>`
     - NEP and FourPhonon Wigner_Park 3ph/4ph RTA population and coherence.
   * - :download:`nep-lbte-3ph-rta-4ph.yaml <../../examples/nep-lbte-3ph-rta-4ph.yaml>`
     - NEP, iterative 3ph LBTE with 4ph RTA scattering.
   * - :download:`nep-lbte-3ph-lbte-4ph.yaml <../../examples/nep-lbte-3ph-lbte-4ph.yaml>`
     - NEP, iterative LBTE for both 3ph and 4ph scattering.
   * - :download:`nep-qha-rta-3ph.yaml <../../examples/nep-qha-rta-3ph.yaml>`
     - NEP, new FC2/FC3 at each QHA equilibrium volume, then 3ph RTA.
   * - :download:`nep-qha-sscha-rta-3ph.yaml <../../examples/nep-qha-sscha-rta-3ph.yaml>`
     - NEP, SSCHA FC2 at QHA volume and temperature, then 3ph RTA.
   * - :download:`tdbte.yaml <../../examples/tdbte.yaml>`
     - Generate NEP FC2/FC3 and matching metadata, then evolve phonons in time.

Use ``nepkappa info examples/<name>.yaml`` to check an input without
running it. Ordinary workflow inputs use ``nepkappa run``; plotting and
reporting are separate commands. The TD-BTE route needs completed matching
FC2/FC3 and displacement metadata, but not a prior conductivity
result. The direct ``tdbte`` command reuses these files, while ``run``
generates them before time evolution.

Three-phonon-only RTA and LBTE use ``kappa.method``. FourPhonon inputs set
``kappa.method-3ph`` and ``kappa.method-4ph`` independently. Supported pairs
are RTA/RTA, LBTE/RTA, and LBTE/LBTE; Wigner_Park requires RTA/RTA. Additional downloadable harmonic, fixed-volume SSCHA and artifact-reuse
inputs are provided in the corresponding manual chapters.

The six-section layout also supports individual calculation commands. For example,
``nepkappa fc2 examples/vasp-rta-3ph.yaml`` generates only FC2, and
``nepkappa fc2fc3 examples/vasp-rta-3ph.yaml`` generates FC2 and FC3.
Run ``nepkappa relax`` first when the input requests relaxation.

External VASP, Fourthorder, and FourPhonon installations and site-specific
Slurm settings must be configured before execution. The calculation can be
submitted through an external batch script, or optional ``parallel`` settings
can let NEP-kappa create stage-specific jobs.

Manual chapter inputs
---------------------

The following inputs follow the numbered manual. Download them to each
chapter's working directory, not the repository root. They are displayed from
the same YAML files used for the downloads.

.. list-table::
   :header-rows: 1
   :widths: 20 45 35

   * - Chapter
     - Download
     - Calculation command
   * - :doc:`starting`
     - :download:`quickstart.yaml <inputs/quickstart.yaml>`
     - ``run``
   * - :doc:`harmonic_temperature`
     - :download:`harmonic.yaml <inputs/harmonic.yaml>`
     - ``relax``, then ``fc2``
   * - QHA
     - :download:`qha.yaml <inputs/qha.yaml>`
     - ``fc2``
   * - Fixed-volume SSCHA
     - :download:`sscha.yaml <inputs/sscha.yaml>`
     - ``run``
   * - QHA+SSCHA
     - :download:`qha-sscha.yaml <inputs/qha-sscha.yaml>`
     - ``run``
   * - :doc:`thermal_transport`
     - :download:`lbte.yaml <inputs/lbte.yaml>`
     - ``kappa``
   * - Four-phonon comparison
     - :download:`fourphonon.yaml <inputs/fourphonon.yaml>`
     - ``run``
   * - Wigner with 3ph+4ph
     - :download:`wigner-4ph.yaml <inputs/wigner-4ph.yaml>`
     - ``kappa``
   * - :doc:`tdbte`
     - :download:`tdbte.yaml <inputs/tdbte.yaml>`
     - ``tdbte``

Generic FourPhonon templates use the default phono3py FC3 backend, a 12×12×12
mesh and 300 K. The chapter comparison instead uses Thirdorder FC3, an 8×8×8
mesh, three temperatures and ``compare-3ph: true``. Generic QHA/SSCHA templates
include transport; Chapter 5 focuses on phonon properties. Generic TD-BTE
creates force constants and uses a 5×5×5 dynamics grid; Chapter 7 reuses force
constants on a 6×6×6 grid. Use the chapter inputs to follow the manual's workflow.
