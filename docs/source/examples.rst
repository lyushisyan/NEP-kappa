Choose an example
===================

The public catalog contains nine six-section static 3C-SiC inputs and one
separate time-dependent BTE input. Run from the repository root. The YAML
files are workflow demonstrations, not converged material results. See
``examples/README.md`` for prerequisites and resource guidance.

.. list-table::
   :header-rows: 1
   :widths: 38 62

   * - YAML file
     - Use
   * - ``vasp-rta-3ph.yaml``
     - VASP force calculations and phono3py three-phonon RTA.
   * - ``nep-rta-wigner-3ph.yaml``
     - NEP, phono3py three-phonon RTA, and SMM19 Wigner transport.
   * - ``nep-lbte-wigner-3ph.yaml``
     - NEP, phono3py three-phonon LBTE, and SMM19 Wigner transport.
   * - ``nep-rta-3ph-4ph.yaml``
     - NEP and FourPhonon combined three-plus-four-phonon RTA.
   * - ``nep-rta-wigner-3ph-4ph.yaml``
     - NEP and FourPhonon Wigner_Park 3ph/4ph RTA population and coherence.
   * - ``nep-lbte-3ph-rta-4ph.yaml``
     - NEP, iterative 3ph LBTE with 4ph RTA scattering.
   * - ``nep-lbte-3ph-lbte-4ph.yaml``
     - NEP, iterative LBTE for both 3ph and 4ph scattering.
   * - ``nep-qha-rta-3ph.yaml``
     - NEP, new FC2/FC3 at each QHA equilibrium volume, then 3ph RTA.
   * - ``nep-qha-sscha-rta-3ph.yaml``
     - NEP, SSCHA FC2 at QHA volume and temperature, then 3ph RTA.
   * - ``tdbte.yaml``
     - Population dynamics from existing matching FC2/FC3 and metadata.

Use ``nepkappa info examples/<name>.yaml`` to check an input without
running it. Ordinary workflow inputs use ``nepkappa run``; plotting and
reporting are separate commands. The TD-BTE route needs completed force
constants and does not regenerate them.

Three-phonon-only RTA and LBTE use ``kappa.method``. FourPhonon inputs set
``kappa.method-3ph`` and ``kappa.method-4ph`` independently. Supported pairs
are RTA/RTA, LBTE/RTA, and LBTE/LBTE; Wigner_Park requires RTA/RTA. Standalone
SSCHA is not part of this catalog.

The six-section layout also supports individual stage commands. For example,
``nepkappa stage fc2 examples/vasp-rta-3ph.yaml`` generates only FC2, and
``nepkappa stage fc2fc3 examples/vasp-rta-3ph.yaml`` generates FC2 and FC3.
Run ``nepkappa stage relax`` first when the input requests relaxation.

External VASP, Fourthorder, and FourPhonon installations and site-specific
Slurm settings must be configured before execution. The calculation can be
submitted through an external batch script, or optional ``parallel`` settings
can let NEP-kappa create stage-specific jobs.
