Choose an example
===================

Use one template per workflow. All public examples live in ``examples/``;
run from the repository root, and copy an input before adapting it.

.. list-table::
   :header-rows: 1
   :widths: 36 64

   * - YAML file
     - Use
   * - ``nep-rta.yaml``
     - Bulk Si NEP RTA; set ``kappa.method: lbte`` for serial LBTE.
   * - ``vasp-rta.yaml``
     - VASP relaxation and FC2/FC3; configure executable and POTCAR paths.
   * - ``mace-rta.yaml``
     - MACE on CPU; install the optional backend and set a checkpoint.
   * - ``wigner.yaml``
     - phono3py SMM19 Wigner transport.
   * - ``qha.yaml``
     - Isotropic volume-scan QHA.
   * - ``sscha.yaml``
     - Phonopy stochastic SSCHA and temperature-dependent transport.
   * - ``qha-sscha.yaml``
     - SSCHA at QHA volumes, with independent three/four-phonon switches.
   * - ``bas-three-four-phonon.yaml``
     - Both pure three-phonon and combined three-plus-four-phonon conductivity.
   * - ``nep-hiphive.yaml``
     - FC2/FC3 fitting with HiPhive.
   * - ``film.yaml``
     - Film geometry with explicit effective thickness.
   * - ``thirdorder.yaml``
     - FC2/FC3 generation for ShengBTE; stage-only workflow.
   * - ``vasp-fc4.yaml``
     - VASP/Fourthorder FC4 generation; stage-only workflow.
   * - ``slurm-vasp.yaml``
     - VASP force arrays and workflow continuation.
   * - ``slurm-lbte.yaml``
     - Distributed phono3py LBTE.
   * - ``fourphonon.yaml``
     - FourPhonon workflow with Slurm transport.
   * - ``compare.yaml``
     - Two or more completed model results; use ``nepkappa compare``.
   * - ``converge-qmesh.yaml``
     - q-mesh study; use ``nepkappa converge``.

All ordinary workflow inputs use ``nepkappa run examples/<name>.yaml``.
For pre-existing force constants, use ``kappa`` or ``kappa4`` instead.
Inspect stage-only configurations with ``validate --for <stage>``.
Comparison and convergence inputs have separate schemas and are not accepted
by the workflow ``validate`` command.

Small variants
----------------

RTA and serial LBTE differ by ``kappa.method``. VASP/NEP and finite
displacement/HiPhive are independent choices: combine the relevant calculator
and fitting sections rather than maintaining every combination as a full file.
HiPhive cutoffs must fit the actual structure and supercell. A film requires
its physical thickness and appropriate non-periodic mesh direction.

The previous numbered filenames have been replaced by descriptive names.
Duplicate RTA/LBTE, VASP/HiPhive, film/HiPhive, VASP/Thirdorder, and two-model
comparison variants were consolidated. Parser compatibility is retained in
tests where needed.

Slurm examples start with ``submit: false`` and require your environment and
resource settings. Earlier local stages in a complete workflow can still run.
For a configuration-only check, use ``validate`` or ``info``.
