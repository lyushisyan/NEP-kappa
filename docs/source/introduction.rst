Overview
========================

NEP-kappa uses a single YAML configuration to connect atomic-force calculations, force-constant extraction, phonon properties and thermal transport.

- **Calculation backends**: NEP, VASP, MACE and other ASE-compatible calculators.
- **Integrated phonon tools**: Phonopy, phono3py and FourPhonon for phonon properties, anharmonic scattering and lattice thermal conductivity.
- **Nonequilibrium phonon dynamics**: time-dependent BTE simulations of phonon populations, energy relaxation and thermalization.

Calculation workflow
------------------------

.. graphviz::
   :align: center
   :caption: NEP-kappa calculation workflow

   digraph NEPkappa {
       graph [rankdir=TB, bgcolor="transparent", pad="0.20",
              nodesep="0.35", ranksep="0.48", splines=polyline];

       node [shape=box, style="rounded,filled",
             color="#A7B8CF", fillcolor="#F3F6FA",
             fontname="sans-serif", fontsize=12,
             margin="0.18,0.12"];

       edge [color="#64748B", penwidth=1.3, arrowsize=0.65];

       structure [label="Structure relaxation"];
       forces [label="Energies and atomic forces\nNEP / VASP / MACE / ASE"];
       ifcs [label="Force-constant generation\nFC2 / FC3 / FC4"];

       structure -> forces -> ifcs;

       thermal [label="Phonon properties and thermal transport\nPhonopy / Phono3py / ShengBTE / FourPhonon", fillcolor="#DFEAF8"];
       td [label="Nonequilibrium phonon dynamics\nNEP-kappa TD-BTE", fillcolor="#DFEAF8"];

       { rank=same; thermal; td; }

       ifcs -> thermal;
       ifcs -> td;

       qha [label="QHA / SSCHA\n(optional)",
            style="rounded,dashed,filled", fillcolor="#FFFFFF"];
       phonon [label="Dispersion / DOS / group velocities / heat capacity"];
       kappa [label="Three- / four-phonon scattering\nConductivity (RTA / LBTE / Wigner)"];
       thermal_result [label="Save results and plot", fillcolor="#E9F4EF"];

       thermal -> qha [style=dashed];
       qha -> phonon [style=dashed];
       thermal -> phonon;
       phonon -> kappa -> thermal_result;

       excitation [label="Initial populations and excitation"];
       kernel [label="Phonon collision kernel and integration"];
       solver [label="Solve time-dependent BTE\nPopulation evolution / energy thermalization"];
       td_result [label="Save results and plot",
                  fillcolor="#E9F4EF"];

       td -> excitation -> kernel -> solver -> td_result;

       { rank=same; qha; excitation; }
       { rank=same; phonon; kernel; }
       { rank=same; kappa; solver; }
       { rank=same; thermal_result; td_result; }
   }

Common commands
------------------------

Replace ``input.yaml`` below with the path to your input configuration.
Typically, check the configuration with ``info``, then use ``run`` to execute the enabled stages.

.. list-table::
   :header-rows: 1
   :widths: 14 36 46
   :width: 100%

   * - Category
     - Command
     - Purpose
   * - Input check
     - ``nepkappa info input.yaml``
     - Validate and display settings without starting calculations.
   * - Workflow
     - ``nepkappa run input.yaml``
     - Execute the configured stages in order.
   * - Relaxation
     - ``nepkappa relax input.yaml``
     - Relax the crystal structure.
   * - Force constants
     - ``nepkappa fc2 input.yaml``
     - Generate second-order force constants.
   * - Force constants
     - ``nepkappa fc2fc3 input.yaml``
     - Generate second- and third-order force constants.
   * - Force constants
     - ``nepkappa fc2fc3fc4 input.yaml``
     - Generate second-, third- and fourth-order force constants in order.
   * - Thermal transport
     - ``nepkappa kappa input.yaml``
     - Calculate conductivity using the configured transport route.
   * - Nonequilibrium dynamics
     - ``nepkappa tdbte input.yaml``
     - Solve the time-dependent BTE for phonon-population evolution.
   * - Post-processing
     - ``nepkappa plot input.yaml``
     - Read and plot existing results.
   * - Post-processing
     - ``nepkappa report input.yaml``
     - Summarize saved results in Markdown and YAML reports.
