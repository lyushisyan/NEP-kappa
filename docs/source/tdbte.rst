Experimental time-dependent BTE
========================================

The ``tdbte`` stage evolves spatially homogeneous phonon occupations at fixed
frequencies using an entropy-interpolated three-phonon energy-shell operator.
Physical relaxation rates have not been independently validated.
It does not model laser absorption, electrons, coherent phonons, spatial drift,
four-phonon scattering, or frequency renormalization during time evolution.

Running TD-BTE
----------------

Supply a prebuilt energy-shell kernel and its JSON metadata. Automatic kernel
construction from FC2/FC3 is not implemented, and a phono3py LBTE matrix is not
an accepted input.

Use ``examples/tdbte.yaml``. From the repository root::

   nepkappa validate examples/tdbte.yaml --for tdbte
   nepkappa info examples/tdbte.yaml --for tdbte
   nepkappa run examples/tdbte.yaml
   nepkappa report calculations/3C-SiC/tdbte/results

The direct command ``nepkappa tdbte input.yaml`` runs the same stage. The example
has a custom plan containing only ``tdbte`` and does not generate forces or load
a calculator. Supply a real kernel before execution. Validation checks options,
not kernel existence; execution loads and audits the artifact. Relative paths
resolve from the launch directory, not the YAML directory.

Input format
--------------

* ``experimental: true`` is mandatory.
* ``temperature`` is the initial equilibrium temperature in kelvin.
* ``branches`` explicitly lists zero-based branch IDs to excite at every q.
* ``excitation`` multiplies selected occupations by ``1 + excitation``; it is
  not a temperature increment, laser fluence, or equal energy per mode.
* ``duration-ps``, ``max-step-ps`` and ``samples`` control integration/output.
  A small max step is important for the equilibrium control of this explicit
  solver. Tolerances are currently fixed at rtol=1e-11, atol=1e-14.

``kernel.npz`` contains numeric arrays (object arrays are not accepted):

* ``frequency_THz``: shape (M,), nonnegative frequencies;
* ``zero_modes``: boolean (M,), only fixed zero-energy translations;
* ``parent``: integer (N,);
* ``daughter1``, ``daughter2``: integer (N,4);
* ``weights``: nonnegative (N,4), each row sums to one;
* ``coefficients``: nonnegative (N,), including all integration, mesh,
  permutation and physical-unit factors, giving dn/dt in ps^-1.

The adjacent ``kernel.json`` must contain ``model`` (provenance description),
``coefficient_time_unit: "ps"``, ``mode_branch_indices`` (M integers), and
``uniform_mode_weight`` (1/Nq for a full uniform grid). Symmetry-reduced grids
with unequal weights are unsupported. Include source hashes and generation
settings in the sidecar for reproducibility. The loader checks exact discrete
energy-shell consistency, but cannot independently establish the correctness of
the supplied physical prefactor or provenance claims.

Numerics and diagnostics
----------------------------------------

For every event, the same barycentric weights interpolate frequency and
entropy variable eta=log[(1+n)/n] and deposit the daughter occupation updates.
This enforces discrete energy conservation and Bose detailed balance. The
nonlinear forward-minus-reverse flux is evaluated in log form to avoid
intermediate overflow. Occupations are not silently clipped.

Outputs are written to ``result-dir/tdbte/``. The stage refuses to overwrite an
existing directory. ``trajectories.npz`` contains the excited trajectory and an
independently integrated unexcited control. ``audit.json`` records numerical
checks, kernel/metadata hashes, tolerances, excitation and limitations.
``branch-energy.csv`` gives energy above the initial Bose state in eV per
primitive cell, including the explicitly supplied uniform mesh weight.
``report`` includes the numerical audit and experimental warning.

Passing the energy, equilibrium and sampled entropy checks is necessary but
not sufficient for physical validation. A global error in the collision
prefactor can pass all of these tests while changing relaxation times.
Independent normalization, permutation-counting and linearized-operator
benchmarks remain required. Plotting interpolation is not part of the solver;
the manuscript Figure6 plotting scripts remain separate research tools.
