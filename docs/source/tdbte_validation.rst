TD-BTE collision-kernel validation gate
=======================================

Status
---------

The independent 2026-11-07 audit passes discrete numerical invariants and a
matched-quadrature scattering-prefactor control, but does not validate the
production kernel's absolute rates. Its sampled off-shell small-q discrepancy
is an auxiliary diagnostic, not by itself a failure of on-shell dynamics.
Following review, this check is optional and does not gate software development.
Stages 2--4 are implemented in :doc:`tdbte_spatial`, with separate numerical
and physical-validation status. The homogeneous integration and all physical
collision coefficients remain unchanged.

The audit is implemented in ``tools/validate_tdbte_kernel.py``. It writes
checksummed production kernels, mode-ordered rate arrays and a JSON report
under a fresh output directory. It uses the existing nonlinear collision RHS;
it introduces no RTA, Callaway or linearized evolution solver. The independent
Jacobian action is a derivative check only.

Rate conventions and normalization
-------------------------------------

Phono3py's public ``gamma`` is the imaginary self-energy (half linewidth),
reported in THz of ordinary frequency. For positive calculated values,

.. math::

   \tau_\lambda^{-1}\ [\mathrm{ps}^{-1}]
   = 2(2\pi)\,\mathrm{gamma}_\lambda\ [\mathrm{THz}].

The factors two and :math:`2\pi` have different origins: half linewidth versus
population lifetime, and ordinary versus angular frequency. See the
`phono3py output definition <https://phonopy.github.io/phono3py/input-output-files.html#gamma>`_.

In the production builder, native interaction strengths already include
:math:`1/N_q`. Both daughter orderings are retained, and an event coefficient is

.. math::

   c_e=4\pi U\,S_e\,w_e^{\mathrm{surface}},

where :math:`U` is the native imaginary-self-energy conversion factor.
Surface weights integrate a delta function of frequency in THz, in mesh-cell
coordinates. The native mesh normalization must not be applied again. The
separate :math:`1/N_q` in reported energy is the observable's Brillouin-zone
average. No fitted scale or extra daughter factor is applied.

The audit checks the full unshifted grid ordering, maps the native
:math:`q_0+q_1+q_2=G` convention to decay parent :math:`q_p=-q_0`, and verifies
the recorded addresses and frequencies against a separately loaded native
reference. Modes are indexed as ``q_index * number_of_branches + branch``;
frequencies alone are never used to match modes. Degenerate linewidths are
averaged using the same reporting convention as the public native accessor.

Three different rate observables must be kept separate:

* A tagged-leg linewidth holds the other two legs fixed. If two daughters are
  the same mode, it counts each tagged occurrence separately.
* The diagonal of the occupation Jacobian combines all occurrences of a mode
  before differentiating. Repeated daughters produce cross terms. With
  interpolated events, deposition and interpolation weights also enter the
  diagonal; it is not directly a native linewidth.
* The parent-decay loss coefficient uses only the parent leg of
  :math:`p\leftrightarrow b+c`, :math:`\sum_e c_e(1+n_b^*+n_c^*)`.
  The corresponding independent native reference is the class-2 decay
  contribution, converted by :math:`4\pi`. This tests one identifiable channel,
  not the complete nonlinear operator's absolute rates.

The prefactor control reconstructs Gaussian-weighted tagged events and compares
them to an independent ``ImagSelfEnergy`` call with identical FCs, mesh,
temperature, sigma, cutoff and physical constants. These Gaussian events are
used only for this normalization diagnostic; they are never propagated and
never substituted for production energy-shell events.

The production solver uses SciPy physical constants, while the installed native
backend uses its own constants. The report records both. This small difference
is retained in the production-shell comparison; the matched Gaussian control
uses native constants to isolate prefactors without changing the solver.

Numerical and physical checks
--------------------------------

The packaged entropy-interpolated shell operator is a weak discretization.
With event incidence vectors :math:`B_e`, equilibrium reverse fluxes
:math:`R_e` and :math:`D_\lambda=n_\lambda(1+n_\lambda)`, its derivative is
:math:`J=-BRB^T D^{-1}`. The audit combines repeated mode indices, assembles
the action in bounded event blocks, and compares it with central differences
of the actual nonlinear collision RHS.

It also checks Bose equilibria at 100, 300 and 1000 K, random positive
nonequilibrium occupations at three amplitudes, energy conservation, entropy
production and a finite-amplitude response involving multiple interacting
legs. These identities do not validate an absolute rate: multiplying every
event coefficient by a constant would retain them.

Surface convergence is tested by subdividing the existing degree-two surface
triangles into four and reconstructing the same linear strength and barycentric
interpolants. Triangle groups split between disk chunks are reassembled without
loading the complete kernel. This isolates surface-integration error; it does
not refine the phonon grid, improve branch tracking or change the acoustic
endpoint model.

At Gamma, the builder fixes three rigid translations to zero frequency,
zero entropy variable and zero vertex strength. Other unstable or near-zero
modes are rejected. The direct small-q diagnostic evaluates native Fourier and
normal-mode vertices at explicit non-Gamma triplets. It first verifies their
branch-summed strength against a native grid triplet, then tests two soft-q
directions and all six leg permutations. Branch sums avoid arbitrary gauges
inside degenerate subspaces. No source force constants are symmetrized or
altered by the diagnostic.

The sampled quantity is off shell: it has no energy delta or Bose factors and
is **not a scattering rate**. A growing small-q strength does not establish a
divergent on-shell linewidth or an error in native steady-state transport.
It establishes that these sampled native vertices do not support the
production interpolant's zero-endpoint assumption. This unresolved physical
approximation remains a limitation to investigate in on-shell rates and
convergence; the off-shell check alone does not block the spatial extension.

Reproduce
------------

From the repository root, with the installed project and an existing matching
FC2/FC3 source directory and a local copy of the audit script. The ``tools/``
directory is ignored by Git and is not distributed in a fresh checkout:

.. code-block:: bash

   OMP_NUM_THREADS=1 RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
     python tools/validate_tdbte_kernel.py \
     --source calculations/3C-SiC/manual-stages \
     --output calculations/3C-SiC/tdbte-validation-reproduction \
     --meshes 2 3 4 --surface-levels 0 1 --acoustic-check

The output directory must not exist. Source paths are relative to the launch
directory. No calculator or force jobs are started. The report records hashes
of ``phono3py_disp.yaml``, ``fc2.hdf5`` and ``fc3.hdf5`` and verifies them again
after calculations. Raw mode-ordered rates are in ``mesh-N/rates.npz``.
The present diagnostic deliberately exits with code 1 and
``stage_1_acceptance_passed: false``: the complete production-shell total-rate
and q-mesh validation is still open, even if individual controls pass.
Omit ``--acoustic-check`` to skip the off-shell diagnostic. Its outcome is
recorded separately and does not prevent spatial implementation.

The local source used for the recorded audit is a relaxed two-atom 3C-SiC
primitive cell, from the public material-specific NEP model, with 3x3x3 FC2/FC3
supercells. It contains no NAC data. Its maximum checked FC2 atom-sum residual
is :math:`1.62\times10^{-14}` eV/Angstrom squared, and the checked FC3 residuals
are :math:`3.59\times10^{-14}` eV/Angstrom cubed. These compact-array partial
checks do not prove every symmetry or off-grid interpolation property.

The source is a local calculation artifact, not included in the public
checkout. A different FC dataset may be supplied, but its results constitute
a new benchmark; do not expect the recorded numerical values or hashes to
match. Regenerating FCs from ``examples/tdbte.yaml`` uses that example's
unrelaxed structure and is also a different source.

Recorded results (2026-11-07)
--------------------------------

The audit used phono3py/Phonopy 4.1.0, Python 3.12 and Numba reductions.
Source SHA256 identities are:

.. code-block:: text

   phono3py_disp.yaml  9b14f2aaabf8eb8ed964683d01c571190321395b1668e6da937e4d85ba17e5ba
   fc2.hdf5            e6b89c8395e6d87a88998071f3976f16becca7258784a18bd7f1d0a0d7530f17
   fc3.hdf5            3ec35cadd98d33d9f7e8db80fcbe78893a3746c65d3d97210b47771257f345a6

.. list-table:: Rates and surface refinement at 300 K
   :header-rows: 1
   :widths: 15 28 28 29

   * - q mesh
     - Matched Gaussian total tagged-rate relative L2 error
     - Production parent-decay relative L2 difference from native tetrahedra
     - Relative change after one surface subdivision
   * - 2x2x2
     - 1.18e-15
     - 0.2375
     - 1.18e-4
   * - 3x3x3
     - 4.62e-16
     - 0.04362
     - 1.35e-6
   * - 4x4x4
     - 1.04e-15
     - 0.02612
     - 1.47e-6

These are different integration schemes for the production comparison;
the percentages are discrepancies, not an inferred correction factor.
Surface convergence of the chosen interpolant does not establish agreement
with the native tetrahedron operator or q-mesh convergence. For example, the
Gamma optical parent-decay loss coefficient changes from 0.24827 to 0.50178
to 0.54936 inverse ps as the mesh increases. The corresponding native values
are 0.21230, 0.47695 and 0.53378 inverse ps. These coarse meshes do not establish
a converged absolute rate.

All sampled equilibrium, energy and entropy checks pass. Central-difference
Jacobian-action errors reach approximately 3e-11, and multi-leg finite-amplitude
perturbations show a nonzero departure from a simply doubled collision RHS.

An actual homogeneous DOP853 integration also used the new 3x3x3 kernel at
300 K, with 1% excitation of branches 3--5, 0.1 ps duration, 0.01 ps maximum
step and 11 saved samples. The input passed ``nepkappa info --for tdbte``.
The unmodified solver's audit passed positivity/finite, equilibrium, energy and
sampled entropy checks: maximum energy drift was 5.65e-15 of injected energy,
with zero measured equilibrium redistribution and entropy decrease. This is
a short numerical smoke calculation; it does not establish time-step
convergence or physical accuracy. Its local input and outputs are
``tdbte-validation-20261107-final/homogeneous-smoke.yaml`` and
``tdbte-validation-20261107-final/homogeneous-smoke/tdbte/`` under the same
``calculations/3C-SiC/`` directory.

For the parent (0.2, 0.1, 0.15), a soft leg along (epsilon, 0, 0) gives
branch-summed raw strengths of 3.86e-5, 3.79e-4 and 3.79e-3 as epsilon goes
from 1e-3 to 1e-4 to 1e-5. The log slope is -0.996; the second direction gives
-0.995. All sampled frequencies remain positive. The six leg permutations
agree within 7e-13 relative, and the direct/native grid-triplet strength
reference agrees within 2e-16. The zero-endpoint assumption therefore remains
unsupported even though the checked FC sum residuals are small.

The historical Stage 1 report (recorded before the decision to make the
off-shell test optional) is available as
:download:`validation JSON <_static/tdbte-validation-20261107.json>`.
The large kernels and raw per-mode arrays remain in the local
``calculations/3C-SiC/tdbte-validation-20261107-final/`` directory.

Run the synthetic diagnostic tests separately:

.. code-block:: bash

   OMP_NUM_THREADS=1 RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
     python -m pytest tests/public/test_tdbte_validation.py -q

They verify repeated-leg derivatives, analytic surface moments and triangle
reconstruction across chunk boundaries. They do not replace the physical audit.

The complete repository test command with all three thread controls above
passed 846 cases (including 132 public cases), with 24 warnings and no failures
or skips. The six new public cases test audit mathematics only. The strict
Sphinx HTML build also passed. The physical-audit command's exit code 1 is the
expected incomplete physical-rate acceptance and is separate from these software results.

Remaining gate
-----------------

Establish complete tagged-rate/operator benchmarks and q-mesh, surface and
force-supercell convergence, including the actual on-shell weight of acoustic
channels. Investigate endpoint interpolation where it affects those integrals.
Agreement in a
matched Gaussian prefactor control or a single decay channel is insufficient.
Do not rescale coefficients, delete acoustic events, disable averaging, or
alter source FCs merely to force agreement. These are requirements for
quantitative physical claims; software development and numerical validation of
the spatial solver proceed independently.
