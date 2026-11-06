Time-dependent phonon BTE and nonequilibrium relaxation
===========================================================

TD-BTE evolves spatially homogeneous phonon occupations at fixed frequencies
using three-phonon scattering. This chapter selectively excites the optical
branches of 3C-SiC and follows energy transfer to the acoustic branches.
The model has no spatial transport, external heat bath or continuing excitation.

Workflow and preparation
----------------------------

Reuse matching ``fc2.hdf5``, ``fc3.hdf5`` and ``phono3py_disp.yaml`` from the
quick-start calculation. The metadata supplies the structure, supercells and
primitive mapping. A conductivity result is not needed.

.. code-block:: bash

   mkdir -p ~/work/3C-SiC/tdbte
   cd ~/work/3C-SiC/tdbte
   ls -lh ../results/fc2.hdf5 ../results/fc3.hdf5 ../results/phono3py_disp.yaml

The calculation constructs a three-phonon collision kernel on the full q mesh,
then integrates the occupations. It stores kernel data in disk-backed chunks.
For a calculation starting from atomic forces instead, the six-section
``examples/tdbte.yaml`` uses ``nepkappa run`` to prepare FC2/FC3 first.

Optical excitation and energy relaxation
--------------------------------------------

Input and execution
~~~~~~~~~~~~~~~~~~~~~~~

Save :download:`tdbte.yaml <inputs/tdbte.yaml>` in the working directory:

.. literalinclude:: inputs/tdbte.yaml
   :language: yaml

The 6×6×6 grid has 216 q points. In two-atom 3C-SiC, branches 0–2 are acoustic
and branches 3–5 are optical. Branch indices are zero-based and frequency-ordered
at each q point. ``excitation: 0.01`` increases the selected initial Bose
occupations by 1% relative to 300 K. ``samples: 201`` sets the saved output
times over 20 ps; the adaptive integrator uses ``max-step-ps`` as its step limit.

.. code-block:: bash

   nepkappa info tdbte.yaml --for tdbte
   nepkappa tdbte tdbte.yaml

Use a :ref:`tmux session <transport-tmux>` named ``tdbte`` for a long remote
workstation run. Progress first reports kernel construction and then integration
time. After completion:

.. code-block:: bash

   tail -n 20 results/optical-m666/run.log
   ls -lh results/optical-m666/tdbte

The output directory contains ``trajectories.npz``, ``branch-energy.csv`` and
``audit.json``. Inspect the audit alongside the completion message. Numerical
conservation and equilibrium checks are distinct from validation of absolute
physical relaxation rates; the latter remains incomplete for this implementation.

Plot and interpret
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   nepkappa info tdbte.yaml --for plot
   nepkappa plot tdbte.yaml
   ls results/optical-m666/plots

The combined figure has six panels:

* Top row: dispersion at 0, 10 and 20 ps, colored by excess energy per mode
  relative to the initial 300 K equilibrium distribution. All three share one
  color scale; frequencies remain fixed.
* Bottom left: excess energy of all six branches, in eV per primitive cell.
* Bottom middle: energy-equivalent temperatures of the three acoustic branches.
* Bottom right: energy-equivalent temperatures of the three optical branches.

The mode color is ``h*f*(n(t)-n_Bose(T0))`` in eV. Periodic trilinear
interpolation maps the q-mesh energies onto the band path; it does not calculate
dynamics at additional q points. Branch energies include the uniform 1/Nq
weight. Branch temperatures are obtained by matching each branch's thermal
energy to a Bose distribution at its own frequencies; they are not fits to
mode occupations and do not imply equilibrium within each branch.

In the manual example, optical excess energy decreases as acoustic branches
gain energy. The branch temperatures remain unequal at 20 ps, so this time
window does not establish complete equilibration. With no heat bath, relaxation
does not require the system to return to the initial 300 K state.

``layout: both`` writes ``tdbte_dispersion_1.png`` through
``tdbte_dispersion_3.png``, ``tdbte_branch_energy.png``,
``tdbte_acoustic_temperature.png``, ``tdbte_optical_temperature.png`` and
``tdbte_combined.png``. It also writes ``tdbte_branch_temperature.csv`` and
``tdbte_plot_metadata.json``. Snapshot times must be increasing and inside the
saved interval; the nearest saved times are used. By default, the start,
midpoint and end are selected.

Keep the original FC2 and metadata at their recorded source, or matching copies
in the result directory. The six-panel plot checks input hashes, frequencies,
the numerical audit and kernel metadata. Legacy CSV-only results support
energy plots but not the three dispersion snapshots.

Input choices
-----------------

* ``force-constants`` is the source directory, not a request to compute forces.
* ``mesh`` is the **q grid for the collision integral**, independent of
  ``kappa.mesh``. Supply three integers of at least 2 for 3D, or one integer
  of at least 2 and two ones for the experimental 1D route, e.g. ``[51,1,1]``.
  The builder uses a full unshifted diagonal grid, not an irreducible grid.
  Two-dimensional meshes are not supported.
* ``temperature`` (default 300 K) defines the initial Bose distribution.
* ``branches`` explicitly lists zero-based, frequency-ordered branch IDs at
  every q point. The number of branches is read from the structure. ``[3,4,5]``
  is not a universal optical-mode selection for other materials.
* ``excitation`` (default 0.01) multiplies occupations of those branches by
  ``1 + excitation``. This is a 1% occupation increase, not a 1 K increase,
  a specified laser fluence, or equal injected energy in every mode.
* ``duration-ps`` defaults to 20 ps; ``max-step-ps`` to 0.5 ps; ``samples`` to
  201 saved times. The explicit DOP853 integrator uses rtol=1e-11 and atol=1e-14.
  Check time-step convergence separately from q-grid convergence. ``samples``
  controls saved times, not the internal integration step.

Kernel construction and memory
----------------------------------

The builder evaluates phono3py native shortest-vector three-phonon interaction
strengths one parent q point at a time. It integrates linearly interpolated
energy shells over six Freudenthal tetrahedra per mesh cell using a positive,
symmetric degree-two surface rule. Ordered daughters use the coefficient
``4*pi*U*strength*surface_weight``, where ``U`` is phono3py's imaginary-self-energy
unit conversion. Native strengths already include ``1/Nq``; it is not applied
a second time in the collision integral.

The same barycentric weights interpolate frequencies and
eta=log[(1+n)/n] and deposit daughter occupation updates. This gives discrete
energy conservation and Bose detailed balance. It is not the native phono3py
tetrahedron collision operator and is not an RTA exponential decay.

Events are written in blocks of at most 20,000 to
``result-dir/tdbte-kernel/``. During each derivative evaluation, only one block
is memory-mapped, reduced and closed. No full event collection or dense
mode-by-mode collision matrix is loaded into memory. Numba is installed through
HiPhive's dependencies and is used automatically, without an additional
installation step. Compiled block reductions are compared against the NumPy
reference before integration. If Numba is unavailable, the solver uses NumPy.
The audit records the backend; acceleration does not change the input format
or physical model.

Chunking does not remove the cost of dense q meshes. Native vertex work scales
approximately as Nq² times the cube of the branch count; construction still
holds native strengths for one parent q, of order Nq times that cube. Force
constants and phonons also remain in memory. Uncompressed event arrays require
about 112 bytes per event, plus metadata; tens of millions of events require
several GB of disk. Repeated block reads can become I/O-bound. Use local scratch
storage when possible and benchmark a smaller grid before a large submission.

The 3D builder requires a stable harmonic spectrum with exactly three near-zero
Gamma translations (absolute frequency <= 1e-4 THz), which it fixes at zero.
Other imaginary or near-zero modes are rejected, not clipped. Input force
constants are not silently symmetrized. The fixed-translation boundary uses
eta=0 and zero vertex strength on these legs; this acoustic-limit treatment
needs separate validation. Flat resonant tetrahedra are rejected rather than
omitted. NAC-bearing inputs and generalized/shifted grids are currently
unsupported; Born-charge data must not be discarded just to bypass the check.

One-dimensional ribbons and chains (experimental)
------------------------------------------------------

A mesh such as ``[51,1,1]`` selects a one-dimensional collision integral along
the first reciprocal coordinate; the other two axes must be physically isolated
by vacuum. The mesh selects this route, independently of the geometry factors
used to normalize conductivity plots. All atomic displacement components and
all phonon branches remain present. This remains a spatially homogeneous model:
it does not propagate a heat pulse along the ribbon or impose end reservoirs.

Reuse existing, matching force constants with the following illustrative input.
Replace the source directory and choose branch indices from the actual ribbon
spectrum; branch 4 is not a universal optical excitation. The mesh and duration
are starting values, not converged settings.

.. code-block:: yaml

   tdbte:
     force-constants: calculations/graphene-ribbon/force-constants
     mesh: [51, 1, 1]
     temperature: 300
     branches: [4]
     excitation: 0.01
     duration-ps: 20
     max-step-ps: 0.5
     samples: 201
   output:
     result-dir: calculations/graphene-ribbon/tdbte-1d

Validate with ``nepkappa info input.yaml --for tdbte`` and execute with
``nepkappa tdbte input.yaml`` from the directory relative to which the paths
are defined. A six-section input can also generate FC2/FC3 first using Tersoff
or another supported calculator. Use supercells along the periodic direction
and check that interactions do not cross the vacuum. This implementation does
not infer vacuum sufficiency from the q mesh.

The 1D builder replaces tetrahedra with intervals. It integrates simple roots
of the piecewise-linear energy mismatch, using the inverse absolute slope and
half weight at shared endpoints. Vertex strengths and entropy variables use
the same two endpoint weights, with adjoint deposition for energy conservation
and detailed balance. Native strengths contain ``1/Nq``; there is no additional
mesh or cross-sectional-area factor in the occupation equation.

Gamma may contain three rigid translations or those translations plus an axial
rotation. The near-zero eigenvector subspace must match the corresponding
mass-weighted rigid motions to a projection tolerance of ``1e-5``. The frequency
threshold remains ``1e-4 THz``. Atom positions must form a contiguous ribbon
across the vacuum for the rotation test. Extra soft modes, imaginary modes and
non-rigid Gamma null spaces are rejected. Inputs are not silently repaired.

An identically zero *linear vertex interpolant* contributes zero and is counted
in ``zero_strength_intervals`` in the manifest and audit. This includes first
acoustic cells whose endpoints have fixed rigid-mode legs; it is an unresolved
low-frequency approximation, not proof of physically absent scattering inside
the interval. Flat resonant intervals with nonzero interpolated coupling raise
an error; no automatic broadening is applied. Near-tangent roots, branch
crossings, rotational sum rules and the acoustic boundary require independent
convergence and physical validation. Absolute graphene-nanoribbon relaxation
rates have not been benchmarked; ``validated_for_physical_dynamics`` stays false.

Plots default to Gamma--X along the active axis. Custom paths must stay on that
axis. Acoustic/optical temperature panels use the Gamma classification stored
in the kernel, which can include four acoustic branches. Indices elsewhere
remain frequency ordered and do not track polarization through crossings.

Reuse and failure handling
------------------------------

A completed build publishes ``tdbte-kernel/manifest.json`` with the mesh,
branch IDs, source hashes, backend version, normalization, boundary assumptions
and checksums of all mode/event files. Missing, duplicated, altered or incomplete
chunks are rejected before propagation. Treat kernel files as immutable while
a run is active.

To reuse a completed kernel for a different excitation or temperature, replace
``force-constants`` and ``mesh`` with::

   kernel: calculations/SiC/time-bte/tdbte-kernel/manifest.json

Use a fresh ``output.result-dir``. The mesh belongs to the reused kernel;
specifying ``mesh`` alongside ``kernel`` is an error. If the automatic route
finds its own completed cache, it checks source hashes, mesh and builder/backend
versions before reuse. An incomplete build is not automatically resumed;
``build-failure.json`` records the error when available. Preserve that directory
for diagnosis and select a new output directory after correcting the cause.
Existing ``result-dir/tdbte/`` results are never overwritten.

The previous single ``kernel.npz`` plus adjacent ``kernel.json`` interface also
works. Legacy research chunk manifests are not the versioned packaged format:
rebuild them from their recorded force constants instead of relabeling files.
For a single-file kernel, required numeric arrays are:

* ``frequency_THz``: (M,), nonnegative frequencies;
* ``zero_modes``: boolean (M,), marking only fixed zero-energy translations;
* ``parent``: integer (N,);
* ``daughter1``, ``daughter2``: integer (N,V);
* ``weights``: nonnegative (N,V), each row sums to one; V=2 for 1D intervals
  or V=4 for 3D tetrahedra;
* ``coefficients``: nonnegative (N,), with all integration and physical factors
  included so that dn/dt is in ps^-1.

The sidecar must supply ``model``, ``coefficient_time_unit: "ps"``,
``mode_branch_indices`` (M integers), and ``uniform_mode_weight`` (1/Nq).
Unequally weighted reduced grids are unsupported. The reader verifies discrete
shell consistency; it cannot verify the provenance or normalization claims of
an arbitrary external kernel.

New chunk manifests use ``nepkappa-tdbte-shell-v2``. Existing v1 manifests remain
readable; automatic caches from the previous builder must be rebuilt in a new
directory because builder identity is part of the cache check.

Outputs and interpretation
------------------------------

``result-dir/tdbte/`` contains:

* ``trajectories.npz``: occupations, frequencies, saved times, branch IDs, entropy
  and an independently integrated unexcited control;
* ``branch-energy.csv``: energy above the initial Bose state, in eV per primitive
  cell, including the uniform 1/Nq weight;
* ``audit.json``: kernel identity, source hashes, backend, integration settings,
  injected energy, energy drift, equilibrium drift and sampled entropy checks.

Occupations are never silently clipped. A failed audit raises an error and
marks diagnostic outputs as failed. ``nepkappa report`` includes the audit.
Passing these checks does not prove mesh convergence or absolute-rate accuracy;
even an incorrect global prefactor could conserve energy and produce entropy.
Kernel parity tests, independent rate/operator benchmarks, acoustic-boundary
checks and mesh/surface-quadrature convergence remain scientifically necessary.

Frequency/q-window excitation, prescribed-energy injection, time-dependent
pumps, checkpoint restart and dedicated Slurm submission are not implemented.
An outer batch job can run the same CLI; see the
:ref:`Slurm submission example <transport-slurm>`.
