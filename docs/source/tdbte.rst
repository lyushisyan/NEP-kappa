Time-dependent BTE
==================

The ``tdbte`` stage evolves spatially homogeneous phonon occupations at fixed
frequencies using a nonlinear three-phonon energy-shell operator. It can build
the operator from existing force constants and solve it in disk-backed chunks.
No calculator, potential, or separately prepared collision matrix is needed.

This is a population-dynamics model, not a laser-absorption calculation. It
does not include electrons, coherent phonons, spatial drift, a heat bath,
four-phonon collisions, or frequency changes during evolution. Numerical audits
do not by themselves validate absolute physical relaxation rates.

From force constants to dynamics
--------------------------------

Put the matching ``phono3py_disp.yaml``, ``fc2.hdf5``, and ``fc3.hdf5`` in one
directory. The YAML supplies the structure, primitive mapping and supercells;
two HDF5 arrays without this metadata are insufficient. Use the same structure
and force-constant convention for all three files.

For a two-atom primitive cell such as 3C-SiC::

   tdbte:
     force-constants: calculations/SiC/fc
     mesh: [5, 5, 5]
     temperature: 300
     branches: [3, 4, 5]
     excitation: 0.01
     duration-ps: 10
   output:
     result-dir: calculations/SiC/time-bte

From the launch directory::

   nepkappa validate input.yaml --for tdbte
   nepkappa run input.yaml
   nepkappa report calculations/SiC/time-bte

``nepkappa run input.yaml`` automatically selects the dynamic TD-BTE route
when the file contains ``tdbte`` and ``output`` sections only.
``nepkappa tdbte input.yaml`` runs the same stage. Validation checks options without constructing a kernel
or starting integration; file contents are checked at execution. Relative
paths resolve from the launch directory, not from the YAML file.
``examples/tdbte.yaml`` is the complete template. Its mesh is a starting value,
not a convergence recommendation.

Input choices
--------------

* ``force-constants`` is the source directory, not a request to compute forces.
* ``mesh`` is the **q grid for the collision integral**, independent of
  ``kappa.mesh``. Supply three integers of at least 2. The builder uses the full
  unshifted, diagonal three-dimensional grid, not an irreducible grid.
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

There is no ``experimental`` switch. Do not include it in new inputs.

Kernel construction and memory
------------------------------

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
mode-by-mode collision matrix is loaded into memory. Optional acceleration::

   pip install 'nepkappa[tdbte]'

When Numba is installed, block reductions are compiled and compared against
the NumPy reference before integration. Otherwise the same solver uses NumPy.
The audit records the backend. This optional dependency does not change the
input format or physical model.

Chunking does not remove the cost of dense q meshes. Native vertex work scales
approximately as Nq² times the cube of the branch count; construction still
holds native strengths for one parent q, of order Nq times that cube. Force
constants and phonons also remain in memory. Uncompressed event arrays require
about 112 bytes per event, plus metadata; tens of millions of events require
several GB of disk. Repeated block reads can become I/O-bound. Use local scratch
storage when possible and benchmark a smaller grid before a large submission.

The builder requires a stable harmonic spectrum with exactly three near-zero
Gamma translations (absolute frequency <= 1e-4 THz), which it fixes at zero.
Other imaginary or near-zero modes are rejected, not clipped. Input force
constants are not silently symmetrized. The fixed-translation boundary uses
eta=0 and zero vertex strength on these legs; this acoustic-limit treatment
needs separate validation. Flat resonant tetrahedra are rejected rather than
omitted. NAC-bearing inputs and generalized/shifted grids are currently
unsupported; Born-charge data must not be discarded just to bypass the check.

Reuse and failure handling
--------------------------

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
* ``daughter1``, ``daughter2``: integer (N,4);
* ``weights``: nonnegative (N,4), each row sums to one;
* ``coefficients``: nonnegative (N,), with all integration and physical factors
  included so that dn/dt is in ps^-1.

The sidecar must supply ``model``, ``coefficient_time_unit: "ps"``,
``mode_branch_indices`` (M integers), and ``uniform_mode_weight`` (1/Nq).
Unequally weighted reduced grids are unsupported. The reader verifies discrete
shell consistency; it cannot verify the provenance or normalization claims of
an arbitrary external kernel.

Outputs and interpretation
--------------------------

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
An outer batch job can run the same CLI. ``nepkappa plot`` does not yet draw
TD-BTE trajectories; the manuscript's Figure6 panels and channel-resolved
postprocessing remain separate from the packaged plotting workflow.
