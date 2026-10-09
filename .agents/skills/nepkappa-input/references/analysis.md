# Existing-result plots and time-dependent dynamics

Read this reference only for plotting or TD-BTE requests. The matching source
and documentation remain authoritative; do not invent YAML keys from a paper
figure's labels.

## Plotting: determine capability from artifacts

Inspect `src/nepkappa/plot.py` and the plotting section of
`docs/source/input_files.rst`. Inventory the result directory before drafting:

| Available data | Supported built-in figures |
| --- | --- |
| `fc2.hdf5` and matching phonon metadata | Dispersion, DOS, volume heat capacity, group velocity |
| Above plus compatible transport HDF5 | Add thermal conductivity and transport-derived properties |
| Transport data with `gamma` | Add relaxation time / scattering rate |
| Transport data with `mode_kappa` | Add cumulative conductivity |

Metadata is searched in order: `phono3py_disp.yaml`, `phonopy.yaml`,
`phonopy_disp.yaml`. Bare FC2 is insufficient. Confirm that cell, species,
FC2 supercell, primitive mapping, and NAC data belong together; do not borrow
metadata from an unrelated calculation. Missing artifacts are not permission
to regenerate force constants.

Use `plot.layout: separate`, `combined`, or `both`; four harmonic panels form
a 2-by-2 combined figure. An FC2-only request needs no potential or FC3.
Set `output.result-dir` to the existing result directory; built-in plots are
written beneath its `plots/`. Distinguish input validation from authorized
plot generation, which can replace existing plot files.

Check mesh and temperature selection. An explicitly requested mesh must match
the transport filename; multiple candidates without an explicit mesh are
ambiguous. The current reader rejects unsupported split-grid/Wigner-only
schemas and extra broadening axes. Do not rename a file or delete it to hide
an incompatibility. Harmonic plotting is a fallback when no transport file
exists, not a recovery path for a corrupt transport file.

Transport plotting chooses the nearest available temperature with a warning;
report the actual temperature. Harmonic heat capacity uses `kappa.temps`
(`[T]` or `[start, stop, step]`) and a mesh integration, not a new transport
calculation. Imaginary modes require a stability warning, not absolute values
or a claim of stable thermodynamics. A visually dense band path does not
establish q-mesh convergence.

The publication Figure1–6 scripts are separate research tools: do not promise that
`nepkappa plot` reproduces their custom panels, Wigner bubbles, pair-frequency
maps, or TD-BTE snapshots. Check the plotting API before promising arbitrary
panel-selection keys; `layout` controls arrangement, not a figure whitelist.

## Time-dependent BTE: scope and prerequisites

Read `docs/source/tdbte.rst`, `examples/tdbte.yaml`, and `src/nepkappa/tdbte.py`
before preparing this route. It evolves spatially homogeneous populations at
fixed frequencies using a three-phonon energy-shell operator. It is not a
laser-absorption, electron-phonon, coherent-phonon, or spatial transport model.

- Choose exactly one route. For automatic construction, set
  `tdbte.force-constants` to a directory containing matching
  `phono3py_disp.yaml`, `fc2.hdf5`, and `fc3.hdf5`, and set `tdbte.mesh` to three
  integers >= 2. The source is reused; no calculator or new force jobs are needed.
  For reuse, set `tdbte.kernel` to a numeric NPZ with its JSON sidecar, or a
  completed packaged `tdbte-kernel/manifest.json`. Omit `mesh` on this route.
  Do not add the removed `experimental` key or substitute an LBTE matrix.
- Automatic construction uses the full unshifted diagonal q grid. The branch
  count comes from the structure, not a hardcoded SiC convention. NAC-bearing
  inputs, unstable modes and extra near-zero modes are rejected. Do not discard
  Born charges or alter force constants just to bypass a check.
- The builder and solver stream checked chunks; no manual concatenation is
  required. Old Figure6 research manifests are not the packaged format: rebuild
  from their recorded sources instead of changing a filename. `kappa.mesh`
  does not control TD-BTE. Dense grids can still be expensive in time and disk;
  chunking is not a convergence claim or permission to launch a large job.
- Check sidecar provenance, ps time unit, zero-based branch IDs, and uniform
  full-grid mode weight. Unequally weighted irreducible meshes are unsupported.
  Never assume branches `[3, 4, 5]` describe every material's optical modes.
- `temperature` defines the initial Bose state. `excitation: 0.01` increases
  selected occupations by 1%, not by 1 K, 1% total energy, or a laser fluence.
  Explain which branches at which q points are excited (currently every q).
- Use `--for tdbte` for stage validation or a custom plan `[tdbte]` for `run`.
  Validation checks options, not source existence or physical normalization.
  Do not run the ODE just to validate an input.
- Use a fresh output directory: an existing `result-dir/tdbte/` is refused.
  A completed automatic cache is reused only with matching input hashes, mesh
  and backend/builder versions; partial builds are not resumed automatically.
  `report` includes the audit; standard `plot` does not render the trajectory.
- Separate energy conservation, equilibrium-control, and entropy checks from
  physical validation of collision prefactors and absolute relaxation rates.
  Passing numerical checks does not establish a physically validated result.

If a kernel is missing but matching FC2/FC3 and metadata are available, prepare
the automatic route. If both sources are absent, deliver a labeled draft and
prerequisite list, not a runnable claim. Do not substitute an RTA exponential decay.
