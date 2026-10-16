# Examples

This directory contains nine phonon/transport 3C-SiC inputs and one time-dependent
BTE input. Run commands from the repository root. The bundled supercells and
meshes demonstrate the workflow; they are not converged research settings.
Keep a separate `output.result-dir` for each calculation.

| Input | Calculation | Before running |
| --- | --- | --- |
| [vasp-rta-3ph.yaml](vasp-rta-3ph.yaml) | VASP forces, phono3py 3ph RTA | Set the VASP executable and licensed POTCAR path |
| [nep-rta-wigner-3ph.yaml](nep-rta-wigner-3ph.yaml) | NEP, phono3py 3ph RTA + SMM19 Wigner | Check the SiC model and q-mesh |
| [nep-lbte-wigner-3ph.yaml](nep-lbte-wigner-3ph.yaml) | NEP, phono3py 3ph LBTE + SMM19 Wigner | Plan LBTE memory and convergence |
| [nep-rta-3ph-4ph.yaml](nep-rta-3ph-4ph.yaml) | NEP, FourPhonon 3ph RTA + 4ph RTA | Install Fourthorder and FourPhonon |
| [nep-rta-wigner-3ph-4ph.yaml](nep-rta-wigner-3ph-4ph.yaml) | NEP, FourPhonon 3ph/4ph RTA + Wigner coherence | Set the Wigner_Park executable |
| [nep-lbte-3ph-rta-4ph.yaml](nep-lbte-3ph-rta-4ph.yaml) | NEP, 3ph LBTE + 4ph RTA | Install FourPhonon CPU |
| [nep-lbte-3ph-lbte-4ph.yaml](nep-lbte-3ph-lbte-4ph.yaml) | NEP, 3ph LBTE + 4ph LBTE | Plan the full iterative solve |
| [nep-qha-rta-3ph.yaml](nep-qha-rta-3ph.yaml) | NEP, QHA volume at each temperature + 3ph RTA | Check volume and temperature ranges |
| [nep-qha-sscha-rta-3ph.yaml](nep-qha-sscha-rta-3ph.yaml) | NEP, QHA-volume SSCHA + 3ph RTA | Converge SSCHA sampling |
| [tdbte.yaml](tdbte.yaml) | NEP FC2/FC3, then time-dependent BTE | Check the SiC model, supercells, and q-mesh |

For a first read-only check:

```bash
nepkappa info examples/nep-rta-wigner-3ph.yaml
```

To calculate, copy an input and set a distinct result directory. `nepkappa run`
executes the selected phonon/transport or TD-BTE route; `nepkappa plot` and `nepkappa report`
operate on completed results. `info` does not run simulations or verify
external executables. The six-section TD-BTE example generates its own matching
FC2/FC3 and phono3py metadata; see [the TD-BTE guide](../docs/source/tdbte.rst).

The phonon/transport inputs have six sections: `structure`, `calculator`,
`force-constant`, `kappa`, `plot`, and `output`. TD-BTE has the same six sections
with `tdbte` in place of `kappa`. To run only one calculation step from a
six-section file, use `nepkappa fc2` or `nepkappa fc2fc3`.
If `structure.relaxation: true`, run `nepkappa relax` before a standalone
force-constant command. The 4ph inputs set
`kappa.method-3ph` and `kappa.method-4ph` separately. FourPhonon supports
RTA/RTA, LBTE/RTA, and LBTE/LBTE; 3ph RTA with 4ph LBTE is unavailable.
Wigner_Park requires RTA/RTA. The phono3py-only Wigner examples use SMM19.
The chapter inputs include a standalone SSCHA example (see below).
The QHA transport example uses `force-constant.qha.run-transport: true` to request
QHA-volume conductivity. The QHA+SSCHA example sets
`force-constant.sscha.run-transport: true`, which also uses the QHA volume.
In the QHA+SSCHA example, every `kappa.temps` point must appear in `sscha.temps`,
and both must lie in the QHA range. SSCHA may sample more temperatures than
transport; QHA equilibrium volumes are interpolated between QHA grid points.

TD-BTE starts from matching `phono3py_disp.yaml`, `fc2.hdf5`, and `fc3.hdf5`.
`nepkappa run examples/tdbte.yaml` generates them and then builds its own
three-phonon kernel. A preceding conductivity run is not required. To reuse
existing files without new force calculations, use a compact input containing
`tdbte.force-constants` with `tdbte.mesh` (or `tdbte.kernel` without a mesh),
plus the excitation settings and `output`. After TD-BTE finishes,
`nepkappa plot examples/tdbte.yaml` generates six panels: three dispersions
colored by mode excess energy at 0, 10 and 20 ps, branch excess energies,
acoustic branch temperatures, and optical branch temperatures.

Server-specific Slurm resources can stay in an external batch script. Add the
optional top-level `parallel.force-constant` or `parallel.kappa` section for
stage-specific execution resources. FourPhonon MPI/OpenMP settings also use
`parallel.kappa` when running locally or inside an external batch allocation. For SiC, quantitative
work also needs convergence checks and, where relevant, validated Born charges
and dielectric data for non-analytic corrections. Generated FourPhonon CONTROL
files disable those corrections; supply a validated custom CONTROL if needed.

Structure files under `structures/` support inputs and documentation. The
material-specific NEP is in `../potentials/3C-SiC/`.

## Plotting and comparisons

The RTA phono3py examples set `plot.tau: nu`, so the scattering-rate figure
shows N and U together when both channels are saved. Wigner plots show
particle, coherence and total conductivity. FourPhonon overlays 3ph and 4ph
rates in one scattering panel and plots the available conductivity solutions.

The FourPhonon templates in this directory calculate the combined 3ph+4ph
case. For a matched three-phonon reference, set
`force-constant.four-phonon.compare-3ph: true` in an ordinary non-Wigner input
with generated CONTROL. It runs both cases sequentially with the same
FourPhonon executable. The chapter's `fourphonon.yaml` already enables this
comparison and uses Thirdorder FC3 with Fourthorder FC4.

## Inputs matching the user manual

The numbered documentation chapters provide downloadable inputs in
[`docs/source/inputs/`](../docs/source/inputs/). Save them in the working directory
specified by the chapter; their relative paths differ from the repository-root
paths used by this directory's templates.

| Chapter | Inputs | Intended stage |
| --- | --- | --- |
| 3. Quick start | [quickstart.yaml](../docs/source/inputs/quickstart.yaml) | `run`, then `plot` |
| 5. Harmonic phonons | [harmonic.yaml](../docs/source/inputs/harmonic.yaml) | `relax`, then `fc2` and `plot` |
| 5. QHA | [qha.yaml](../docs/source/inputs/qha.yaml) | `fc2`, then `plot` |
| 5. SSCHA | [sscha.yaml](../docs/source/inputs/sscha.yaml), [qha-sscha.yaml](../docs/source/inputs/qha-sscha.yaml) | `run`, then `plot` |
| 6. LBTE | [lbte.yaml](../docs/source/inputs/lbte.yaml) | `kappa` reusing FC2/FC3, then `plot` |
| 6. Four-phonon comparison | [fourphonon.yaml](../docs/source/inputs/fourphonon.yaml) | `run`, then `plot` |
| 6. Wigner with 3ph+4ph | [wigner-4ph.yaml](../docs/source/inputs/wigner-4ph.yaml) | `kappa` reusing FC2/FC3/FC4, then `plot` |
| 7. TD-BTE | [tdbte.yaml](../docs/source/inputs/tdbte.yaml) | `tdbte` reusing FC2/FC3, then `plot` |

The generic and chapter inputs serve different demonstrations:

- Generic FourPhonon templates use the default phono3py FC3 backend and a
  12×12×12 transport mesh at 300 K. The chapter comparison uses Thirdorder,
  an 8×8×8 mesh and 300, 600 and 900 K. Both routes are supported; their
  force-constant cutoffs and solver settings must be matched for comparisons.
- Generic QHA/SSCHA templates include thermal conductivity. Chapter 5 inputs
  disable transport to focus on phonons and temperature effects.
- Generic TD-BTE generates its own force constants on a 5×5×5 dynamics mesh.
  Chapter 7 reuses force constants on a 6×6×6 dynamics mesh. Both integrate
  to 20 ps and use the six-panel plot.
- The chapter quick start uses `plot.tau: total` and `plot.kappa: x`; the
  generic three-phonon RTA template uses `nu` and the diagonal average (`all`).

Use the chapter inputs when following the manual. These differences are
intentional example choices, not recommended convergence settings.
