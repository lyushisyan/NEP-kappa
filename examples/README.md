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
No standalone SSCHA input is included.
The QHA-only example uses `force-constant.qha.run-transport: true` to request
QHA-volume conductivity. The QHA+SSCHA example sets
`force-constant.sscha.run-transport: true`, which also uses the QHA volume.
In the QHA+SSCHA example, every `kappa.temps` point must appear in `sscha.temps`,
and both must lie in the QHA range. SSCHA may sample more temperatures than
transport; QHA equilibrium volumes are interpolated between QHA grid points.

TD-BTE starts from matching `phono3py_disp.yaml`, `fc2.hdf5`, and `fc3.hdf5`.
`nepkappa run examples/tdbte.yaml` generates them and then builds its own
three-phonon kernel. A preceding conductivity run is not required. To reuse
existing files without new force calculations, use a compact input containing
only `tdbte.force-constants` (or `tdbte.kernel`) and `output`. After TD-BTE
finishes, `nepkappa plot examples/tdbte.yaml` draws branch and total excess
energy against time.

Server-specific Slurm resources can stay in an external batch script. Add the
optional top-level `parallel.force-constant` or `parallel.kappa` section only
when NEP-kappa should create stage-specific Slurm jobs. For SiC, quantitative
work also needs convergence checks and, where relevant, validated Born charges
and dielectric data for non-analytic corrections. Generated FourPhonon CONTROL
files disable those corrections; supply a validated custom CONTROL if needed.

Structure files under `structures/` support inputs and documentation. The
material-specific NEP is in `../potentials/3C-SiC/`.

The RTA phono3py examples set `plot.tau: nu`, so the scattering-rate figure
shows N and U together when the completed HDF5 contains both channels.
Wigner examples plot particle, coherence, and total conductivity when those
components are present. FourPhonon results generate separate 3ph and 4ph
scattering-rate figures and compare the conductivity solutions saved by the
run. A 3ph-only comparison curve requires a matching phono3py result in the
same result directory; it is not reconstructed from 3ph+4ph.
