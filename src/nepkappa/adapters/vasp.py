"""VASP file, execution, force, and relaxation adapter."""

from __future__ import annotations

from pathlib import Path
import shlex
import shutil

from ase.io import read, write


class VaspBackend:
    """Own all VASP-specific I/O and process orchestration.

    The adapter intentionally receives a process runner. This keeps Slurm and
    terminal logging outside VASP-specific code and makes the backend testable
    without launching VASP.
    """

    def __init__(
        self,
        config,
        *,
        output_dir,
        force_root,
        relax_root,
        relaxed_structure,
        run_command,
    ):
        self.config = config
        self.output_dir = Path(output_dir)
        self.force_root = Path(force_root)
        self.relax_root = Path(relax_root)
        self.relaxed_structure = Path(relaxed_structure)
        self.run_command = run_command

    def command(self):
        """Resolve the VASP command from YAML or common server locations."""
        if getattr(self.config, "vasp_command", None):
            return shlex.split(self.config.vasp_command)
        if getattr(self.config, "vasp_path", None):
            return [self.config.vasp_path]
        detected = self.detect_path()
        if detected is not None:
            print(f"  - Detected VASP executable: {detected}")
            return [str(detected)]
        raise FileNotFoundError(
            "VASP executable not found. Set calculator.vasp_path or "
            "calculator.vasp_command in the YAML file."
        )
    def detect_path(self):
        """Find a likely VASP executable on common local/server paths."""
        candidates = [
            Path("/root/software/vasp.6.4.3/bin/vasp_std"),
            Path("/root/software/vasp.6.4.3/build/std/vasp"),
        ]
        for candidate in candidates:
            try:
                if candidate.is_file():
                    return candidate
            except OSError:
                continue
        try:
            discovered = sorted(Path("/root/software").glob("vasp*/bin/vasp_std"))
        except OSError:
            discovered = []
        for candidate in discovered:
            try:
                if candidate.is_file():
                    return candidate
            except OSError:
                continue
        resolved = shutil.which("vasp_std")
        return Path(resolved) if resolved else None

    def resolve_potcar(self, atoms):
        """Resolve a POTCAR file or assemble one from a potential directory."""
        configured = getattr(self.config, "potcar_path", None)
        symbols = list(dict.fromkeys(atoms.get_chemical_symbols()))
        if configured:
            path = Path(configured)
            if path.is_file():
                return path, (path,)
            if path.is_dir():
                return self.assemble_potcar(path, symbols)
            raise FileNotFoundError(f"POTCAR path not found: {path}")

        detected = self.detect_potcar(symbols)
        if detected is not None:
            print(f"  - Detected POTCAR source: {detected}")
            return detected
        raise FileNotFoundError(
            "POTCAR not found. Set calculator.potcar_path in the YAML file."
        )

    def detect_potcar(self, symbols):
        """Find a likely POTCAR source under common local locations."""
        if len(symbols) == 1:
            symbol = symbols[0]
            candidates = [
                Path(f"/root/software/potpaw_PBE.64/{symbol}/POTCAR"),
                Path(f"/root/software/vasp.6.4.3/testsuite/POTCARS/POTCAR.{symbol}"),
                Path(f"/root/software/vasp.6.4.3/potpaw_PBE/{symbol}/POTCAR"),
                Path(f"/root/software/potpaw_PBE/{symbol}/POTCAR"),
            ]
            for candidate in candidates:
                try:
                    if candidate.is_file():
                        return candidate, (candidate,)
                except OSError:
                    continue
        for library in (
            Path("/root/software/potpaw_PBE.64"),
            Path("/root/software/vasp.6.4.3/potpaw_PBE"),
            Path("/root/software/potpaw_PBE"),
            Path("/root/software/vasp.6.4.3/testsuite/POTCARS"),
        ):
            try:
                if library.is_dir():
                    return self.assemble_potcar(library, symbols)
            except OSError:
                continue
        return None

    def assemble_potcar(self, library, symbols):
        """Return one POTCAR or create a correctly joined multi-species file."""
        library = Path(library)
        sources = []
        chunks = []
        for symbol in symbols:
            source = self.find_potcar(library, symbol)
            if source is None:
                raise FileNotFoundError(f"missing POTCAR for {symbol} under {library}")
            sources.append(source)
            chunks.append(source.read_bytes())
        if len(sources) == 1:
            return sources[0], tuple(sources)

        combined = self.output_dir / "POTCAR.combined"
        combined.parent.mkdir(parents=True, exist_ok=True)
        combined.write_bytes(
            b"".join(chunk.rstrip(b"\r\n") + b"\n" for chunk in chunks)
        )
        (self.output_dir / "POTCAR.spec").write_text(
            "\n".join(str(source) for source in sources) + "\n",
            encoding="utf-8",
        )
        return combined, tuple(sources)

    @staticmethod
    def find_potcar(library, symbol):
        candidates = [
            library / symbol / "POTCAR",
            library / f"{symbol}_sv" / "POTCAR",
            library / f"{symbol}_pv" / "POTCAR",
            library / f"POTCAR.{symbol}",
            library / "POTCAR",
        ]
        return next((path for path in candidates if path.is_file()), None)

    def write_run_inputs(self, run_dir, atoms, incar_overrides=None, system=None):
        """Write POSCAR, INCAR, KPOINTS, POTCAR, and POTCAR.spec."""
        run_dir = Path(run_dir)
        run_dir.mkdir(parents=True, exist_ok=True)
        write(str(run_dir / "POSCAR"), atoms, format="vasp", direct=True, vasp5=True)
        kwargs = self.combined_kwargs(incar_overrides, system=system)
        self.write_incar(run_dir / "INCAR", kwargs)
        self.write_kpoints(run_dir / "KPOINTS", kwargs)
        potcar, sources = self.resolve_potcar(atoms)
        shutil.copyfile(potcar, run_dir / "POTCAR")
        (run_dir / "POTCAR.spec").write_text(
            "\n".join(str(source) for source in sources) + "\n",
            encoding="utf-8",
        )

    def run_forces(self, atoms, label, index):
        run_dir = self.force_root / label.lower() / f"{index:05d}"
        return self.run_forces_in_dir(run_dir, atoms, label, index)

    def run_forces_in_dir(self, run_dir, atoms, label, index):
        self.write_run_inputs(run_dir, atoms, system="NEP-kappa VASP single point")
        print(f"    VASP {label.upper()} #{index}: {run_dir}")
        return_code = self.run_command(self.command(), cwd=run_dir)
        if return_code != 0:
            raise RuntimeError(f"VASP failed in {run_dir} with return code {return_code}")
        return self.read_forces(run_dir)

    def run_order_job(self, atoms, job_dir, label, *, consumer):
        """Run a VASP job whose full vasprun.xml is consumed externally."""
        job_dir = Path(job_dir)
        self.write_run_inputs(
            job_dir,
            atoms,
            system=f"NEP-kappa VASP {label} single point",
        )
        print(f"    VASP forces ({label}): {job_dir}")
        return_code = self.run_command(self.command(), cwd=job_dir)
        if return_code != 0:
            raise RuntimeError(
                f"VASP failed in {job_dir} with return code {return_code}"
            )
        output = job_dir / "vasprun.xml"
        if not output.is_file() or output.stat().st_size == 0:
            raise RuntimeError(f"{consumer} requires {output}; VASP did not write it.")
        return output

    def combined_kwargs(self, overrides=None, system=None):
        kwargs = dict(getattr(self.config, "vasp_kwargs", {}) or {})
        defaults = {
            "system": system or "PESMaker single point",
            "gga": "PE",
            "lreal": "Auto",
            "ibrion": -1,
            "nsw": 0,
            "algo": "Normal",
            "ediff": 1.0e-6,
            "sigma": 0.02,
            "ismear": 0,
            "prec": "Accurate",
            "nelm": 150,
            "lwave": False,
            "lcharg": False,
        }
        for key, value in defaults.items():
            kwargs.setdefault(key, value)
        if overrides:
            kwargs.update(overrides)
        return kwargs

    def write_incar(self, path, kwargs=None):
        kwargs = self.combined_kwargs() if kwargs is None else kwargs
        skip = {"kpts", "gamma", "txt", "directory", "command", "xc"}
        with Path(path).open("w", encoding="utf-8") as handle:
            for key, value in kwargs.items():
                if key.lower() not in skip:
                    handle.write(f"{key.upper()} = {self.format_value(value)}\n")

    def write_kpoints(self, path, kwargs=None):
        path = Path(path)
        kwargs = self.combined_kwargs() if kwargs is None else kwargs
        if "kspacing" in {str(key).lower() for key in kwargs}:
            if path.exists():
                path.unlink()
            return
        kpoints = kwargs.get("kpts", [1, 1, 1])
        if isinstance(kpoints, int):
            kpoints = [kpoints] * 3
        if len(kpoints) != 3:
            raise ValueError("vasp_kwargs.kpts must be an integer or a length-3 list")
        mode = "Gamma" if kwargs.get("gamma", True) else "Monkhorst-Pack"
        path.write_text(
            "Automatic mesh\n"
            "0\n"
            f"{mode}\n"
            f"{int(kpoints[0])} {int(kpoints[1])} {int(kpoints[2])}\n"
            "0 0 0\n",
            encoding="utf-8",
        )

    @classmethod
    def format_value(cls, value):
        if isinstance(value, bool):
            return ".TRUE." if value else ".FALSE."
        if isinstance(value, (list, tuple)):
            return " ".join(cls.format_value(item) for item in value)
        return str(value)

    @staticmethod
    def read_forces(run_dir):
        run_dir = Path(run_dir)
        last_error = None
        for filename in ("vasprun.xml", "OUTCAR"):
            path = run_dir / filename
            if not path.exists():
                continue
            try:
                return read(str(path), index=-1).get_forces()
            except Exception as exc:
                last_error = exc
        raise RuntimeError(
            f"Could not read VASP forces from {run_dir / 'vasprun.xml'} "
            f"or {run_dir / 'OUTCAR'}"
        ) from last_error

    @staticmethod
    def default_relax_stages():
        return {
            "coarse": {
                "system": "NEP-kappa VASP coarse relaxation",
                "nsw": 80,
                "ibrion": 2,
                "isif": 3,
                "ediff": 1.0e-5,
                "ediffg": -0.05,
                "prec": "Normal",
            },
            "fine": {
                "system": "NEP-kappa VASP fine relaxation",
                "nsw": 150,
                "ibrion": 2,
                "isif": 3,
                "ediff": 1.0e-6,
                "ediffg": -0.01,
                "prec": "Accurate",
            },
        }

    def relax_stages(self):
        stages = getattr(self.config, "vasp_relax_stages", None)
        if not stages:
            stages = self.default_relax_stages()
        if not isinstance(stages, dict):
            raise ValueError("vasp_relax.stages must be a mapping")
        ordered = []
        for name in ("coarse", "fine"):
            if name in stages:
                ordered.append((name, self._validated_relax_stage(name, stages[name])))
        for name, value in stages.items():
            if name not in {"coarse", "fine"}:
                ordered.append((name, self._validated_relax_stage(name, value)))
        return ordered

    @staticmethod
    def _validated_relax_stage(name, value):
        value = value or {}
        if not isinstance(value, dict):
            raise ValueError(f"vasp_relax.stages.{name} must be a mapping")
        return value

    @staticmethod
    def read_relaxed_structure(run_dir):
        run_dir = Path(run_dir)
        last_error = None
        for filename in ("CONTCAR", "vasprun.xml", "OUTCAR"):
            path = run_dir / filename
            if not path.exists() or path.stat().st_size == 0:
                continue
            try:
                return read(str(path), index=-1)
            except Exception as exc:
                last_error = exc
        raise RuntimeError(
            f"Could not read relaxed VASP structure from {run_dir}"
        ) from last_error

    def run_relax_stage(self, atoms, stage_name, stage_index, overrides):
        run_dir = self.relax_root / f"{stage_index:02d}-{stage_name}"
        system = overrides.get("system", f"NEP-kappa VASP {stage_name} relaxation")
        self.write_run_inputs(run_dir, atoms, incar_overrides=overrides, system=system)
        print(f"    VASP relax {stage_name}: {run_dir}")
        return_code = self.run_command(self.command(), cwd=run_dir)
        if return_code != 0:
            raise RuntimeError(
                f"VASP relaxation failed in {run_dir} with return code {return_code}"
            )
        relaxed = self.read_relaxed_structure(run_dir)
        write(str(run_dir / "POSCAR_next"), relaxed, format="vasp", direct=True, vasp5=True)
        return relaxed

    def relax(self, atoms):
        print("  - Relaxing structure using VASP")
        current = atoms.copy()
        for index, (name, overrides) in enumerate(self.relax_stages(), 1):
            print(f"  - VASP relaxation stage {index}: {name}")
            current = self.run_relax_stage(current, name, index, overrides)
        self.relaxed_structure.parent.mkdir(parents=True, exist_ok=True)
        write(
            str(self.relaxed_structure),
            current,
            format="vasp",
            direct=True,
            vasp5=True,
        )
        print(f"  - Relaxed structure saved to {self.relaxed_structure}")
        return current
