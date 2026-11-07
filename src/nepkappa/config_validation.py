"""Scientific and command-specific input checks; no calculation is launched."""

import difflib
import math
import re
from pathlib import Path

from nepkappa.config_schema import normalize_yaml_key


def validate_temps(temps):
    """Validate the shared conductivity temperature convention."""
    if not isinstance(temps, (list, tuple)):
        return "kappa.temps must be a list of temperatures."
    count = len(temps)
    if count not in (1, 3):
        return (
            "kappa.temps expects [T] or [Tmin, Tmax, step]; "
            f"got {count} values: {temps}"
        )
    if any(not math.isfinite(value) for value in temps):
        return "kappa.temps values must be finite."
    if count == 1:
        if temps[0] <= 0:
            return "kappa.temps requires T > 0 K for conductivity."
        return None
    minimum, maximum, step = temps
    if minimum <= 0 or maximum < minimum or step <= 0:
        return "kappa.temps requires 0 < Tmin <= Tmax and step > 0."
    intervals = (maximum - minimum) / step
    if not math.isclose(intervals, round(intervals), rel_tol=1.0e-10, abs_tol=1.0e-10):
        return "kappa.temps step must divide Tmax - Tmin exactly."
    return None



def expand_transport_temperatures(temps):
    """Expand a validated single temperature or inclusive regular grid."""
    if len(temps) == 1:
        return [float(temps[0])]
    minimum, maximum, step = (float(value) for value in temps)
    intervals = round((maximum - minimum) / step)
    return [minimum + index * step for index in range(intervals + 1)]



def validate_temperature_relationships(args):
    """Check the QHA, SSCHA, and conductivity temperature sets together."""
    if not (args.qha_enabled or args.scph_enabled):
        return None
    qha_minimum = qha_maximum = None
    if args.qha_enabled:
        qha_minimum, qha_maximum, _ = args.qha_temps
    if args.scph_enabled:
        sscha_minimum, sscha_maximum, sscha_step = args.scph_temps
        if args.qha_enabled and (
            sscha_minimum < qha_minimum or sscha_maximum > qha_maximum
        ):
            return "force-constant.sscha.temps must lie within force-constant.qha.temps."
    if not args.kappa_enabled:
        return None
    requested = expand_transport_temperatures(args.temps)
    if args.qha_enabled and (
        requested[0] < qha_minimum or requested[-1] > qha_maximum
    ):
        return "kappa.temps must lie within force-constant.qha.temps."
    if args.scph_enabled:
        sscha_count = round((sscha_maximum - sscha_minimum) / sscha_step)
        sscha_points = [
            sscha_minimum + index * sscha_step for index in range(sscha_count + 1)
        ]
        if any(
            not any(math.isclose(value, point, rel_tol=1.0e-10, abs_tol=1.0e-8)
                    for point in sscha_points)
            for value in requested
        ):
            return "Every kappa.temps point must be a force-constant.sscha.temps point."
    return None



def validate_qha(args):
    """Validate isotropic quasi-harmonic approximation settings."""
    if not getattr(args, "qha_enabled", False):
        return None
    ratios = args.qha_volume_ratios
    if not isinstance(ratios, (list, tuple)) or len(ratios) < 5:
        return "force-constant.qha.volume-ratios requires at least 5 values."
    if any(not math.isfinite(value) or value <= 0 for value in ratios):
        return "force-constant.qha.volume-ratios values must be finite and positive."
    if any(right <= left for left, right in zip(ratios, ratios[1:])):
        return "force-constant.qha.volume-ratios must be strictly ascending and unique."
    for name, values in (("dim-fc2", args.qha_dim_fc2), ("mesh", args.qha_mesh)):
        if values is not None and (
            len(values) != 3 or any(int(value) <= 0 for value in values)
        ):
            return f"force-constant.qha.{name} must contain three positive integers."
    if not isinstance(args.qha_temps, (list, tuple)) or len(args.qha_temps) != 3:
        return "force-constant.qha.temps must contain exactly [minimum, maximum, step]."
    tmin, tmax, tstep = args.qha_temps
    if not all(math.isfinite(value) for value in (tmin, tmax, tstep)):
        return "force-constant.qha.temps values must be finite."
    if tmin < 0 or tmax <= tmin or tstep <= 0:
        return "force-constant.qha.temps must satisfy 0 <= minimum < maximum and step > 0."
    intervals = (tmax - tmin) / tstep
    if not math.isclose(intervals, round(intervals), rel_tol=1.0e-10, abs_tol=1.0e-10):
        return "force-constant.qha.temps step must divide the requested temperature interval."
    if not math.isfinite(args.qha_pressure):
        return "force-constant.qha.pressure must be finite."
    if not math.isfinite(args.qha_imaginary_frequency_tolerance):
        return "force-constant.qha.imaginary-frequency-tolerance must be finite."
    if not math.isfinite(args.qha_cutoff_frequency) or args.qha_cutoff_frequency < 0:
        return "force-constant.qha.cutoff-frequency must be finite and non-negative."
    if args.qha_relax_fmax <= 0:
        return "force-constant.qha.relax-fmax must be positive."
    if args.qha_relax_steps <= 0:
        return "force-constant.qha.relax-steps must be positive."
    if args.qha_displacement_distance <= 0:
        return "force-constant.qha.displacement-distance must be positive."
    if not isinstance(args.qha_vasp_relax_kwargs, dict):
        return "force-constant.qha.vasp-relax-kwargs must be a mapping."
    if not isinstance(args.qha_vasp_static_kwargs, dict):
        return "force-constant.qha.vasp-static-kwargs must be a mapping."
    if args.dimensionality != 3:
        return "QHA currently supports only structure.dimensionality: 3."
    return None



def validate_bubble(args):
    """Validate on-shell bubble postprocessing, independently of ASE sampling."""
    if not args.scph_enabled:
        return "bubble requires an scph section describing the completed SSCHA results."
    if args.dimensionality != 3:
        return "bubble currently supports only structure.dimensionality: 3."
    workdir = Path(args.scph_workdir)
    if workdir.is_absolute() or ".." in workdir.parts or not workdir.parts:
        return "force-constant.sscha.workdir must be a non-empty relative path inside output.result-dir."
    if any(value <= 0 for value in args.scph_bubble_mesh):
        return "force-constant.sscha.bubble-mesh must contain three positive integers."
    epsilons = args.scph_bubble_epsilons
    if not epsilons or any(not math.isfinite(v) or v <= 0 for v in epsilons):
        return "force-constant.sscha.bubble-epsilons must contain finite positive widths in THz."
    if len(set(epsilons)) != len(epsilons):
        return "force-constant.sscha.bubble-epsilons must not contain duplicates."
    points = args.scph_bubble_grid_points
    if points is not None and (not points or min(points) < 0 or len(set(points)) != len(points)):
        return "force-constant.sscha.bubble-grid-points must contain distinct non-negative BZ grid indices."
    if not math.isfinite(args.scph_cutoff_frequency) or args.scph_cutoff_frequency < 0:
        return "force-constant.sscha.cutoff-frequency must be finite and non-negative."
    return _validate_scph_temperatures(args)



def validate_scph(args):
    """Validate the Phonopy stochastic self-consistent phonon workflow."""
    if not getattr(args, "scph_enabled", False):
        return None
    if args.scph_run_transport and str(
        (args.lbte_parallel or {}).get("backend", "none")
    ).lower() == "slurm":
        return (
            "SSCHA and QHA-SSCHA require synchronous transport; nested Slurm RTA/LBTE "
            "is not supported, including submit: false. Set kappa.parallel.backend: "
            "none and run the whole workflow inside one allocation, or run direct "
            "Slurm RTA/LBTE separately on completed force constants."
        )
    if args.dimensionality != 3:
        return "SCPH currently supports only structure.dimensionality: 3."
    workdir = Path(args.scph_workdir)
    if workdir.is_absolute() or ".." in workdir.parts:
        return "force-constant.sscha.workdir must be a relative path inside output.result-dir."
    if args.calculator == "vasp":
        return (
            "scph requires an ASE calculator such as nep or mace; "
            "the command-driven VASP backend is not supported."
        )
    if args.scph_snapshots <= 0 or args.scph_iterations <= 0:
        return "force-constant.sscha.snapshots and force-constant.sscha.iterations must be positive."
    if not 0 <= args.scph_transient < args.scph_iterations:
        return "force-constant.sscha.transient must satisfy 0 <= transient < iterations."
    if len(args.scph_sscha_mesh) != 3 or any(
        value <= 0 for value in args.scph_sscha_mesh
    ):
        return "force-constant.sscha.sscha-mesh must contain three positive integers."
    if (
        not math.isfinite(args.scph_cutoff_frequency)
        or args.scph_cutoff_frequency < 0
    ):
        return "force-constant.sscha.cutoff-frequency must be finite and non-negative."
    if args.scph_run_transport:
        temps_error = validate_temps(args.temps)
        if temps_error:
            return temps_error
    return _validate_scph_temperatures(args)



def qha_sscha_transport_flags(args):
    """Return the three- and four-phonon switches for QHA-volume SSCHA."""
    configured_three_phonon = getattr(args, "qha_sscha_three_phonon", None)
    if configured_three_phonon is None:
        three_phonon = bool(getattr(args, "scph_run_transport", False))
    else:
        three_phonon = bool(configured_three_phonon)
    four_phonon = bool(getattr(args, "qha_sscha_four_phonon", False))
    return three_phonon, four_phonon



def validate_qha_sscha(args):
    """Validate coupled QHA+SSCHA transport selection."""
    if not getattr(args, "qha_enabled", False) or not getattr(
        args, "scph_enabled", False
    ):
        return "QHA-volume SSCHA requires force-constant.qha and force-constant.sscha."
    three_phonon, four_phonon = qha_sscha_transport_flags(args)
    if (
        getattr(args, "qha_sscha_three_phonon", None) is False
        and getattr(args, "scph_run_transport", False)
    ):
        return "QHA-volume SSCHA transport switches disagree."
    if three_phonon:
        temps_error = validate_temps(args.temps)
        if temps_error:
            return temps_error
    if four_phonon:
        if not args.fp_enabled:
            return (
                "QHA-volume SSCHA four-phonon transport requires force-constant.four-phonon."
            )
        if args.fp_harmonic_format != "shengbte":
            return (
                "QHA+SSCHA four-phonon transport requires "
                "force-constant.four-phonon.harmonic-format: shengbte."
            )
        backend = str((args.fp_parallel or {}).get("backend", "none")).lower()
        if backend == "slurm":
            return (
                "qha-sscha does not support nested FourPhonon Slurm jobs; submit "
                "the whole run as one batch job with force-constant.four-phonon.parallel.backend: none."
            )
    return None



def _validate_scph_temperatures(args):
    """Validate the shared SCPH/SSCHA temperature-range convention."""
    if len(args.scph_temps) != 3:
        return "force-constant.sscha.temps must contain [minimum, maximum, step]."
    tmin, tmax, tstep = args.scph_temps
    if not all(math.isfinite(value) for value in (tmin, tmax, tstep)):
        return "force-constant.sscha.temps values must be finite."
    if tmin < 0 or tmax < tmin or tstep <= 0:
        return "force-constant.sscha.temps must satisfy 0 <= minimum <= maximum and step > 0."
    intervals = (tmax - tmin) / tstep
    if not math.isclose(intervals, round(intervals), rel_tol=1.0e-10, abs_tol=1.0e-10):
        return "force-constant.sscha.temps step must divide the requested interval."
    return None



def validate_fourphonon(args):
    """Validate FourPhonon 3ph+4ph transport settings."""
    if not getattr(args, "fp_enabled", False):
        return None
    if getattr(args, "fp_compare_3ph", False):
        if args.fp_control or args.fp_wigner or args.fp_only_harmonic:
            return (
                "force-constant.four-phonon.compare-3ph requires generated CONTROL "
                "and ordinary transport (no custom control, Wigner or only-harmonic)."
            )
    args.fp_mesh = list(args.fp_mesh or args.mesh)
    args.fp_temps = list(args.fp_temps or args.temps)
    args.fp_scell = list(args.fp_scell or args.dim_fc2)
    if args.dimensionality != 3:
        return "FourPhonon transport currently supports only dimensionality: 3."
    for name, values in (("mesh", args.fp_mesh), ("scell", args.fp_scell)):
        if len(values) != 3 or any(int(value) <= 0 for value in values):
            return f"force-constant.four-phonon.{name} must contain three positive integers."
    temps_error = validate_temps(args.fp_temps)
    if temps_error:
        return temps_error.replace("kappa.temps", "force-constant.four-phonon.temps")
    if any(not math.isfinite(value) or value <= 0 for value in args.fp_temps):
        return "force-constant.four-phonon.temps values must be finite and positive."
    if len(args.fp_temps) == 3:
        tmin, tmax, tstep = args.fp_temps
        if tmax < tmin or tstep <= 0:
            return "force-constant.four-phonon.temps must satisfy minimum <= maximum and step > 0."
    if not math.isfinite(args.fp_scalebroad) or args.fp_scalebroad <= 0:
        return "force-constant.four-phonon.scalebroad must be finite and positive."
    if args.fp_mpi_processes <= 0 or args.fp_omp_threads <= 0:
        return "fourphonon MPI processes and OpenMP threads must be positive."
    if not args.fp_command:
        return "force-constant.four-phonon.command must not be empty."
    if not isinstance(args.fp_mpi_launcher, list) or not all(
        isinstance(token, str) for token in args.fp_mpi_launcher
    ):
        return "force-constant.four-phonon.mpi-launcher must be a list of strings."
    if "{nproc}" not in args.fp_mpi_launcher:
        return "force-constant.four-phonon.mpi-launcher must contain a {nproc} token."
    workdir = Path(args.fp_workdir)
    if workdir.is_absolute() or ".." in workdir.parts:
        return "force-constant.four-phonon.workdir must be relative to output.result-dir."
    if args.fp_nonanalytic and not args.fp_control:
        return (
            "force-constant.four-phonon.nonanalytic requires force-constant.four-phonon.control with dielectric "
            "and Born-charge data."
        )
    if args.fp_harmonic_format == "espresso" and not args.fp_control:
        return "force-constant.four-phonon.harmonic-format: espresso requires force-constant.four-phonon.control."
    samples = {
        "sample-3ph": args.fp_sample_3ph,
        "sample-3ph-phase-space": args.fp_sample_3ph_phase_space,
        "sample-4ph": args.fp_sample_4ph,
        "sample-4ph-phase-space": args.fp_sample_4ph_phase_space,
    }
    for name, value in samples.items():
        if value != -1 and value <= 0:
            return f"force-constant.four-phonon.{name} must be -1 or a positive integer."
    if args.fp_solver != "rta" and args.fp_sample_3ph > 0:
        return "Three-phonon sampling requires force-constant.four-phonon.solver: rta."
    if args.fp_solver == "full-iterative" and args.fp_sample_4ph > 0:
        return "Four-phonon sampling is incompatible with full-iterative solving."
    if args.fp_wigner:
        if args.fp_solver != "rta":
            return "force-constant.four-phonon.wigner currently requires force-constant.four-phonon.solver: rta."
        if args.fp_only_harmonic:
            return "force-constant.four-phonon.wigner requires a transport calculation, not only-harmonic."
        if args.fp_harmonic_format != "shengbte":
            return "force-constant.four-phonon.wigner currently requires harmonic-format: shengbte."
        if any(value != -1 for value in samples.values()):
            return (
                "force-constant.four-phonon.wigner requires all sample-* settings to be -1; "
                "the Wigner_Park branch does not support sampling."
            )

    settings = args.fp_parallel or {}
    if not isinstance(settings, dict):
        return "force-constant.four-phonon.parallel must be a mapping."
    settings = {normalize_yaml_key(key): value for key, value in settings.items()}
    args.fp_parallel = settings
    valid = {
        "backend", "partition", "account", "time", "memory", "nodes",
        "ntasks", "cpus_per_task", "job_name", "preamble", "extra_sbatch",
        "submit",
    }
    unknown = sorted(set(settings) - valid)
    if unknown:
        return f"Unknown force-constant.four-phonon.parallel key '{unknown[0].replace('_', '-')}'."
    backend = str(settings.get("backend", "none")).lower()
    if backend not in {"none", "slurm"}:
        return "force-constant.four-phonon.parallel.backend must be 'none' or 'slurm'."
    for key in ("nodes", "ntasks", "cpus_per_task"):
        if settings.get(key) is not None:
            try:
                if int(settings[key]) <= 0:
                    raise ValueError
            except (TypeError, ValueError):
                return f"force-constant.four-phonon.parallel.{key} must be a positive integer."
    for key in ("preamble", "extra_sbatch"):
        value = settings.get(key, [])
        if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
            return f"force-constant.four-phonon.parallel.{key} must be a list of strings."
    if "submit" in settings and not isinstance(settings["submit"], bool):
        return "force-constant.four-phonon.parallel.submit must be true or false."
    return None



def validate_force_constant_backends(args):
    """Validate thirdorder/Fourthorder cutoff conventions."""
    cutoffs = [("Fourthorder cutoff-fc4", args.cutoff_fc4)]
    if args.fc3_backend == "thirdorder":
        cutoffs.insert(0, ("thirdorder cutoff-fc3", args.cutoff_fc3))
    for backend, value in cutoffs:
        if not math.isfinite(value):
            return f"{backend} must be finite."
        if value == 0:
            return f"{backend} must not be zero."
        if value < 0 and not float(value).is_integer():
            return (
                f"{backend} must be a negative integer neighbor shell or a "
                "positive cutoff distance in nm."
            )
    if args.fc3_backend == "phono3py" and args.cutoff_fc3 is not None:
        if not math.isfinite(args.cutoff_fc3) or args.cutoff_fc3 <= 0:
            return "force-constant.cutoff-fc3 must be a positive distance in Angstrom for phono3py."
    return None



def validate_force_parallel(args):
    """Validate optional Slurm settings for displaced-structure force jobs."""
    settings = getattr(args, "force_parallel", {}) or {}
    if not isinstance(settings, dict):
        return "force-constant.parallel must be a mapping."
    settings = {normalize_yaml_key(key): value for key, value in settings.items()}
    args.force_parallel = settings
    valid_keys = {
        "backend",
        "jobs",
        "max_concurrent",
        "partition",
        "account",
        "time",
        "memory",
        "nodes",
        "ntasks",
        "cpus_per_task",
        "collect_time",
        "collect_memory",
        "collect_cpus_per_task",
        "job_name",
        "preamble",
        "extra_sbatch",
        "submit",
    }
    unknown = sorted(set(settings) - valid_keys)
    if unknown:
        key = unknown[0]
        matches = difflib.get_close_matches(key, sorted(valid_keys), n=1, cutoff=0.6)
        suggestion = f" Did you mean '{matches[0].replace('_', '-')}'?" if matches else ""
        return (
            f"Unknown force-constant.parallel key '{key.replace('_', '-')}'."
            f"{suggestion}"
        )

    backend = str(settings.get("backend", "none")).lower()
    if backend not in {"none", "slurm"}:
        return "force-constant.parallel.backend must be 'none' or 'slurm'."
    if backend == "none":
        return None
    if args.fc3_backend == "thirdorder":
        return (
            "force-constant.parallel currently supports the native phono3py "
            "and HiPhive routes, not fc3-backend: thirdorder."
        )
    try:
        values = [int(settings.get("jobs", 64))]
        for key in (
            "max_concurrent",
            "nodes",
            "ntasks",
            "cpus_per_task",
            "collect_cpus_per_task",
        ):
            if settings.get(key) is not None:
                values.append(int(settings[key]))
    except (TypeError, ValueError):
        return "Force Slurm job and resource counts must be integers."
    if min(values) < 1:
        return "Force Slurm job and resource counts must be positive."
    for key in ("preamble", "extra_sbatch"):
        value = settings.get(key, [])
        if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
            return f"force-constant.parallel.{key} must be a list of strings."
    if "submit" in settings and not isinstance(settings["submit"], bool):
        return "force-constant.parallel.submit must be true or false."
    return None



def validate_effective_geometry(args):
    """Validate low-dimensional effective-volume settings."""
    if args.effective_thickness is not None and args.effective_thickness <= 0:
        return "--effective_thickness must be positive."
    if args.effective_area is not None and args.effective_area <= 0:
        return "--effective_area must be positive."
    if args.dimensionality == 2 and args.effective_thickness is None:
        return (
            "structure.effective_thickness is required when "
            "structure.dimensionality is 2."
        )
    if args.dimensionality == 1 and args.effective_area is None:
        return (
            "structure.effective_area is required when "
            "structure.dimensionality is 1."
        )
    return None



def validate_calculator(args):
    """Validate calculator-specific settings."""
    args.calculator = str(args.calculator).strip().lower()
    if not args.calculator:
        return "calculator.name must not be empty."
    if args.calculator == "nep" and not args.nep_model:
        return "calculator.model is required when calculator.name is 'nep' (legacy alias: calculator.nep-model)."
    if args.nep_batch_size < 1:
        return "calculator.nep-batch-size must be a positive integer."
    if not str(args.nep_command).strip():
        return "calculator.nep-command must name a nep executable."
    if args.calculator != "nep" and (
        args.nep_force_backend != "cpu" or args.nep_batch_size != 128
        or args.nep_command != "nep"
    ):
        return "calculator.nep-* options require calculator.name: nep."
    if not isinstance(args.calculator_kwargs, dict):
        return "calculator.kwargs must be a mapping."
    if not isinstance(args.calculator_model_files, list) or not all(
        isinstance(path, str) for path in args.calculator_model_files
    ):
        return "calculator.model-files must be a list of paths."
    if args.calculator == "ase" and not args.calculator_factory:
        return "calculator.factory is required when calculator.name is 'ase'."
    if args.calculator in {"nep", "vasp", "mace"} and args.calculator_factory:
        return "calculator.factory is only valid for ASE or plugin calculators."
    if args.calculator == "mace" and not (
        args.mace_model or args.mace_foundation
    ):
        return (
            "calculator.model is required for a local MACE checkpoint, or set "
            "calculator.foundation (for example 'mp')."
        )
    if args.calculator == "mace":
        args.mace_device = args.mace_device or "cpu"
        args.mace_dtype = args.mace_dtype or "float64"
    if args.calculator == "mace" and args.mace_foundation:
        foundation = str(args.mace_foundation).replace("-", "_")
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", foundation):
            return "calculator.foundation must be a simple name such as 'mp'."
    if args.calculator != "mace" and any(
        value is not None
        for value in (
            args.mace_model,
            args.mace_foundation,
            args.mace_device,
            args.mace_dtype,
        )
    ):
        return (
            "calculator.foundation, device, and dtype are native MACE "
            "options and require calculator.name: mace."
        )
    if args.calculator_factory and not any(
        separator in str(args.calculator_factory) for separator in (":", ".")
    ):
        return (
            "calculator.factory must use 'module:callable' or "
            "'module.callable' syntax."
        )
    return None



def validate_lbte_parallel(args):
    """Validate optional Slurm settings for distributed RTA/LBTE calculations."""
    settings = getattr(args, "lbte_parallel", {}) or {}
    if not isinstance(settings, dict):
        return "kappa.parallel must be a mapping."
    settings = {normalize_yaml_key(key): value for key, value in settings.items()}
    args.lbte_parallel = settings

    valid_keys = {
        "backend",
        "jobs",
        "max_concurrent",
        "partition",
        "account",
        "time",
        "memory",
        "cpus_per_task",
        "collect_time",
        "collect_memory",
        "collect_cpus_per_task",
        "job_name",
        "preamble",
        "extra_sbatch",
        "submit",
        # Retained only to produce the existing targeted validation message.
        "nodes",
        "ntasks",
    }
    unknown = sorted(set(settings) - valid_keys)
    if unknown:
        key = unknown[0]
        matches = difflib.get_close_matches(key, sorted(valid_keys), n=1, cutoff=0.6)
        suggestion = f" Did you mean '{matches[0].replace('_', '-')}'?" if matches else ""
        return f"Unknown kappa.parallel key '{key.replace('_', '-')}'.{suggestion}"

    backend = str(settings.get("backend", "none")).lower()
    if backend not in {"none", "slurm"}:
        return "kappa.parallel.backend must be 'none' or 'slurm'."
    if backend == "none":
        return None
    if getattr(args, "qha_volumes", False):
        return (
            "QHA-volume transport requires synchronous transport; nested Slurm "
            "RTA is not supported, including submit: false. Use "
            "parallel.kappa.backend: none inside the temperature loop."
        )

    try:
        values = [int(settings.get("jobs", 32))]
        for key in ("cpus_per_task", "collect_cpus_per_task", "max_concurrent"):
            if settings.get(key) is not None:
                values.append(int(settings[key]))
        nodes = int(settings.get("nodes", 1))
        ntasks = int(settings.get("ntasks", 1))
    except (TypeError, ValueError):
        return "Slurm job counts and CPU settings must be integers."
    if min(values) < 1:
        return "Slurm job counts and CPU settings must be positive."
    if nodes != 1 or ntasks != 1:
        return (
            "kappa.parallel.nodes and ntasks must be 1; use jobs to distribute "
            "RTA/LBTE grid points across the Slurm array."
        )

    for key in ("preamble", "extra_sbatch"):
        value = settings.get(key, [])
        if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
            return f"kappa.parallel.{key} must be a list of strings."
    if "submit" in settings and not isinstance(settings["submit"], bool):
        return "kappa.parallel.submit must be true or false."
    return None



def validate_workflow_args(args, command, parser):
    """Validate and normalize the selected scientific operation in place."""
    targets = set(args.workflow_steps) if command == "run" else {command}
    if args.kappa_force_constants is not None and (
        command in {None, "run", "qha-kappa", "qha-sscha", "scph", "kappa4"}
    ):
        parser.error(
            "kappa.force-constants is for direct phono3py kappa reuse. "
            "Omit it when generating force constants with run or coupled transport."
        )
    if args.kappa_engine == "fourphonon" and "kappa" in targets:
        parser.error("kappa.engine: fourphonon cannot run the phono3py kappa stage.")
    if args.kappa_engine == "phono3py" and "kappa4" in targets:
        parser.error("kappa.engine: phono3py cannot run the FourPhonon kappa4 stage.")
    if command == "run" and args.qha_volumes and "qha-kappa" not in targets:
        parser.error("force-constant.qha.run-transport requires QHA-volume transport.")
    calculator_commands = {
        "relax",
        "fc2",
        "fc2fc3",
        "fc4",
        "qha",
        "qha-kappa",
        "scph",
        "qha-sscha",
    }
    if command is None or targets & calculator_commands:
        calculator_error = validate_calculator(args)
        if calculator_error:
            parser.error(calculator_error)
        if args.nep_force_backend == "gpu-batch" and (
            args.qha_enabled or args.scph_enabled or args.qha_sscha_enabled
            or args.fp_enabled or args.fc3_backend == "thirdorder" or "fc4" in targets
        ):
            parser.error(
                "calculator.nep-force-backend: gpu-batch supports native phono3py "
                "FC2/FC3 and HiPhive force batches only; QHA, SSCHA, Thirdorder "
                "and Fourthorder are not supported. Use cpu for those workflows."
            )
    if args.kappa_enabled and (command is None or command == "run" or targets & {
        "kappa", "qha-kappa", "qha-sscha", "scph", "kappa4"
    }):
        temps_error = validate_temps(args.temps)
        if temps_error:
            parser.error(temps_error)
    if not math.isfinite(args.plot_scattering_rate_min) or args.plot_scattering_rate_min < 0:
        parser.error("plot.scattering-rate-min must be finite and nonnegative (ps^-1).")
    if args.plot_tdbte_times is not None and (
        any(not math.isfinite(t) or t < 0 for t in args.plot_tdbte_times)
        or any(b <= a for a, b in zip(args.plot_tdbte_times, args.plot_tdbte_times[1:]))
    ):
        parser.error("plot.tdbte-times requires three finite, nonnegative, increasing times in ps.")
    geometry_error = validate_effective_geometry(args)
    if geometry_error:
        parser.error(geometry_error)
    if command is None or targets & {"kappa", "qha-kappa"} or args.scph_run_transport:
        parallel_error = validate_lbte_parallel(args)
        if parallel_error:
            parser.error(parallel_error)
    if command is None or targets & {"fc2", "fc2fc3", "fc4", "qha-sscha", "qha-kappa"}:
        force_constant_error = validate_force_constant_backends(args)
        if force_constant_error:
            parser.error(force_constant_error)
    if command is None or targets & {"fc2", "fc2fc3", "qha-sscha", "qha-kappa"}:
        force_parallel_error = validate_force_parallel(args)
        if force_parallel_error:
            parser.error(force_parallel_error)
    if "qha-sscha" in targets and str(
        (args.force_parallel or {}).get("backend", "none")
    ).lower() == "slurm":
        parser.error(
            "qha-sscha does not support nested force-constant Slurm arrays; "
            "submit the whole run as one batch job with parallel.backend: none."
        )
    if "qha-kappa" in targets:
        if not args.qha_enabled or not args.qha_volumes:
            parser.error("QHA-volume transport requires force-constant.qha.run-transport: true.")
        if args.method != "rta" or args.kappa_engine == "fourphonon":
            parser.error("qha-kappa currently supports phono3py RTA only.")
        if str((args.force_parallel or {}).get("backend", "none")).lower() == "slurm":
            parser.error("qha-kappa runs its force calculations inside one allocation; use parallel.force-constant.backend: none.")
    if "qha-sscha" in targets:
        qha_sscha_error = validate_qha_sscha(args)
        if qha_sscha_error:
            parser.error(qha_sscha_error)
    if command is None or targets & {"qha", "qha-sscha", "qha-kappa"}:
        qha_error = validate_qha(args)
        if qha_error:
            parser.error(qha_error)
    if command is None or targets & {"scph", "qha-sscha"}:
        scph_error = validate_scph(args)
        if scph_error:
            parser.error(scph_error)
    if command is None or command == "run" or targets & {"qha-kappa", "qha-sscha"}:
        relation_error = validate_temperature_relationships(args)
        if relation_error:
            parser.error(relation_error)
    if "bubble" in targets or args.scph_bubble:
        bubble_error = validate_bubble(args)
        if bubble_error:
            parser.error(bubble_error)
    if "tdbte" in targets:
        from nepkappa.tdbte import validate_options
        try:
            validate_options(args)
        except ValueError as exc:
            parser.error(str(exc))
    if "kappa4" in targets and not args.fp_enabled:
        parser.error("The selected command or workflow requires a fourphonon section.")
    if command is None or "kappa4" in targets or (
        "qha-sscha" in targets and qha_sscha_transport_flags(args)[1]
    ):
        fourphonon_error = validate_fourphonon(args)
        if fourphonon_error:
            parser.error(fourphonon_error)
