"""Configuration parsing utilities for NEP-kappa."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

try:
    import yaml
except ImportError:  # pragma: no cover - exercised only when PyYAML is absent.
    yaml = None


def str2bool(value):
    """Parse common boolean spellings accepted by the CLI."""
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in ("true", "1", "yes", "y", "t"):
        return True
    if normalized in ("false", "0", "no", "n", "f"):
        return False
    raise argparse.ArgumentTypeError(f"Invalid boolean value: '{value}'")


def json_dict(value):
    """Parse a JSON object used for calculator-specific keyword arguments."""
    if isinstance(value, dict):
        return value
    try:
        data = json.loads(value)
    except json.JSONDecodeError as exc:
        raise argparse.ArgumentTypeError(
            f"Expected a JSON object, got: {value}"
        ) from exc
    if not isinstance(data, dict):
        raise argparse.ArgumentTypeError("Expected a JSON object.")
    return data


def json_mapping(value):
    """Parse a JSON mapping, preserving nested stage dictionaries."""
    return json_dict(value)


def json_list(value):
    """Parse a JSON list used for nested YAML list options."""
    if isinstance(value, list):
        return value
    try:
        data = json.loads(value)
    except json.JSONDecodeError as exc:
        raise argparse.ArgumentTypeError(
            f"Expected a JSON list, got: {value}"
        ) from exc
    if not isinstance(data, list):
        raise argparse.ArgumentTypeError("Expected a JSON list.")
    return data


def initialise_parser() -> argparse.ArgumentParser:
    """Build the internal workflow argument parser used by YAML input files."""
    parser = argparse.ArgumentParser(
        prog="nepkappa run",
        description="NEP + HiPhive + phono3py workflow",
    )

    parser.add_argument("--poscar", default="POSCAR", help="Input structure file")
    parser.add_argument(
        "--dimensionality",
        type=int,
        choices=[1, 2, 3],
        default=3,
        help="Material dimensionality used for plot post-processing",
    )
    parser.add_argument(
        "--effective_thickness",
        "--effective-thickness",
        dest="effective_thickness",
        type=float,
        default=None,
        help="[2D] Effective film thickness in Angstrom",
    )
    parser.add_argument(
        "--effective_area",
        "--effective-area",
        dest="effective_area",
        type=float,
        default=None,
        help="[1D] Effective nanowire cross-sectional area in Angstrom^2",
    )
    parser.add_argument(
        "--vacuum_axis",
        "--vacuum-axis",
        dest="vacuum_axis",
        choices=["x", "y", "z"],
        default="z",
        help="[2D] Non-periodic film axis used for effective-thickness correction",
    )
    parser.add_argument(
        "--periodic_axis",
        "--periodic-axis",
        dest="periodic_axis",
        choices=["x", "y", "z"],
        default="z",
        help="[1D] Nanowire periodic axis used for effective-area correction",
    )
    parser.add_argument("--nep_model", default=None, help="Path to NEP model file")
    parser.add_argument(
        "--calculator",
        choices=["nep", "vasp"],
        default="nep",
        help="Force calculator used for force-constant generation",
    )
    parser.add_argument(
        "--vasp_command",
        default=None,
        help="Command used to run VASP, e.g. 'mpirun -np 64 vasp_std'",
    )
    parser.add_argument(
        "--vasp_path",
        default=None,
        help="Path to VASP executable, e.g. /root/software/vasp.6.4.3/bin/vasp_std",
    )
    parser.add_argument(
        "--potcar_path",
        default=None,
        help="Path to a POTCAR file or potential library directory",
    )
    parser.add_argument(
        "--vasp_workdir",
        default="vasp-runs",
        help="Directory under result_dir for VASP force calculations",
    )
    parser.add_argument(
        "--vasp_kwargs",
        type=json_dict,
        default={},
        help="JSON object with INCAR/KPOINTS options for VASP",
    )
    parser.add_argument(
        "--vasp_relax_workdir",
        default="vasp-relax",
        help="Directory under result_dir for VASP structure relaxation",
    )
    parser.add_argument(
        "--vasp_relax_stages",
        type=json_mapping,
        default=None,
        help="JSON object with coarse/fine VASP relaxation INCAR overrides",
    )
    parser.add_argument(
        "--do_relax",
        type=str2bool,
        nargs="?",
        const=True,
        default=False,
        help="Relax structure",
    )
    parser.add_argument(
        "--dim",
        type=int,
        nargs=3,
        default=None,
        help="[Deprecated] Supercell dimension used for both FC2 and FC3",
    )
    parser.add_argument(
        "--dim_fc2",
        "--dim-fc2",
        dest="dim_fc2",
        type=int,
        nargs=3,
        default=None,
        help="Phonon supercell dimension for FC2",
    )
    parser.add_argument(
        "--dim_fc3",
        "--dim-fc3",
        dest="dim_fc3",
        type=int,
        nargs=3,
        default=None,
        help="Supercell dimension for FC3",
    )

    parser.add_argument(
        "--use_hiphive",
        type=str2bool,
        nargs="?",
        const=True,
        default=False,
        help="Use HiPhive",
    )
    parser.add_argument(
        "--n_structures",
        type=int,
        default=50,
        help="[HiPhive] Number of structures",
    )
    parser.add_argument(
        "--rattle_std", type=float, default=0.03, help="[HiPhive] Rattle std"
    )
    parser.add_argument(
        "--cutoffs", type=float, nargs="+", default=[5.0], help="[HiPhive] Cutoffs"
    )
    parser.add_argument(
        "--min_dist",
        type=float,
        default=2.0,
        help="[HiPhive] Minimum atomic distance for rattling",
    )
    parser.add_argument("--mesh", type=int, nargs=3, default=[21, 21, 1])
    parser.add_argument(
        "--temps",
        type=float,
        nargs="+",
        default=[100, 1000, 100],
        help="One value or tmin tmax tstep",
    )
    parser.add_argument("--method", choices=["lbte", "rta"], default="lbte")
    parser.add_argument(
        "--isotope",
        type=str2bool,
        nargs="?",
        const=True,
        default=False,
        help="[Kappa] Include isotope scattering",
    )
    parser.add_argument(
        "--bfmp",
        type=float,
        default=1.0e6,
        help="[Kappa] Boundary mean free path in micrometer",
    )
    parser.add_argument(
        "--wigner",
        type=str2bool,
        nargs="?",
        const=True,
        default=False,
        help="Use Wigner transport via phono3py-wte (--tt wte)",
    )
    parser.add_argument(
        "--progress",
        type=str2bool,
        nargs="?",
        const=True,
        default=True,
        help="Show progress bars and timing summaries",
    )
    parser.add_argument(
        "--result_dir", default="result", help="Directory for outputs and run.log"
    )
    parser.add_argument(
        "--plot_layout",
        choices=["separate", "combined", "both"],
        default="separate",
        help="[Plot] Write separate figures, a combined figure, or both",
    )
    parser.add_argument(
        "--plot_path",
        choices=["seekpath", "custom"],
        default="seekpath",
        help="[Plot] High-symmetry path source for dispersion",
    )
    parser.add_argument(
        "--plot_path_points",
        type=json_dict,
        default=None,
        help="[Plot] Custom high-symmetry q-points",
    )
    parser.add_argument(
        "--plot_path_segments",
        type=json_list,
        default=None,
        help="[Plot] Custom high-symmetry path segments",
    )
    parser.add_argument(
        "--plot_tau",
        choices=["total", "normal", "umklapp", "all"],
        default="total",
        help="[Plot] Relaxation-time channel",
    )
    parser.add_argument(
        "--plot_kappa",
        choices=["x", "y", "z", "all"],
        default="all",
        help="[Plot] Thermal-conductivity tensor component",
    )
    parser.add_argument(
        "--plot_temperature",
        type=float,
        default=300.0,
        help="[Plot] Temperature used for relaxation-time plots",
    )
    parser.add_argument("--plot_dpi", type=int, default=300, help="[Plot] Figure DPI")

    return parser


def parse_input_file(filename):
    """Convert a YAML input file to workflow CLI tokens."""
    path = Path(filename)
    if not path.exists():
        return []
    if path.suffix.lower() not in (".yaml", ".yml"):
        raise ValueError("NEP-kappa input files must be YAML files ending in .yaml or .yml.")
    return parse_yaml_input_file(filename)


def parse_yaml_input_file(filename):
    """Convert structured YAML input into workflow CLI tokens."""
    if yaml is None:
        raise RuntimeError(
            "YAML input files require PyYAML. Install it with: pip install PyYAML"
        )
    with open(filename, "r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError("YAML input must be a mapping of sections and options.")

    flat = {}
    sections = yaml_input_sections()
    for section in sections:
        section_data = data.get(section, {})
        if section_data is None:
            continue
        if not isinstance(section_data, dict):
            raise ValueError(f"YAML section '{section}' must be a mapping.")
        for key, value in section_data.items():
            normalized_key = normalize_yaml_key(key)
            if section == "relaxation":
                if normalized_key == "enabled":
                    flat["do_relax"] = value
                elif normalized_key == "workdir":
                    flat["vasp_relax_workdir"] = value
                elif normalized_key == "stages":
                    flat["vasp_relax_stages"] = value
                else:
                    flat[normalized_key] = value
                continue
            if section == "force-constant":
                if normalized_key == "workdir":
                    flat["vasp_workdir"] = value
                elif normalized_key == "vasp_kwargs":
                    flat["vasp_kwargs"] = value
                else:
                    flat[normalized_key] = value
                continue
            if section == "plot":
                flat[f"plot_{normalized_key}"] = value
                continue
            flat[normalized_key] = value

    for key, value in data.items():
        normalized = normalize_yaml_key(key)
        if key not in sections and normalized not in sections:
            flat[normalized] = value

    aliases = {
        "name": "calculator",
    }
    for alias, canonical in aliases.items():
        if alias in flat and canonical not in flat:
            flat[canonical] = flat[alias]

    args = []
    for key in yaml_arg_order():
        if key in flat:
            append_arg(args, key, flat[key])
    return args


def yaml_input_sections():
    """Return supported YAML sections and their documented keys."""
    return {
        "structure": {
            "poscar",
            "dimensionality",
            "effective_thickness",
            "effective_area",
            "vacuum_axis",
            "periodic_axis",
        },
        "calculator": {
            "calculator",
            "name",
            "nep_model",
            "vasp_command",
            "vasp_path",
            "potcar_path",
        },
        "relaxation": {
            "enabled",
            "workdir",
            "stages",
        },
        "force-constant": {
            "dim",
            "dim_fc2",
            "dim_fc3",
            "use_hiphive",
            "n_structures",
            "rattle_std",
            "cutoffs",
            "min_dist",
            "workdir",
            "vasp_kwargs",
        },
        "kappa": {"mesh", "temps", "method", "isotope", "bfmp", "wigner"},
        "plot": {
            "layout",
            "path",
            "path_points",
            "path_segments",
            "tau",
            "kappa",
            "temperature",
            "dpi",
        },
        "output": {"progress", "result_dir"},
    }


def yaml_arg_order():
    """Return a stable option order for parsed YAML values."""
    return [
        "poscar",
        "dimensionality",
        "effective_thickness",
        "effective_area",
        "vacuum_axis",
        "periodic_axis",
        "nep_model",
        "calculator",
        "vasp_command",
        "vasp_path",
        "potcar_path",
        "vasp_workdir",
        "vasp_kwargs",
        "vasp_relax_workdir",
        "vasp_relax_stages",
        "do_relax",
        "dim",
        "dim_fc2",
        "dim_fc3",
        "use_hiphive",
        "n_structures",
        "rattle_std",
        "cutoffs",
        "min_dist",
        "mesh",
        "temps",
        "method",
        "isotope",
        "bfmp",
        "wigner",
        "progress",
        "result_dir",
        "plot_layout",
        "plot_path",
        "plot_path_points",
        "plot_path_segments",
        "plot_tau",
        "plot_kappa",
        "plot_temperature",
        "plot_dpi",
    ]


def normalize_yaml_key(key):
    """Normalize YAML keys to argparse option names."""
    return str(key).strip().replace("-", "_")


def append_arg(args, key, value):
    """Append one parsed YAML option to an argparse token list."""
    if value is None:
        return
    args.append(f"--{key}")
    if isinstance(value, bool):
        args.append(str(value).lower())
    elif isinstance(value, dict):
        args.append(json.dumps(value))
    elif key in {"plot_path_segments"}:
        args.append(json.dumps(value))
    elif isinstance(value, (list, tuple)):
        args.extend(str(item) for item in value)
    else:
        args.append(str(value))


def parse_workflow_args(config_path):
    """Parse workflow configuration from a YAML file."""
    parser = initialise_parser()
    if not Path(config_path).exists():
        parser.error(f"input file not found: {config_path}")
    tokens = parse_input_file(config_path)
    args = parser.parse_args(tokens)
    resolve_force_constant_dimensions(args)
    calculator_error = validate_calculator(args)
    if calculator_error:
        parser.error(calculator_error)
    temps_error = validate_temps(args.temps)
    if temps_error:
        parser.error(temps_error)
    geometry_error = validate_effective_geometry(args)
    if geometry_error:
        parser.error(geometry_error)
    return args


def parse_compare_args(config_path):
    """Parse a DFT-vs-NEP comparison YAML file."""
    if yaml is None:
        raise RuntimeError(
            "YAML input files require PyYAML. Install it with: pip install PyYAML"
        )
    path = Path(config_path)
    if not path.exists():
        raise FileNotFoundError(f"input file not found: {config_path}")
    if path.suffix.lower() not in (".yaml", ".yml"):
        raise ValueError("NEP-kappa compare files must end in .yaml or .yml.")

    with open(path, "r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError("Compare YAML input must be a mapping.")

    reference = require_mapping(data, "reference")
    candidate = require_mapping(data, "candidate")
    compare = require_mapping(data, "compare")
    plot = data.get("plot", {}) or {}
    structure = data.get("structure", {}) or {}
    if not isinstance(plot, dict):
        raise ValueError("YAML section 'plot' must be a mapping.")
    if not isinstance(structure, dict):
        raise ValueError("YAML section 'structure' must be a mapping.")
    reference = normalize_mapping_keys(reference)
    candidate = normalize_mapping_keys(candidate)
    compare = normalize_mapping_keys(compare)
    plot = normalize_mapping_keys(plot)
    structure = normalize_mapping_keys(structure)

    args = SimpleNamespace(
        dft_dir=required_section_value(reference, "reference", "dft_dir"),
        nep_dir=required_section_value(candidate, "candidate", "nep_dir"),
        compare_dir=required_section_value(compare, "compare", "compare_dir"),
        reference_label=str(reference.get("label", "DFT")),
        candidate_label=str(candidate.get("label", "NEP")),
        dimensionality=int(structure.get("dimensionality", 3)),
        effective_thickness=structure.get("effective_thickness", None),
        effective_area=structure.get("effective_area", None),
        vacuum_axis=str(structure.get("vacuum_axis", "z")),
        periodic_axis=str(structure.get("periodic_axis", "z")),
        plot_layout=str(plot.get("layout", "separate")),
        plot_path=str(plot.get("path", "seekpath")),
        plot_path_points=plot.get("path_points", None),
        plot_path_segments=plot.get("path_segments", None),
        plot_tau=str(plot.get("tau", "total")),
        plot_kappa=str(plot.get("kappa", "all")),
        plot_temperature=float(plot.get("temperature", 300.0)),
        plot_dpi=int(plot.get("dpi", 300)),
        mesh=None,
    )
    validate_compare_args(args)
    return args


def require_mapping(data, section):
    """Return a required YAML section mapping."""
    section_data = data.get(section)
    if not isinstance(section_data, dict):
        raise ValueError(f"YAML section '{section}' is required and must be a mapping.")
    return section_data


def normalize_mapping_keys(data):
    """Return a copy of a mapping with YAML keys normalized to option names."""
    return {normalize_yaml_key(key): value for key, value in data.items()}


def required_section_value(section_data, section, key):
    """Return a required section value."""
    value = section_data.get(key)
    if value in (None, ""):
        raise ValueError(f"YAML section '{section}' requires '{key}'.")
    return str(value)


def validate_compare_args(args):
    """Validate compare command options."""
    if args.dimensionality not in (1, 2, 3):
        raise ValueError("structure.dimensionality must be 1, 2, or 3.")
    geometry_error = validate_effective_geometry(args)
    if geometry_error:
        raise ValueError(geometry_error)
    if args.plot_layout not in ("separate", "combined", "both"):
        raise ValueError("plot.layout must be separate, combined, or both.")
    if args.plot_path not in ("seekpath", "custom"):
        raise ValueError("plot.path must be seekpath or custom.")
    if args.plot_tau not in ("total", "normal", "umklapp", "all"):
        raise ValueError("plot.tau must be total, normal, umklapp, or all.")
    if args.plot_kappa not in ("x", "y", "z", "all"):
        raise ValueError("plot.kappa must be x, y, z, or all.")
    if args.plot_dpi <= 0:
        raise ValueError("plot.dpi must be positive.")


def format_compare_config(args):
    """Return a readable multi-line comparison configuration summary."""
    lines = ["Running comparison with configuration:"]
    for key in [
        "dft_dir",
        "nep_dir",
        "compare_dir",
        "reference_label",
        "candidate_label",
        "dimensionality",
        "effective_thickness",
        "effective_area",
        "vacuum_axis",
        "periodic_axis",
        "plot_layout",
        "plot_path",
        "plot_tau",
        "plot_kappa",
        "plot_temperature",
        "plot_dpi",
    ]:
        value = getattr(args, key)
        if value is None:
            continue
        if args.dimensionality != 2 and key in {"effective_thickness", "vacuum_axis"}:
            continue
        if args.dimensionality != 1 and key in {"effective_area", "periodic_axis"}:
            continue
        lines.append(f"  {key:<18} : {value}")
    return "\n".join(lines)


def resolve_force_constant_dimensions(args):
    """Resolve deprecated dim into explicit FC2 and FC3 supercell dimensions."""
    legacy_dim = args.dim
    default_dim = [4, 4, 1]
    if args.dim_fc2 is None:
        args.dim_fc2 = list(legacy_dim or default_dim)
    if args.dim_fc3 is None:
        args.dim_fc3 = list(legacy_dim or default_dim)


def iter_display_args(args):
    """Iterate over user-facing config values, hiding inactive route settings."""
    compatibility_only = {"dim"}
    hiphive_only = {"n_structures", "rattle_std", "cutoffs", "min_dist"}
    vasp_only = {
        "vasp_command",
        "vasp_path",
        "potcar_path",
        "vasp_workdir",
        "vasp_kwargs",
    }
    vasp_relax_only = {
        "vasp_relax_workdir",
        "vasp_relax_stages",
    }
    film_only = {"effective_thickness", "vacuum_axis"}
    wire_only = {"effective_area", "periodic_axis"}
    for arg, value in vars(args).items():
        if value is None:
            continue
        if arg in compatibility_only:
            continue
        if args.dimensionality != 2 and arg in film_only:
            continue
        if args.dimensionality != 1 and arg in wire_only:
            continue
        if args.calculator == "vasp" and arg == "nep_model" and value is None:
            continue
        if not args.use_hiphive and arg in hiphive_only:
            continue
        if args.calculator != "vasp" and arg in vasp_only:
            continue
        if (args.calculator != "vasp" or not args.do_relax) and arg in vasp_relax_only:
            continue
        yield arg, value


def format_config(args):
    """Return a readable multi-line configuration summary."""
    lines = ["Running Workflow with configuration:"]
    for arg, value in iter_display_args(args):
        lines.append(f"  {arg:<15} : {value}")
    return "\n".join(lines)


def validate_temps(temps):
    """Validate phono3py temperature input shape."""
    if not isinstance(temps, (list, tuple)):
        return "--temps must be a list of floats."
    count = len(temps)
    if count not in (1, 3):
        return (
            f"--temps expects 1 value (single temperature) or 3 values "
            f"(tmin tmax tstep), got {count}: {temps}"
        )
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
    if args.calculator == "nep" and not args.nep_model:
        return "--nep_model is required when --calculator nep"
    return None
