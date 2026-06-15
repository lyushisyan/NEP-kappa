"""Configuration parsing utilities for NEP-kappa."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

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


def initialise_parser() -> argparse.ArgumentParser:
    """Build the internal workflow argument parser used by YAML input files."""
    parser = argparse.ArgumentParser(
        prog="nepkappa run",
        description="NEP + HiPhive + phono3py workflow",
    )

    parser.add_argument("--poscar", default="POSCAR", help="Input structure file")
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
        "--dim", type=int, nargs=3, default=[4, 4, 1], help="Supercell dimension"
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
    parser.add_argument(
        "--fc_calculator",
        choices=["traditional", "symfc", "alm"],
        default="traditional",
        help="[Finite displacement] Force constants solver",
    )
    parser.add_argument(
        "--fc_calculator_options",
        default=None,
        help="[Finite displacement] Options passed to the force constants solver",
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
        "structure": {"poscar"},
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
            "use_hiphive",
            "n_structures",
            "rattle_std",
            "cutoffs",
            "min_dist",
            "fc_calculator",
            "fc_calculator_options",
            "workdir",
            "vasp_kwargs",
        },
        "kappa": {"mesh", "temps", "method", "wigner"},
        "output": {"progress", "result_dir"},
    }


def yaml_arg_order():
    """Return a stable option order for parsed YAML values."""
    return [
        "poscar",
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
        "use_hiphive",
        "n_structures",
        "rattle_std",
        "cutoffs",
        "min_dist",
        "fc_calculator",
        "fc_calculator_options",
        "mesh",
        "temps",
        "method",
        "wigner",
        "progress",
        "result_dir",
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
    calculator_error = validate_calculator(args)
    if calculator_error:
        parser.error(calculator_error)
    temps_error = validate_temps(args.temps)
    if temps_error:
        parser.error(temps_error)
    return args


def iter_display_args(args):
    """Iterate over user-facing config values, hiding inactive route settings."""
    hiphive_only = {"n_structures", "rattle_std", "cutoffs", "min_dist"}
    finite_disp_only = {
        "fc_calculator",
        "fc_calculator_options",
    }
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
    for arg, value in vars(args).items():
        if value is None:
            continue
        if args.calculator == "vasp" and arg == "nep_model" and value is None:
            continue
        if not args.use_hiphive and arg in hiphive_only:
            continue
        if args.use_hiphive and arg in finite_disp_only:
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


def validate_calculator(args):
    """Validate calculator-specific settings."""
    if args.calculator == "nep" and not args.nep_model:
        return "--nep_model is required when --calculator nep"
    return None
