"""Accepted YAML keys and scalar conversions shared by parsing and planning."""

import argparse
import difflib
import json

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



def normalize_yaml_key(key):
    """Normalize YAML keys to argparse option names."""
    return str(key).strip().replace("-", "_")



def raise_unknown_yaml_key(key, location, valid_keys):
    """Raise a helpful error for an unsupported YAML option."""
    normalized_valid = sorted({str(item).replace("_", "-") for item in valid_keys})
    matches = difflib.get_close_matches(
        str(key).replace("_", "-"), normalized_valid, n=1, cutoff=0.6
    )
    suggestion = f" Did you mean '{matches[0]}'?" if matches else ""
    raise ValueError(f"Unknown YAML key '{key}' in {location}.{suggestion}")



def yaml_input_sections():
    """Return supported YAML sections and their documented keys."""
    return {
        "workflow": {"steps"},
        "structure": {
            "poscar",
            "relaxation",
            "dimensionality",
            "effective_thickness",
            "effective_area",
            "vacuum_axis",
            "periodic_axis",
        },
        "calculator": {
            "name",
            "factory",
            "kwargs",
            "model_files",
            "model",
            "foundation",
            "device",
            "dtype",
            "nep_model",
            "nep_force_backend",
            "nep_batch_size",
            "nep_command",
            "vasp_command",
            "vasp_path",
            "potcar_path",
        },
        "force-constant": {
            "dim_fc2",
            "dim_fc3",
            "fc3_backend",
            "cutoff_fc3",
            "thirdorder_command",
            "fc3_workdir",
            "dim_fc4",
            "cutoff_fc4",
            "fourthorder_command",
            "fc4_workdir",
            "use_hiphive",
            "compact_fc",
            "format",
            "n_structures",
            "rattle_std",
            "cutoffs",
            "min_dist",
            "workdir",
            "vasp_kwargs",
            "parallel",
            "qha",
            "sscha",
            "four_phonon",
        },
        "kappa": {
            "enabled",  # Internal field derived from force-constant switches.
            "force_constants",
            "engine",
            "qha_volumes",  # Internal QHA route flag.
            "mesh",
            "temps",
            "method",
            "method_3ph",
            "method_4ph",
            "command",
            "isotope",
            "bfmp",
            "wigner",
            "parallel",
        },
        "fourphonon": {
            "enabled",
            "command",
            "workdir",
            "control",
            "compare_3ph",
            "fc2",
            "fc3",
            "fc4",
            "harmonic_format",
            "mesh",
            "temps",
            "scell",
            "solver",
            "wigner",
            "scalebroad",
            "isotopes",
            "nonanalytic",
            "only_harmonic",
            "sample_3ph",
            "sample_3ph_phase_space",
            "sample_4ph",
            "sample_4ph_phase_space",
            "mpi_launcher",
            "mpi_processes",
            "omp_threads",
            "omp_stacksize",
            "parallel",
        },
        "plot": {
            "layout",
            "path",
            "path_points",
            "path_segments",
            "tau",
            "kappa",
            "temperature",
            "dpi",
            "scattering_rate_min",
            "tdbte_times",
        },
        "qha": {
            "equilibrium_phonons",
            "enabled",
            "run_transport",
            "volume_ratios",
            "dim_fc2",
            "mesh",
            "temps",
            "eos",
            "pressure",
            "relax_internal",
            "relax_fmax",
            "relax_steps",
            "displacement_distance",
            "imaginary_frequency_tolerance",
            "cutoff_frequency",
            "vasp_relax_kwargs",
            "vasp_static_kwargs",
        },
        "scph": {
            "enabled",
            "workdir",
            "temps",
            "initial_fc2",
            "born",
            "snapshots",
            "iterations",
            "transient",
            "sscha_mesh",
            "random_seed",
            "cutoff_frequency",
            "fc_calculator",
            "fc_calculator_options",
            "save_datasets",
            "run_transport",
            "transport_fc3",
            "transport_metadata",
            "bubble",
            "bubble_mesh",
            "bubble_epsilons",
            "bubble_grid_points",
        },
        "qha-sscha": {
            "three_phonon",
            "four_phonon",
        },
        "tdbte": {"kernel", "force_constants", "mesh", "temperature", "duration_ps", "max_step_ps", "excitation", "branches", "samples"},
        "output": {"progress", "result_dir"},
    }


def load_yaml_strict(stream):
    """Load YAML while rejecting duplicate mapping keys at every level."""
    if yaml is None:
        raise RuntimeError("YAML input files require PyYAML. Install it with: pip install PyYAML")

    class UniqueKeyLoader(yaml.SafeLoader):
        pass

    def construct_unique_mapping(loader, node, deep=False):
        mapping = {}
        for key_node, value_node in node.value:
            key = loader.construct_object(key_node, deep=deep)
            if key in mapping:
                raise ValueError(
                    f"Duplicate YAML key '{key}' at line {key_node.start_mark.line + 1}."
                )
            mapping[key] = loader.construct_object(value_node, deep=deep)
        return mapping

    UniqueKeyLoader.add_constructor(
        yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
        construct_unique_mapping,
    )
    return yaml.load(stream, Loader=UniqueKeyLoader)
