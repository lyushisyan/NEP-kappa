"""Configuration parsing utilities for NEP-kappa."""

from __future__ import annotations

import argparse
import difflib
import json
import math
import re
from pathlib import Path

from nepkappa.config_models import WorkflowConfig

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

    parser.add_argument(
        "--workflow_preset",
        choices=[
            "three-phonon",
            "four-phonon",
            "qha",
            "qha-kappa",
            "scph",
            "qha-sscha",
            "custom",
        ],
        default="three-phonon",
        help="High-level workflow selected by `nepkappa run`",
    )
    parser.add_argument(
        "--workflow_steps",
        nargs="+",
        default=None,
        help="Advanced custom stage list used with workflow.preset: custom",
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
        default="nep",
        help="Force calculator: nep, vasp, mace, ase, or an installed plugin name",
    )
    parser.add_argument(
        "--calculator_factory",
        "--calculator-factory",
        dest="calculator_factory",
        default=None,
        help="Import path for a generic ASE Calculator class or factory",
    )
    parser.add_argument(
        "--calculator_kwargs",
        "--calculator-kwargs",
        dest="calculator_kwargs",
        type=json_dict,
        default={},
        help="JSON object passed to the external ASE Calculator factory",
    )
    parser.add_argument(
        "--calculator_model_files",
        "--calculator-model-files",
        dest="calculator_model_files",
        nargs="+",
        default=[],
        help="Model files included explicitly in cache and provenance hashes",
    )
    parser.add_argument(
        "--mace_model",
        "--mace-model",
        dest="mace_model",
        default=None,
        help="[MACE] Local checkpoint path or foundation-model identifier",
    )
    parser.add_argument(
        "--mace_foundation",
        "--mace-foundation",
        dest="mace_foundation",
        default=None,
        help="[MACE] Foundation factory suffix, e.g. mp for mace_mp",
    )
    parser.add_argument(
        "--mace_device",
        "--mace-device",
        dest="mace_device",
        default=None,
        help="[MACE] Inference device, e.g. cpu, cuda, or mps",
    )
    parser.add_argument(
        "--mace_dtype",
        "--mace-dtype",
        dest="mace_dtype",
        choices=["float32", "float64"],
        default=None,
        help="[MACE] Default floating-point precision",
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
        "--fc3_backend",
        "--fc3-backend",
        dest="fc3_backend",
        choices=["phono3py", "thirdorder"],
        default="phono3py",
        help="Backend used to generate FC3",
    )
    parser.add_argument(
        "--cutoff_fc3",
        "--cutoff-fc3",
        dest="cutoff_fc3",
        type=float,
        default=None,
        help="FC3 cutoff: phono3py positive Angstrom; thirdorder negative shell or positive nm",
    )
    parser.add_argument(
        "--thirdorder_command",
        "--thirdorder-command",
        dest="thirdorder_command",
        default=None,
        help="[thirdorder] Command used to run thirdorder_vasp.py",
    )
    parser.add_argument(
        "--fc3_workdir",
        "--fc3-workdir",
        dest="fc3_workdir",
        default="fc3-thirdorder-runs",
        help="[thirdorder] Directory under result_dir for FC3 calculations",
    )
    parser.add_argument(
        "--dim_fc4",
        "--dim-fc4",
        dest="dim_fc4",
        type=int,
        nargs=3,
        default=None,
        help="[FourPhonon] Supercell dimension for FC4",
    )
    parser.add_argument(
        "--cutoff_fc4",
        "--cutoff-fc4",
        dest="cutoff_fc4",
        type=float,
        default=-2.0,
        help="[FourPhonon] FC4 cutoff passed to Fourthorder_vasp.py",
    )
    parser.add_argument(
        "--fourthorder_command",
        "--fourthorder-command",
        dest="fourthorder_command",
        default=None,
        help="[FourPhonon] Command used to run Fourthorder_vasp.py",
    )
    parser.add_argument(
        "--fc4_workdir",
        "--fc4-workdir",
        dest="fc4_workdir",
        default="fc4-runs",
        help="[FourPhonon] Directory under result_dir for FC4 calculations",
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
        "--compact_fc",
        "--compact-fc",
        dest="compact_fc",
        type=str2bool,
        nargs="?",
        const=True,
        default=True,
        help="[Force constants] Write compact FC2/FC3 arrays; false writes full arrays",
    )
    parser.add_argument(
        "--fc_format",
        "--fc-format",
        "--format",
        dest="fc_format",
        choices=["phono3py", "shengbte", "both"],
        default="phono3py",
        help="[Force constants] Export format for FC2/FC3",
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
        "--kappa_engine", choices=["phono3py", "fourphonon"], default=None,
        help="Transport engine selected by kappa.engine",
    )
    parser.add_argument(
        "--qha_volumes", type=str2bool, nargs="?", const=True, default=False,
        help="Recompute FC2/FC3 and 3ph conductivity at each QHA volume",
    )
    parser.add_argument(
        "--kappa_command",
        "--kappa-command",
        dest="kappa_command",
        default=None,
        help="[Kappa] Full custom phono3py command run in output.result_dir",
    )
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
        help="Use Wigner transport via phono3py SMM19 (--tt smm19)",
    )
    parser.add_argument(
        "--lbte_parallel",
        type=json_dict,
        default={},
        help="[LBTE] Slurm parallel calculation settings",
    )
    parser.add_argument(
        "--force_parallel",
        type=json_dict,
        default={},
        help="[FC] Slurm array settings for displaced-structure forces",
    )
    parser.add_argument("--fp_enabled", type=str2bool, nargs="?", const=True, default=False)
    parser.add_argument("--fp_command", default="ShengBTE_cpu")
    parser.add_argument("--fp_workdir", default="fourphonon")
    parser.add_argument("--fp_control", default=None)
    parser.add_argument("--fp_fc2", default=None)
    parser.add_argument("--fp_fc3", default=None)
    parser.add_argument("--fp_fc4", default=None)
    parser.add_argument(
        "--fp_harmonic_format",
        choices=["shengbte", "espresso"],
        default="shengbte",
    )
    parser.add_argument("--fp_mesh", type=int, nargs=3, default=None)
    parser.add_argument("--fp_temps", type=float, nargs="+", default=None)
    parser.add_argument("--fp_scell", type=int, nargs=3, default=None)
    parser.add_argument(
        "--fp_solver",
        choices=["rta", "3ph-iterative", "full-iterative"],
        default="rta",
    )
    parser.add_argument(
        "--fp_wigner", type=str2bool, nargs="?", const=True, default=False,
        help="Use FourPhonon Wigner_Park for 3+4ph population and coherence transport",
    )
    parser.add_argument("--fp_scalebroad", type=float, default=1.0)
    parser.add_argument("--fp_isotopes", type=str2bool, nargs="?", const=True, default=False)
    parser.add_argument("--fp_nonanalytic", type=str2bool, nargs="?", const=True, default=False)
    parser.add_argument("--fp_only_harmonic", type=str2bool, nargs="?", const=True, default=False)
    parser.add_argument("--fp_sample_3ph", type=int, default=-1)
    parser.add_argument("--fp_sample_3ph_phase_space", type=int, default=-1)
    parser.add_argument("--fp_sample_4ph", type=int, default=-1)
    parser.add_argument("--fp_sample_4ph_phase_space", type=int, default=-1)
    parser.add_argument(
        "--fp_mpi_launcher",
        type=json_list,
        default=["mpirun", "-np", "{nproc}"],
    )
    parser.add_argument("--fp_mpi_processes", type=int, default=1)
    parser.add_argument("--fp_omp_threads", type=int, default=1)
    parser.add_argument("--fp_omp_stacksize", default="1G")
    parser.add_argument("--fp_parallel", type=json_dict, default={})
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
        choices=["total", "normal", "umklapp", "nu", "all"],
        default="total",
        help="[Plot] Relaxation-time channel (nu plots N and U together)",
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
    parser.add_argument(
        "--qha_enabled",
        type=str2bool,
        nargs="?",
        const=True,
        default=False,
        help="Enable quasi-harmonic approximation settings",
    )
    parser.add_argument(
        "--qha_volume_ratios",
        type=float,
        nargs="+",
        default=[0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06],
        help="QHA ratios relative to the input cell volume, in ascending order",
    )
    parser.add_argument(
        "--qha_dim_fc2",
        type=int,
        nargs=3,
        default=None,
        help="QHA harmonic-force-constant supercell dimension",
    )
    parser.add_argument(
        "--qha_mesh",
        type=int,
        nargs=3,
        default=[31, 31, 31],
        help="QHA phonon mesh",
    )
    parser.add_argument(
        "--qha_temps",
        type=float,
        nargs=3,
        default=[0, 1100, 10],
        help="QHA temperature range: minimum maximum step in K",
    )
    parser.add_argument(
        "--qha_eos",
        choices=["vinet", "birch_murnaghan", "murnaghan"],
        default="vinet",
        help="QHA equation of state",
    )
    parser.add_argument(
        "--qha_pressure",
        type=float,
        default=0.0,
        help="QHA external pressure in GPa",
    )
    parser.add_argument(
        "--qha_relax_internal",
        type=str2bool,
        nargs="?",
        const=True,
        default=True,
        help="Relax internal coordinates at each fixed QHA cell",
    )
    parser.add_argument("--qha_relax_fmax", type=float, default=1.0e-3)
    parser.add_argument("--qha_relax_steps", type=int, default=200)
    parser.add_argument(
        "--qha_displacement_distance",
        type=float,
        default=0.01,
        help="QHA FC2 finite-displacement distance in Angstrom",
    )
    parser.add_argument(
        "--qha_imaginary_frequency_tolerance",
        type=float,
        default=-0.1,
        help="Lowest allowed QHA mesh frequency in THz",
    )
    parser.add_argument(
        "--qha_cutoff_frequency",
        type=float,
        default=1.0e-4,
        help="Exclude QHA thermal modes below this frequency in THz",
    )
    parser.add_argument(
        "--qha_vasp_relax_kwargs",
        type=json_mapping,
        default={},
        help="VASP INCAR overrides for fixed-cell QHA relaxation",
    )
    parser.add_argument(
        "--qha_vasp_static_kwargs",
        type=json_mapping,
        default={},
        help="VASP INCAR overrides for QHA static energies",
    )
    parser.add_argument("--scph_enabled", type=str2bool, nargs="?", const=True, default=False)
    parser.add_argument("--scph_workdir", default="phonopy-sscha")
    parser.add_argument("--scph_temps", type=float, nargs=3, default=[0, 1000, 50])
    parser.add_argument("--scph_initial_fc2", default=None)
    parser.add_argument("--scph_born", default=None)
    parser.add_argument("--scph_snapshots", type=int, default=1000)
    parser.add_argument("--scph_iterations", type=int, default=10)
    parser.add_argument("--scph_transient", type=int, default=1)
    parser.add_argument("--scph_sscha_mesh", type=int, nargs=3, default=[20, 20, 20])
    parser.add_argument("--scph_random_seed", type=int, default=42)
    parser.add_argument("--scph_cutoff_frequency", type=float, default=0.01)
    parser.add_argument(
        "--scph_fc_calculator",
        choices=["symfc", "traditional"],
        default="symfc",
    )
    parser.add_argument("--scph_fc_calculator_options", default=None)
    parser.add_argument(
        "--scph_save_datasets", type=str2bool, nargs="?", const=True, default=False
    )
    parser.add_argument(
        "--scph_run_transport", type=str2bool, nargs="?", const=True, default=False
    )
    parser.add_argument("--scph_transport_fc3", default=None)
    parser.add_argument("--scph_transport_metadata", default=None)
    parser.add_argument("--scph_bubble", type=str2bool, nargs="?", const=True, default=False)
    parser.add_argument("--scph_bubble_mesh", type=int, nargs=3, default=[9, 9, 9])
    parser.add_argument("--scph_bubble_epsilons", type=float, nargs="+", default=[0.05, 0.1])
    parser.add_argument("--scph_bubble_grid_points", type=int, nargs="+", default=None)
    parser.add_argument(
        "--qha_sscha_enabled", type=str2bool, nargs="?", const=True, default=False
    )
    parser.add_argument(
        "--qha_sscha_three_phonon",
        type=str2bool,
        nargs="?",
        const=True,
        default=None,
        help="Run phono3py three-phonon transport after QHA+SSCHA",
    )
    parser.add_argument(
        "--qha_sscha_four_phonon",
        type=str2bool,
        nargs="?",
        const=True,
        default=False,
        help="Run FourPhonon three+four-phonon transport after QHA+SSCHA",
    )

    parser.add_argument("--tdbte_kernel", default=None)
    parser.add_argument("--tdbte_force_constants", default=None)
    parser.add_argument("--tdbte_mesh", type=int, nargs=3, default=None)
    parser.add_argument("--tdbte_temperature", type=float, default=300.)
    parser.add_argument("--tdbte_duration_ps", type=float, default=20.)
    parser.add_argument("--tdbte_max_step_ps", type=float, default=.5)
    parser.add_argument("--tdbte_excitation", type=float, default=.01)
    parser.add_argument("--tdbte_branches", type=int, nargs="+", default=None)
    parser.add_argument("--tdbte_samples", type=int, default=201)
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
        data = load_yaml_strict(handle) or {}
    if not isinstance(data, dict):
        raise ValueError("YAML input must be a mapping of sections and options.")
    validate_input_shape(data)
    data = expand_input_modes(data)

    flat = {}
    sections = yaml_input_sections()
    section_aliases = {
        normalize_yaml_key(section): section for section in sections
    }
    for raw_section, section_data in data.items():
        normalized_section = normalize_yaml_key(raw_section)
        section = section_aliases.get(normalized_section)
        if section is None:
            raise_unknown_yaml_key(str(raw_section), "top level", sections)

        if section_data is None:
            continue
        if not isinstance(section_data, dict):
            raise ValueError(f"YAML section '{raw_section}' must be a mapping.")
        valid_keys = sections[section]
        seen_keys = set()
        for key, value in section_data.items():
            normalized_key = normalize_yaml_key(key)
            if normalized_key in seen_keys:
                raise ValueError(
                    f"Duplicate YAML key '{key}' after '-'/'_' normalization "
                    f"in section '{section}'."
                )
            seen_keys.add(normalized_key)
            if normalized_key not in valid_keys:
                raise_unknown_yaml_key(str(key), f"section '{section}'", valid_keys)
            if section == "structure" and normalized_key == "relaxation":
                if isinstance(value, dict):
                    options = {
                        normalize_yaml_key(name): setting
                        for name, setting in value.items()
                    }
                    if len(options) != len(value):
                        raise ValueError(
                            "Duplicate structure.relaxation option after "
                            "'-'/'_' normalization."
                        )
                    unknown = options.keys() - {"enabled", "workdir", "stages"}
                    if unknown:
                        raise_unknown_yaml_key(
                            sorted(unknown)[0], "structure.relaxation",
                            {"enabled", "workdir", "stages"},
                        )
                    flat["do_relax"] = options.get("enabled", True)
                    if "workdir" in options:
                        flat["vasp_relax_workdir"] = options["workdir"]
                    if "stages" in options:
                        flat["vasp_relax_stages"] = options["stages"]
                elif isinstance(value, bool):
                    flat["do_relax"] = value
                else:
                    raise ValueError(
                        "structure.relaxation must be true, false, or a mapping."
                    )
                continue
            if section == "workflow":
                flat[f"workflow_{normalized_key}"] = value
                continue
            if section == "calculator":
                if normalized_key == "name":
                    flat["calculator"] = value
                elif normalized_key == "factory":
                    flat["calculator_factory"] = value
                elif normalized_key == "kwargs":
                    if not isinstance(value, dict):
                        raise ValueError("calculator.kwargs must be a mapping.")
                    flat["calculator_kwargs"] = value
                elif normalized_key == "model_files":
                    if not isinstance(value, list) or not all(
                        isinstance(path, str) for path in value
                    ):
                        raise ValueError(
                            "calculator.model-files must be a list of paths."
                        )
                    flat["calculator_model_files"] = value
                elif normalized_key == "model":
                    flat["mace_model"] = value
                elif normalized_key == "foundation":
                    flat["mace_foundation"] = value
                elif normalized_key == "device":
                    flat["mace_device"] = value
                elif normalized_key == "dtype":
                    flat["mace_dtype"] = value
                else:
                    flat[normalized_key] = value
                continue
            if section == "force-constant":
                if normalized_key == "workdir":
                    flat["vasp_workdir"] = value
                elif normalized_key == "vasp_kwargs":
                    flat["vasp_kwargs"] = value
                elif normalized_key == "format":
                    flat["fc_format"] = value
                elif normalized_key == "parallel":
                    flat["force_parallel"] = value
                else:
                    flat[normalized_key] = value
                continue
            if section == "kappa" and normalized_key == "command":
                flat["kappa_command"] = value
                continue
            if section == "kappa" and normalized_key == "engine":
                flat["kappa_engine"] = value
                continue
            if section == "kappa" and normalized_key == "qha_volumes":
                flat["qha_volumes"] = value
                continue
            if section == "kappa" and normalized_key == "parallel":
                flat["lbte_parallel"] = value
                continue
            if section == "plot":
                flat[f"plot_{normalized_key}"] = value
                continue
            if section == "tdbte":
                flat[f"tdbte_{normalized_key}"] = value
                continue
            if section == "qha":
                if normalized_key == "enabled":
                    flat["qha_enabled"] = value
                    continue
                flat["qha_enabled"] = True
                flat[f"qha_{normalized_key}"] = value
                continue
            if section == "scph":
                if normalized_key == "enabled":
                    flat["scph_enabled"] = value
                    continue
                flat["scph_enabled"] = True
                flat[f"scph_{normalized_key}"] = value
                continue
            if section == "qha-sscha":
                flat["qha_sscha_enabled"] = True
                flat[f"qha_sscha_{normalized_key}"] = value
                continue
            if section == "fourphonon":
                if normalized_key == "enabled":
                    flat["fp_enabled"] = value
                    continue
                flat["fp_enabled"] = True
                flat[f"fp_{normalized_key}"] = value
                continue
            flat[normalized_key] = value

    args = []
    for key in yaml_arg_order():
        if key in flat:
            append_arg(args, key, flat[key])
    return args


def validate_input_shape(data):
    """Accept only the public static or dynamic YAML shape before expansion."""
    static = {
        "structure", "calculator", "force_constant", "kappa", "plot", "output"
    }
    dynamic = {"tdbte", "output"}
    optional = {"parallel"}
    names = {normalize_yaml_key(name): name for name in data}
    if len(names) != len(data):
        raise ValueError("Duplicate top-level section after '-'/'_' normalization.")
    if "tdbte" in names:
        unexpected = names.keys() - dynamic
        if unexpected:
            if unexpected & static:
                raise ValueError(
                    "Dynamic tdbte input cannot mix static sections: "
                    + ", ".join(sorted(unexpected)) + "."
                )
            raise ValueError(
                "Dynamic input accepts only tdbte and output; remove: "
                + ", ".join(sorted(unexpected)) + "."
            )
        missing = dynamic - names.keys()
        if missing:
            raise ValueError(
                "Dynamic input requires: " + ", ".join(sorted(missing)) + "."
            )
    else:
        unexpected = names.keys() - static - optional
        if unexpected:
            raise_unknown_yaml_key(
                str(names[sorted(unexpected)[0]]),
                "top level",
                {
                    "structure", "calculator", "force-constant", "kappa",
                    "plot", "output", "parallel",
                },
            )
        missing = static - names.keys()
        if missing:
            raise ValueError(
                "Static input requires six sections; missing: "
                + ", ".join(name.replace("_", "-") for name in sorted(missing))
                + "."
            )
    for original in names.values():
        if not isinstance(data[original], dict):
            raise ValueError(f"YAML section '{original}' must be a mapping.")


def expand_input_modes(data):
    """Compile the six-section static form or the separate TD-BTE form."""
    data = _expand_static_parallel(data)
    keys = {normalize_yaml_key(key): key for key in data}
    if "tdbte" in keys and "workflow" not in keys:
        static = {"structure", "calculator", "force_constant", "kappa", "plot"}
        mixed = sorted(static & keys.keys())
        if mixed:
            raise ValueError(
                "Dynamic tdbte input cannot mix static sections: "
                + ", ".join(mixed)
                + "."
            )
        return {"workflow": {"preset": "custom", "steps": ["tdbte"]}, **data}

    fc_key = keys.get("force_constant")
    if fc_key is None or not isinstance(data[fc_key], dict):
        return data
    fc = data[fc_key]
    kappa = data.get(keys.get("kappa"), {})
    kappa_values = (
        {normalize_yaml_key(key): value for key, value in kappa.items()}
        if isinstance(kappa, dict)
        else {}
    )
    engine = kappa_values.get("engine")
    if engine is not None and engine not in {"phono3py", "fourphonon"}:
        raise ValueError("kappa.engine must be phono3py or fourphonon.")
    separate_methods = {"method_3ph", "method_4ph"} & kappa_values.keys()
    if separate_methods and not {"method_3ph", "method_4ph"} <= kappa_values.keys():
        raise ValueError(
            "Set both kappa.method-3ph and kappa.method-4ph for four-phonon transport."
        )
    switch_keys = {
        normalize_yaml_key(key): key
        for key in fc
        if normalize_yaml_key(key) in {"qha", "sscha", "four_phonon"}
    }
    if not switch_keys:
        if separate_methods:
            raise ValueError(
                "kappa.method-3ph and kappa.method-4ph require "
                "force-constant.four-phonon: true."
            )
        if engine == "fourphonon":
            raise ValueError(
                "kappa.engine: fourphonon requires force-constant.four-phonon: true."
            )
        return data
    incompatible = {"workflow", "qha", "scph", "qha_sscha", "fourphonon"}
    mixed = sorted(incompatible & keys.keys())
    if mixed:
        raise ValueError(
            "Nested force-constant switches cannot be combined with top-level "
            + ", ".join(mixed)
            + "."
        )

    sections = yaml_input_sections()
    qha_on, qha_options = _nested_static_switch(
        fc.get(switch_keys.get("qha")), "qha", sections["qha"]
    )
    sscha_on, sscha_options = _nested_static_switch(
        fc.get(switch_keys.get("sscha")), "sscha", sections["scph"]
    )
    fc4_keys = {"dim_fc4", "cutoff_fc4", "fourthorder_command", "fc4_workdir"}
    four_on, four_options = _nested_static_switch(
        fc.get(switch_keys.get("four_phonon")),
        "four-phonon",
        sections["fourphonon"] | fc4_keys,
    )
    if engine == "fourphonon" and not four_on:
        raise ValueError(
            "kappa.engine: fourphonon requires force-constant.four-phonon: true."
        )
    if engine == "phono3py" and four_on:
        raise ValueError(
            "force-constant.four-phonon: true requires kappa.engine: fourphonon."
        )
    if separate_methods:
        if not four_on:
            raise ValueError(
                "kappa.method-3ph and kappa.method-4ph require "
                "force-constant.four-phonon: true."
            )
        if "method" in kappa_values:
            raise ValueError(
                "Use kappa.method-3ph and kappa.method-4ph without kappa.method."
            )
        methods = (kappa_values["method_3ph"], kappa_values["method_4ph"])
        solver_by_methods = {
            ("rta", "rta"): "rta",
            ("lbte", "rta"): "3ph-iterative",
            ("lbte", "lbte"): "full-iterative",
        }
        if any(method not in {"rta", "lbte"} for method in methods):
            raise ValueError("kappa.method-3ph and kappa.method-4ph must be rta or lbte.")
        if methods not in solver_by_methods:
            raise ValueError(
                "FourPhonon does not support 3ph RTA with 4ph LBTE; "
                "use 3ph LBTE when 4ph uses LBTE."
            )
        selected_solver = solver_by_methods[methods]
        if "solver" in four_options and four_options["solver"] != selected_solver:
            raise ValueError(
                "kappa.method-3ph/method-4ph conflict with "
                "force-constant.four-phonon.solver."
            )
        four_options["solver"] = selected_solver
    qha_volumes = kappa_values.get("qha_volumes", False)
    if not isinstance(qha_volumes, bool):
        raise ValueError("kappa.qha-volumes must be true or false.")
    if qha_volumes and (not qha_on or sscha_on or four_on):
        raise ValueError(
            "kappa.qha-volumes requires QHA without SSCHA or four-phonon transport."
        )
    if qha_volumes and kappa_values.get("method", "rta") != "rta":
        raise ValueError("kappa.qha-volumes currently requires kappa.method: rta.")
    if four_on and qha_on and not sscha_on:
        raise ValueError(
            "QHA plus four-phonon transport requires SSCHA; "
            "QHA alone does not feed corrected FC2 into four-phonon transport."
        )
    if four_on and sscha_on and not qha_on:
        raise ValueError(
            "Standalone SSCHA plus four-phonon transport is not supported; "
            "enable QHA as well for the QHA-volume SSCHA route."
        )
    if four_on and isinstance(kappa, dict):
        if kappa_values.get("method", "rta") != "rta":
            raise ValueError(
                "Four-phonon static input requires kappa.method: rta; "
                "choose a FourPhonon solver inside force-constant.four-phonon."
            )
        try:
            boundary_mfp = float(kappa_values.get("bfmp", 1.0e6))
        except (TypeError, ValueError):
            boundary_mfp = None  # The normal argparse parser reports invalid input.
        if boundary_mfp is not None and boundary_mfp != 1.0e6:
            raise ValueError(
                "Four-phonon static input does not apply kappa.bfmp; "
                "remove the non-default boundary setting."
            )
        for source, target in (("wigner", "wigner"), ("isotope", "isotopes")):
            if (
                source in kappa_values
                and target in four_options
                and kappa_values[source] != four_options[target]
            ):
                raise ValueError(
                    f"kappa.{source} conflicts with "
                    f"force-constant.four-phonon.{target}."
                )
            if kappa_values.get(source) is True:
                four_options[target] = True
    structure = data.get(keys.get("structure"), {})
    relaxation = structure.get("relaxation") if isinstance(structure, dict) else None
    relax_requested = relaxation is True or (
        isinstance(relaxation, dict) and relaxation.get("enabled", True) is True
    )
    if qha_on and relax_requested:
        raise ValueError(
            "QHA reads structure.poscar directly; set structure.relaxation: false "
            "and provide the intended starting structure."
        )

    result = dict(data)
    fc_detail = {
        key: value for key, value in fc.items()
        if normalize_yaml_key(key) not in switch_keys
    }
    if four_on:
        for key, value in list(four_options.items()):
            if key in fc4_keys:
                if any(normalize_yaml_key(existing) == key for existing in fc_detail):
                    raise ValueError(
                        f"force-constant.four-phonon.{key} duplicates "
                        f"force-constant.{key}."
                    )
                fc_detail[key.replace("_", "-")] = value
                del four_options[key]
        if not any(normalize_yaml_key(key) in {"format", "fc_format"} for key in fc_detail):
            fc_detail["format"] = "both"
        elif not qha_on:
            selected_format = next(
                value for key, value in fc_detail.items()
                if normalize_yaml_key(key) in {"format", "fc_format"}
            )
            if selected_format not in {"both", "shengbte"}:
                raise ValueError(
                    "Four-phonon transport requires force-constant.format: "
                    "both or shengbte."
                )
    result[fc_key] = fc_detail
    if separate_methods:
        result[keys["kappa"]] = {
            key: value for key, value in kappa.items()
            if normalize_yaml_key(key) not in separate_methods
        }
    if qha_on:
        result["qha"] = {"enabled": True, **qha_options}
    if sscha_on:
        result["scph"] = {"enabled": True, **sscha_options}
    if four_on:
        result["fourphonon"] = {"enabled": True, **four_options}
    if qha_on and sscha_on:
        result["qha-sscha"] = {
            "three-phonon": sscha_options.get("run_transport", False),
            "four-phonon": four_on,
        }
        preset = "qha-sscha"
    elif qha_on:
        preset = "qha-kappa" if qha_volumes else "qha"
    elif sscha_on:
        preset = "scph"
    elif four_on:
        preset = "four-phonon"
    else:
        preset = "three-phonon"
    if sscha_on and not qha_on and relax_requested:
        result["workflow"] = {
            "preset": "custom", "steps": ["relax", "fc2fc3", "scph"]
        }
    else:
        result["workflow"] = {"preset": preset}
    return result


def _expand_static_parallel(data):
    """Move optional execution settings to their existing stage-specific slots."""
    keys = {normalize_yaml_key(key): key for key in data}
    parallel_key = keys.get("parallel")
    if parallel_key is None:
        return data
    parallel = data[parallel_key]
    if not isinstance(parallel, dict):
        raise ValueError("parallel must be a mapping.")
    result = {key: value for key, value in data.items() if key != parallel_key}
    fc_key = keys.get("force_constant")
    kappa_key = keys.get("kappa")
    if fc_key is None or kappa_key is None:
        raise ValueError("parallel requires force-constant and kappa sections.")
    fc = dict(result[fc_key])
    kappa = dict(result[kappa_key])
    for raw_target, settings in parallel.items():
        target = normalize_yaml_key(raw_target)
        if target not in {"force_constant", "kappa"}:
            raise_unknown_yaml_key(str(raw_target), "parallel", {"force_constant", "kappa"})
        if not isinstance(settings, dict):
            raise ValueError(f"parallel.{raw_target} must be a mapping.")
        if target == "force_constant":
            if any(normalize_yaml_key(key) == "parallel" for key in fc):
                raise ValueError("parallel.force-constant duplicates force-constant.parallel.")
            fc["parallel"] = settings
            continue
        switches = {normalize_yaml_key(key): key for key in fc}
        four_key = switches.get("four_phonon")
        four = fc.get(four_key) if four_key else None
        four_on = four is True or (isinstance(four, dict) and four.get("enabled", True) is True)
        engine = next((value for key, value in kappa.items() if normalize_yaml_key(key) == "engine"), None)
        if engine == "fourphonon" or (engine is None and four_on):
            if not four_on:
                raise ValueError("parallel.kappa for FourPhonon requires force-constant.four-phonon: true.")
            four = {} if four is True else dict(four)
            special = {"mpi_launcher", "mpi_processes", "omp_threads", "omp_stacksize"}
            scheduler = {}
            for key, value in settings.items():
                normalized = normalize_yaml_key(key)
                if normalized in special:
                    if any(normalize_yaml_key(existing) == normalized for existing in four):
                        raise ValueError(f"parallel.kappa.{key} duplicates force-constant.four-phonon.{key}.")
                    four[key] = value
                else:
                    scheduler[key] = value
            if scheduler:
                if any(normalize_yaml_key(key) == "parallel" for key in four):
                    raise ValueError("parallel.kappa duplicates force-constant.four-phonon.parallel.")
                four["parallel"] = scheduler
            fc[four_key] = four
        else:
            if any(normalize_yaml_key(key) == "parallel" for key in kappa):
                raise ValueError("parallel.kappa duplicates kappa.parallel.")
            kappa["parallel"] = settings
    result[fc_key] = fc
    result[kappa_key] = kappa
    return result


def _nested_static_switch(raw, name, allowed):
    if raw is None:
        return False, {}
    if isinstance(raw, bool):
        return raw, {}
    if not isinstance(raw, dict):
        raise ValueError(f"force-constant.{name} must be true, false, or a mapping.")
    options = {}
    for key, value in raw.items():
        normalized = normalize_yaml_key(key)
        if normalized not in allowed | {"enabled"}:
            raise_unknown_yaml_key(
                str(key), f"force-constant.{name}", allowed | {"enabled"}
            )
        if normalized in options:
            raise ValueError(f"Duplicate force-constant.{name} key '{key}'.")
        options[normalized] = value
    enabled = options.pop("enabled", True)
    if not isinstance(enabled, bool):
        raise ValueError(f"force-constant.{name}.enabled must be true or false.")
    if not enabled and options:
        raise ValueError(
            f"force-constant.{name} is disabled but still configures options."
        )
    return enabled, options


def raise_unknown_yaml_key(key, location, valid_keys):
    """Raise a helpful error for an unsupported YAML option."""
    normalized_valid = sorted({str(item).replace("_", "-") for item in valid_keys})
    matches = difflib.get_close_matches(
        str(key).replace("_", "-"), normalized_valid, n=1, cutoff=0.6
    )
    suggestion = f" Did you mean '{matches[0]}'?" if matches else ""
    raise ValueError(f"Unknown YAML key '{key}' in {location}.{suggestion}")


def load_yaml_strict(stream):
    """Load YAML while rejecting duplicate mapping keys at every level."""
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


def yaml_input_sections():
    """Return supported YAML sections and their documented keys."""
    return {
        "workflow": {"preset", "steps"},
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
            "engine",
            "qha_volumes",
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
        },
        "qha": {
            "enabled",
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


def yaml_arg_order():
    """Return a stable option order for parsed YAML values."""
    return [
        "tdbte_kernel", "tdbte_force_constants", "tdbte_mesh", "tdbte_temperature",
        "tdbte_duration_ps", "tdbte_max_step_ps", "tdbte_excitation",
        "tdbte_branches", "tdbte_samples",
        "workflow_preset",
        "workflow_steps",
        "poscar",
        "dimensionality",
        "effective_thickness",
        "effective_area",
        "vacuum_axis",
        "periodic_axis",
        "nep_model",
        "calculator",
        "calculator_factory",
        "calculator_kwargs",
        "calculator_model_files",
        "mace_model",
        "mace_foundation",
        "mace_device",
        "mace_dtype",
        "vasp_command",
        "vasp_path",
        "potcar_path",
        "vasp_workdir",
        "vasp_kwargs",
        "vasp_relax_workdir",
        "vasp_relax_stages",
        "do_relax",
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
        "fc_format",
        "n_structures",
        "rattle_std",
        "cutoffs",
        "min_dist",
        "force_parallel",
        "mesh",
        "temps",
        "method",
        "kappa_engine",
        "qha_volumes",
        "kappa_command",
        "isotope",
        "bfmp",
        "wigner",
        "lbte_parallel",
        "fp_enabled",
        "fp_command",
        "fp_workdir",
        "fp_control",
        "fp_fc2",
        "fp_fc3",
        "fp_fc4",
        "fp_harmonic_format",
        "fp_mesh",
        "fp_temps",
        "fp_scell",
        "fp_solver",
        "fp_wigner",
        "fp_scalebroad",
        "fp_isotopes",
        "fp_nonanalytic",
        "fp_only_harmonic",
        "fp_sample_3ph",
        "fp_sample_3ph_phase_space",
        "fp_sample_4ph",
        "fp_sample_4ph_phase_space",
        "fp_mpi_launcher",
        "fp_mpi_processes",
        "fp_omp_threads",
        "fp_omp_stacksize",
        "fp_parallel",
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
        "qha_enabled",
        "qha_volume_ratios",
        "qha_dim_fc2",
        "qha_mesh",
        "qha_temps",
        "qha_eos",
        "qha_pressure",
        "qha_relax_internal",
        "qha_relax_fmax",
        "qha_relax_steps",
        "qha_displacement_distance",
        "qha_imaginary_frequency_tolerance",
        "qha_cutoff_frequency",
        "qha_vasp_relax_kwargs",
        "qha_vasp_static_kwargs",
        "scph_enabled",
        "scph_workdir",
        "scph_temps",
        "scph_initial_fc2",
        "scph_born",
        "scph_snapshots",
        "scph_iterations",
        "scph_transient",
        "scph_sscha_mesh",
        "scph_random_seed",
        "scph_cutoff_frequency",
        "scph_fc_calculator",
        "scph_fc_calculator_options",
        "scph_save_datasets",
        "scph_run_transport",
        "scph_transport_fc3",
        "scph_transport_metadata",
        "scph_bubble",
        "scph_bubble_mesh",
        "scph_bubble_epsilons",
        "scph_bubble_grid_points",
        "qha_sscha_enabled",
        "qha_sscha_three_phonon",
        "qha_sscha_four_phonon",
    ]


def normalize_yaml_key(key):
    """Normalize YAML keys to argparse option names."""
    return str(key).strip().replace("-", "_")


def append_arg(args, key, value):
    """Append one parsed YAML option to an argparse token list."""
    if value is None:
        return
    if key == "calculator_model_files" and not value:
        return
    if key == "fp_mpi_launcher":
        args.extend([f"--{key}", json.dumps(value)])
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


WORKFLOW_PRESETS = {
    "three-phonon": ["relax", "fc2fc3", "kappa"],
    "four-phonon": ["relax", "fc2fc3", "fc4", "kappa4"],
    "qha": ["qha"],
    "qha-kappa": ["qha", "qha-kappa"],
    "qha-sscha": ["qha", "qha-sscha"],
}
WORKFLOW_STEPS = {
    "relax", "fc2", "fc2fc3", "fc4", "qha", "qha-kappa", "scph", "qha-sscha",
    "kappa", "kappa4", "plot", "report", "bubble", "tdbte",
}


def resolve_workflow_plan(args, parser=None):
    """Resolve a high-level preset into the ordered stages used by ``run``."""
    preset = args.workflow_preset
    explicit = args.workflow_steps
    error = None
    if preset == "custom":
        if not explicit:
            error = "workflow.preset: custom requires a non-empty workflow.steps list."
        steps = list(explicit or [])
    elif explicit:
        error = "workflow.steps may be used only with workflow.preset: custom."
        steps = []
    elif preset == "scph":
        steps = ["fc2fc3", "scph"]
    else:
        steps = list(WORKFLOW_PRESETS[preset])

    unknown = [step for step in steps if step not in WORKFLOW_STEPS]
    if unknown:
        error = f"Unknown workflow step '{unknown[0]}'."
    elif len(steps) != len(set(steps)):
        error = "workflow.steps must not contain duplicate stages."
    elif "fc2" in steps and "fc2fc3" in steps:
        error = "workflow.steps cannot contain both fc2 and fc2fc3."
    if error:
        if parser is not None:
            parser.error(error)
        raise ValueError(error)
    args.workflow_steps = steps
    return steps


def parse_workflow_args(config_path, command=None):
    """Parse workflow configuration from a YAML file."""
    parser = initialise_parser()
    if not Path(config_path).exists():
        parser.error(f"input file not found: {config_path}")
    tokens = parse_input_file(config_path)
    args = parser.parse_args(tokens, namespace=WorkflowConfig())
    args.config_path = str(Path(config_path).resolve())
    cutoff = args.cutoff_fc3
    if cutoff is None and args.fc3_backend == "thirdorder":
        cutoff = -3.0
    args.cutoff_fc3 = cutoff
    resolve_force_constant_dimensions(args)
    resolve_workflow_plan(args, parser=parser)
    targets = set(args.workflow_steps) if command == "run" else {command}
    if args.kappa_engine == "fourphonon" and "kappa" in targets:
        parser.error("kappa.engine: fourphonon cannot run the phono3py kappa stage.")
    if args.kappa_engine == "phono3py" and "kappa4" in targets:
        parser.error("kappa.engine: phono3py cannot run the FourPhonon kappa4 stage.")
    if command == "run" and args.qha_volumes and "qha-kappa" not in targets:
        parser.error("kappa.qha-volumes: true requires the qha-kappa workflow.")
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
    if command is None or targets & {"kappa", "qha-kappa"}:
        temps_error = validate_temps(args.temps)
        if temps_error:
            parser.error(temps_error)
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
            parser.error("qha-kappa requires force-constant.qha.enabled and kappa.qha-volumes: true.")
        if args.method != "rta" or args.kappa_engine == "fourphonon":
            parser.error("qha-kappa currently supports phono3py RTA only.")
        if str((args.force_parallel or {}).get("backend", "none")).lower() == "slurm":
            parser.error("qha-kappa runs its force calculations inside one allocation; use parallel.force-constant.backend: none.")
        if len(args.temps) == 1:
            requested = args.temps
        else:
            minimum, maximum, step = args.temps
            if minimum < 0 or maximum < minimum or step <= 0:
                parser.error("qha-kappa requires 0 <= minimum <= maximum and step > 0 in kappa.temps.")
            requested = [minimum, maximum]
        if len(args.qha_temps) != 3:
            parser.error("force-constant.qha.temps must contain [minimum, maximum, step].")
        qha_minimum, qha_maximum, _ = args.qha_temps
        if min(requested) < qha_minimum or max(requested) > qha_maximum:
            parser.error("kappa.temps must lie within the QHA temperature range.")
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
    return args


def resolve_force_constant_dimensions(args):
    """Fill omitted FC dimensions from the standard defaults."""
    default_dim = [4, 4, 1]
    if args.dim_fc2 is None:
        args.dim_fc2 = list(default_dim)
    if args.dim_fc3 is None:
        args.dim_fc3 = list(default_dim)
    if args.dim_fc4 is None:
        args.dim_fc4 = list(args.dim_fc3)


def iter_display_args(args, command=None):
    """Iterate over user-facing config values, hiding inactive route settings."""
    compatibility_only = {"config_path", "invoked_command"}
    hiphive_only = {"n_structures", "rattle_std", "cutoffs", "min_dist"}
    thirdorder_only = {"thirdorder_command", "fc3_workdir"}
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
    external_calculator_only = {
        "calculator_factory",
        "calculator_kwargs",
        "calculator_model_files",
    }
    mace_only = {"mace_model", "mace_foundation", "mace_device", "mace_dtype"}
    film_only = {"effective_thickness", "vacuum_axis"}
    wire_only = {"effective_area", "periodic_axis"}
    qha_sscha_only = {
        name for name in vars(args) if name.startswith("qha_sscha_")
    }
    qha_only = {
        name
        for name in vars(args)
        if name.startswith("qha_") and name not in qha_sscha_only
    }
    scph_only = {name for name in vars(args) if name.startswith("scph_")}
    fourphonon_only = {name for name in vars(args) if name.startswith("fp_")}
    common = {
        "poscar",
        "dimensionality",
        "effective_thickness",
        "effective_area",
        "vacuum_axis",
        "periodic_axis",
        "progress",
        "result_dir",
    }
    plan_fields = {"workflow_preset", "workflow_steps"}
    calculator_fields = {
        "nep_model",
        "calculator",
        "calculator_factory",
        "calculator_kwargs",
        "calculator_model_files",
        "mace_model",
        "mace_foundation",
        "mace_device",
        "mace_dtype",
        "vasp_command",
        "vasp_path",
        "potcar_path",
        "vasp_workdir",
        "vasp_kwargs",
        "vasp_relax_workdir",
        "vasp_relax_stages",
        "do_relax",
    }
    fc2_fields = {
        "dim_fc2", "use_hiphive", "compact_fc", "fc_format",
        "force_parallel",
    }
    fc3_fields = {
        "dim_fc3", "fc3_backend", "cutoff_fc3", "thirdorder_command",
        "fc3_workdir",
    }
    fc4_fields = {
        "dim_fc4", "cutoff_fc4", "fourthorder_command", "fc4_workdir",
    }
    kappa_fields = {
        "mesh", "temps", "method", "kappa_engine", "qha_volumes", "kappa_command", "isotope", "bfmp",
        "wigner", "lbte_parallel",
    }
    command_fields = None
    tdbte_fields = {name for name in vars(args) if name.startswith("tdbte_")}
    if command == "tdbte":
        command_fields = {"result_dir"} | plan_fields | tdbte_fields
    elif command == "run":
        command_fields = (
            common | plan_fields | calculator_fields | fc2_fields | fc3_fields
            | fc4_fields | kappa_fields | qha_only | scph_only
            | qha_sscha_only | fourphonon_only | tdbte_fields
        )
    elif command == "relax":
        command_fields = common | calculator_fields
    elif command == "fc2":
        command_fields = common | calculator_fields | fc2_fields | hiphive_only
    elif command == "fc2fc3":
        command_fields = (
            common | calculator_fields | fc2_fields | fc3_fields | hiphive_only
        )
    elif command == "fc4":
        command_fields = common | calculator_fields | fc4_fields
    elif command == "scph":
        command_fields = common | scph_only | calculator_fields
        if args.scph_run_transport:
            command_fields |= kappa_fields
    elif command == "bubble":
        command_fields = common | scph_only | {"workflow_preset"}
    elif command == "qha":
        command_fields = common | calculator_fields | qha_only
    elif command == "qha-sscha":
        command_fields = (
            common | calculator_fields | fc2_fields | fc3_fields | hiphive_only
            | fc4_fields | kappa_fields | qha_only | scph_only
            | qha_sscha_only | fourphonon_only
        )
    elif command == "qha-kappa":
        command_fields = (
            common | calculator_fields | fc2_fields | fc3_fields
            | hiphive_only | kappa_fields | qha_only
        )
    elif command == "kappa4":
        command_fields = common | fourphonon_only
    elif command == "kappa":
        command_fields = common | kappa_fields
    elif command == "plot":
        command_fields = common | {"mesh", "temps"} | {
            name for name in vars(args) if name.startswith("plot_")
        }
    if command == "run" and args.workflow_steps == ["tdbte"]:
        command_fields = {"result_dir"} | plan_fields | tdbte_fields
    for arg, value in vars(args).items():
        if value is None:
            continue
        if command_fields is not None and arg not in command_fields:
            continue
        if arg == "lbte_parallel" and not value:
            continue
        if arg == "force_parallel" and not value:
            continue
        if arg in {"calculator_kwargs", "calculator_model_files"} and not value:
            continue
        if arg in compatibility_only:
            continue
        if args.dimensionality != 2 and arg in film_only:
            continue
        if args.dimensionality != 1 and arg in wire_only:
            continue
        if not args.qha_enabled and arg in qha_only:
            continue
        if not args.scph_enabled and arg in scph_only:
            continue
        if (
            command != "qha-sscha"
            and not args.qha_sscha_enabled
            and arg in qha_sscha_only
        ):
            continue
        if not args.fp_enabled and arg in fourphonon_only:
            continue
        if args.calculator == "vasp" and arg == "nep_model" and value is None:
            continue
        if not args.use_hiphive and arg in hiphive_only:
            continue
        if args.fc3_backend != "thirdorder" and arg in thirdorder_only:
            continue
        if args.calculator != "vasp" and arg in vasp_only:
            continue
        if args.calculator in {"nep", "vasp"} and arg in external_calculator_only:
            continue
        if args.calculator != "mace" and arg in mace_only:
            continue
        if (args.calculator != "vasp" or not args.do_relax) and arg in vasp_relax_only:
            continue
        yield arg, value


def format_config(args, command=None):
    """Return a readable multi-line configuration summary."""
    lines = ["Running Workflow with configuration:"]
    for arg, value in iter_display_args(args, command=command):
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
    return None


def validate_fourphonon(args):
    """Validate FourPhonon 3ph+4ph transport settings."""
    if not getattr(args, "fp_enabled", False):
        return None
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
        return temps_error.replace("--temps", "force-constant.four-phonon.temps")
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
        return "--nep_model is required when --calculator nep"
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
            "calculator.model, foundation, device, and dtype are native MACE "
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
    """Validate optional Slurm settings for distributed LBTE calculations."""
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
    if args.method != "lbte":
        return "kappa.parallel.backend 'slurm' requires kappa.method: lbte."

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
            "LBTE grid points across the Slurm array."
        )

    for key in ("preamble", "extra_sbatch"):
        value = settings.get(key, [])
        if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
            return f"kappa.parallel.{key} must be a list of strings."
    if "submit" in settings and not isinstance(settings["submit"], bool):
        return "kappa.parallel.submit must be true or false."
    return None
