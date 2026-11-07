"""Configuration parsing utilities for NEP-kappa."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from nepkappa.config_models import WorkflowConfig
from nepkappa.config_schema import (
    str2bool, json_dict, json_list, normalize_yaml_key,
    raise_unknown_yaml_key, yaml_input_sections, load_yaml_strict,
)
from nepkappa.config_plan import (
    THERMAL_SECTIONS, DYNAMIC_SECTIONS, STAGE_REQUIRED_SECTIONS,
    expand_input_modes, resolve_workflow_plan,
)
from nepkappa.config_validation import (
    validate_workflow_args, validate_temps, expand_transport_temperatures,
)
from nepkappa.command_sampling import phono3py_sampling

def initialise_parser() -> argparse.ArgumentParser:
    """Build the internal workflow argument parser used by YAML input files."""
    parser = argparse.ArgumentParser(
        prog="nepkappa run",
        description="NEP-kappa calculation stages selected by the input file",
    )

    parser.add_argument(
        "--workflow_steps",
        nargs="+",
        default=None,
        help="Internal ordered stages derived from input switches",
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
    parser.add_argument("--nep_force_backend", choices=["cpu", "gpu-batch"], default="cpu")
    parser.add_argument("--nep_batch_size", type=int, default=128)
    parser.add_argument("--nep_command", default="nep")
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
        type=json_dict,
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
        "--kappa_force_constants", default=None,
        help="Existing FC2/FC3 directory; defaults to output.result-dir",
    )
    parser.add_argument(
        "--kappa_engine", choices=["phono3py", "fourphonon"], default=None,
        help="Transport engine selected by kappa.engine",
    )
    parser.add_argument(
        "--kappa_enabled", type=str2bool, nargs="?", const=True, default=True,
        help="Compute thermal conductivity after the enabled force-constant stages",
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
        help="[Kappa] Slurm parallel RTA/LBTE calculation settings",
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
    parser.add_argument("--fp_compare_3ph", type=str2bool, nargs="?", const=True, default=False)
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
        choices=["total", "phph", "isotope", "normal", "umklapp", "nu", "all"],
        default="total",
        help="[Plot] Scattering rates: total compares phph, isotope and their sum; phph is N+U; nu separates N/U",
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
    parser.add_argument("--plot_tdbte_times", type=float, nargs=3, default=None,
                        help="[Plot] Three TD-BTE snapshot times in ps; default: start, midpoint, end")
    parser.add_argument("--plot_scattering_rate_min", type=float, default=0.0,
                        help="[Plot] Omit scattering rates below this value in ps^-1 (display only)")
    parser.add_argument(
        "--qha_enabled",
        type=str2bool,
        nargs="?",
        const=True,
        default=False,
        help="Enable quasi-harmonic approximation settings",
    )
    parser.add_argument(
        "--qha_run_transport", type=str2bool, nargs="?", const=True,
        default=False, help="Use QHA equilibrium volumes for transport",
    )
    parser.add_argument(
        "--qha_equilibrium_phonons", type=str2bool, nargs="?", const=True,
        default=False, help="Recompute FC2 at QHA equilibrium volumes for phonon plots",
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
        type=json_dict,
        default={},
        help="VASP INCAR overrides for fixed-cell QHA relaxation",
    )
    parser.add_argument(
        "--qha_vasp_static_kwargs",
        type=json_dict,
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


def read_input_file(filename, parser):
    """Read one YAML document; the parser supplies consistent CLI errors."""
    path = Path(filename)
    if not path.exists():
        parser.error(f"input file not found: {filename}")
    if path.suffix.lower() not in (".yaml", ".yml"):
        raise ValueError("NEP-kappa input files must be YAML files ending in .yaml or .yml.")
    with path.open("r", encoding="utf-8") as handle:
        data = load_yaml_strict(handle) or {}
    if not isinstance(data, dict):
        raise ValueError("YAML input must be a mapping of sections and options.")
    return data


def yaml_to_tokens(data, command=None):
    """Adapt YAML to the existing argparse type/default rules.

    This small boundary deliberately preserves legacy scalar coercion, boolean
    spellings and nargs semantics instead of maintaining a second type system.
    """
    validate_input_shape(data, command=command)
    data = expand_input_modes(data, command=command)

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
        if section == "calculator":
            calculator_options = {
                normalize_yaml_key(key): value for key, value in section_data.items()
            }
            calculator_name = str(calculator_options.get("name", "nep")).strip().lower()
            if {"model", "nep_model"} <= calculator_options.keys():
                raise ValueError(
                    "Use only calculator.model; do not also set its legacy alias "
                    "calculator.nep-model."
                )
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
            if section in {"kappa", "tdbte"} and normalized_key == "force_constants":
                if not isinstance(value, str) or not value.strip():
                    raise ValueError(f"{section}.force-constants must be a nonempty directory path.")
                flat[f"{section}_force_constants"] = value
                continue
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
                    if calculator_name not in {"nep", "mace"}:
                        raise ValueError(
                            "calculator.model requires calculator.name: nep or mace. "
                            "For ASE or plugin calculators, pass model options in calculator.kwargs."
                        )
                    flat[f"{calculator_name}_model"] = value
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
            if section == "kappa" and normalized_key == "enabled":
                flat["kappa_enabled"] = value
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
    for key, value in flat.items():
        append_arg(args, key, value)
    return args


def validate_input_shape(data, command=None):
    """Check public sections before deriving a full or single-stage plan."""
    names = {normalize_yaml_key(name): name for name in data}
    if len(names) != len(data):
        raise ValueError("Duplicate top-level section after '-'/'_' normalization.")
    allowed = THERMAL_SECTIONS | DYNAMIC_SECTIONS | {"parallel"}
    unexpected = names.keys() - allowed
    if unexpected:
        raise_unknown_yaml_key(
            str(names[sorted(unexpected)[0]]), "top level",
            {name.replace("_", "-") for name in allowed},
        )
    if {"kappa", "tdbte"} <= names.keys():
        raise ValueError("An input cannot contain both kappa and tdbte sections.")
    for original in names.values():
        if not isinstance(data[original], dict):
            raise ValueError(f"YAML section '{original}' must be a mapping.")
    if command in {None, "run", "info"}:
        required = DYNAMIC_SECTIONS if "tdbte" in names else THERMAL_SECTIONS
        label = "nepkappa run requires six sections"
    else:
        required = STAGE_REQUIRED_SECTIONS[command]
        label = f"nepkappa {command} requires input sections"
    missing = required - names.keys()
    if missing:
        raise ValueError(
            label + "; missing: "
            + ", ".join(name.replace("_", "-") for name in sorted(missing)) + "."
        )
    if not names:
        raise ValueError("Input must contain at least one configuration section.")


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


def _read_config(config_path, command):
    """Parse and resolve shared defaults once, before command validation."""
    parser = initialise_parser()
    data = read_input_file(config_path, parser)
    tokens = yaml_to_tokens(data, command=command)
    args = parser.parse_args(tokens, namespace=WorkflowConfig())
    args.config_path = str(Path(config_path).resolve())
    sampling = phono3py_sampling(args.kappa_command) if args.kappa_command else {}
    for field, values in sampling.items():
        if field == "temps":
            error = validate_temps(values)
            if error:
                raise ValueError(error.replace("kappa.temps", "kappa.command temperatures"))
        if f"--{field}" in tokens:
            configured = getattr(args, field)
            if field == "temps":
                error = validate_temps(configured)
                if error:
                    raise ValueError(error)
                configured = expand_transport_temperatures(configured)
                compared = expand_transport_temperatures(values)
            else:
                compared = values
            if len(configured) != len(compared) or any(
                not math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-10)
                for a, b in zip(configured, compared)
            ):
                raise ValueError(f"kappa.{field} conflicts with kappa.command; "
                                 f"remove kappa.{field} or make the values agree.")
        setattr(args, field, values)
    cutoff = args.cutoff_fc3
    if cutoff is None and args.fc3_backend == "thirdorder":
        cutoff = -3.0
    args.cutoff_fc3 = cutoff
    resolve_force_constant_dimensions(args)
    resolve_workflow_plan(args, parser=parser)
    explicit_mesh = "--mesh" in tokens or "mesh" in sampling
    return parser, data, args, explicit_mesh


def _validate_command(args, command, parser, *, explicit_mesh):
    if command == "plot" and not explicit_mesh:
        # A default calculation mesh must not constrain existing-result plots.
        args.mesh = None
    if command == "report":
        from nepkappa.report import resolve_result_directory

        args.result_dir = str(resolve_result_directory(args.config_path))
    validate_workflow_args(args, command, parser)
    return args


def parse_workflow_args(config_path, command=None):
    """Read, plan and validate a YAML file for one concrete operation."""
    parser, _, args, explicit_mesh = _read_config(config_path, command)
    return _validate_command(args, command, parser, explicit_mesh=explicit_mesh)


def parse_selected_command(config_path, command):
    """Resolve public fc2/kappa routes without rereading or reparsing YAML."""
    from copy import copy
    from nepkappa.command_registry import resolve_calculation_command

    parser, data, args, explicit_mesh = _read_config(config_path, command)
    selected = command
    if command in {"fc2", "kappa"}:
        # Preserve the existing preflight checks and their error precedence.
        # The temporary plot mesh must not replace the calculation default.
        plan = _validate_command(copy(args), "plot", parser, explicit_mesh=explicit_mesh)
        try:
            selected = resolve_calculation_command(command, plan)
        except ValueError as exc:
            raise SystemExit(f"nepkappa {command}: error: {exc}") from None
        if selected != command:
            validate_input_shape(data, command=selected)
    return selected, _validate_command(args, selected, parser, explicit_mesh=explicit_mesh)


def resolve_force_constant_dimensions(args):
    """Fill omitted FC dimensions from the standard defaults."""
    default_dim = [4, 4, 1]
    if args.dim_fc2 is None:
        args.dim_fc2 = list(default_dim)
    if args.dim_fc3 is None:
        args.dim_fc3 = list(default_dim)
    if args.dim_fc4 is None:
        args.dim_fc4 = list(args.dim_fc3)
