"""Human-readable configuration presentation, independent of parsing."""

from nepkappa.config_validation import qha_sscha_transport_flags

def iter_display_args(args, command=None):
    """Iterate over user-facing config values, hiding inactive route settings."""
    if command is None:
        command = "tdbte" if args.workflow_steps == ["tdbte"] else "run"
    # The plan is derived internally from the six-section input, not a user option.
    compatibility_only = {
        "config_path", "invoked_command", "workflow_steps",
        # Internal route flags are derived from the force-constant switches.
        "qha_volumes", "kappa_enabled",
    }
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
    calculator_fields = {
        "nep_model",
        "nep_force_backend", "nep_batch_size", "nep_command",
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
        "mesh", "temps", "method", "kappa_engine", "kappa_command", "isotope", "bfmp",
        "wigner", "lbte_parallel", "kappa_force_constants",
    }
    command_fields = None
    tdbte_fields = {name for name in vars(args) if name.startswith("tdbte_")}
    if command == "report":
        command_fields = {"result_dir"}
    elif command == "tdbte":
        command_fields = {"result_dir"} | tdbte_fields
    elif command == "run":
        command_fields = (
            common | calculator_fields | fc2_fields | fc3_fields
            | fc4_fields | kappa_fields | qha_only | scph_only
            | qha_sscha_only | fourphonon_only
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
        command_fields = common | scph_only
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
        if "tdbte" in args.workflow_steps:
            command_fields = {"result_dir", "plot_layout", "plot_dpi", "plot_tdbte_times",
                              "plot_path", "plot_path_points", "plot_path_segments"}
        elif args.scph_enabled and not (args.scph_run_transport or args.fp_enabled):
            command_fields = common | {"mesh", "temps", "scph_temps", "plot_layout", "plot_dpi",
                                       "plot_path", "plot_path_points", "plot_path_segments"}
        elif args.qha_enabled and not args.scph_enabled and not args.qha_volumes:
            command_fields = {"result_dir", "plot_layout", "plot_dpi", "plot_path",
                              "plot_path_points", "plot_path_segments", "qha_mesh", "qha_temps"}
        else:
            command_fields = common | {"mesh", "temps"} | {
                name for name in vars(args) if name.startswith("plot_")
            }
    if command == "run" and "tdbte" in args.workflow_steps:
        command_fields = {"result_dir"} | tdbte_fields
        if "fc2fc3" in args.workflow_steps:
            command_fields |= (
                common | calculator_fields | fc2_fields | fc3_fields | hiphive_only
                | {"plot_layout", "plot_dpi"}
            )
    if command in {"run", "qha-sscha"}:
        targets = set(args.workflow_steps) if command == "run" else {command}
        coupled_transport = (
            "qha-sscha" in targets and any(qha_sscha_transport_flags(args))
        )
        transport_active = bool(
            targets & {"kappa", "kappa4", "qha-kappa"}
            or ("scph" in targets and args.scph_run_transport)
            or coupled_transport
        )
        fc3_active = bool(
            "fc2fc3" in targets or "qha-kappa" in targets
            or ("qha-sscha" in targets and (coupled_transport or args.scph_bubble))
        )
        if not transport_active:
            command_fields -= kappa_fields
        if not fc3_active:
            command_fields -= fc3_fields
    for arg, value in vars(args).items():
        if value is None:
            continue
        if args.kappa_command and arg in {"method", "isotope", "bfmp", "wigner", "lbte_parallel"}:
            # The custom command owns these settings; parser defaults are inactive.
            continue
        if command_fields is not None and arg not in command_fields:
            continue
        if arg == "lbte_parallel" and not value:
            continue
        if arg == "force_parallel" and not value:
            continue
        if arg.startswith("scph_bubble") and not args.scph_bubble and command != "bubble":
            continue
        if arg == "scph_save_datasets" and not value:
            continue
        if arg in {"calculator_kwargs", "calculator_model_files"} and not value:
            continue
        if arg in compatibility_only:
            continue
        if arg == "qha_run_transport" and command == "qha":
            continue
        if arg == "scph_run_transport" and args.fp_enabled:
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
        if not args.fp_enabled and arg in fc4_fields and command == "run":
            continue
        if args.calculator == "vasp" and arg == "nep_model" and value is None:
            continue
        if arg in {"nep_force_backend", "nep_batch_size", "nep_command"}:
            if args.calculator != "nep":
                continue
            if arg != "nep_force_backend" and args.nep_force_backend == "cpu":
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
    lines = ["Input configuration:"]
    for arg, value in iter_display_args(args, command=command):
        lines.append(f"  {arg:<15} : {value}")
    return "\n".join(lines)
