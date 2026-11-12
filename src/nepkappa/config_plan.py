"""Derive scientific stages from input switches without executing calculations."""

from nepkappa.config_schema import normalize_yaml_key, raise_unknown_yaml_key, str2bool, yaml_input_sections
from nepkappa.command_registry import calculation_steps


THERMAL_SECTIONS = {
    "structure", "calculator", "force_constant", "kappa", "plot", "output"
}


DYNAMIC_SECTIONS = (THERMAL_SECTIONS - {"kappa"}) | {"tdbte"}


STAGE_REQUIRED_SECTIONS = {
    "relax": {"structure", "calculator"},
    "fc2": {"structure", "calculator", "force_constant"},
    "fc2fc3": {"structure", "calculator", "force_constant"},
    "fc2fc3fc4": {"structure", "calculator", "force_constant"},
    "fc4": {"structure", "calculator", "force_constant"},
    "qha": {"structure", "calculator", "force_constant"},
    "scph": {"structure", "calculator", "force_constant"},
    "qha-kappa": {"structure", "calculator", "force_constant", "kappa"},
    "qha-sscha": {"structure", "calculator", "force_constant", "kappa"},
    "kappa": {"kappa"},
    "kappa4": {"kappa", "force_constant"},
    "tdbte": {"tdbte"},
    "plot": set(),
    "report": {"output"},
    "bubble": {"force_constant"},
}


WORKFLOW_STEPS = {
    "relax", "fc2", "fc2fc3", "fc4", "qha", "qha-kappa", "scph", "qha-sscha",
    "kappa", "kappa4", "plot", "report", "bubble", "tdbte",
}



def expand_input_modes(data, command=None):
    """Derive ordered stages for thermal, full TD-BTE, or artifact reuse inputs."""
    data = dict(data)
    # Defaults are internal: users need not write unrelated empty sections.
    keys = {normalize_yaml_key(key): key for key in data}
    explicit_kappa = "kappa" in keys
    if "tdbte" not in keys:
        for name in THERMAL_SECTIONS - keys.keys():
            data[name.replace("_", "-")] = {}
    data = _expand_static_parallel(data)
    keys = {normalize_yaml_key(key): key for key in data}
    if "tdbte" in keys and "workflow" not in keys:
        plot = data.get(keys.get("plot"), {})
        unsupported_plot = {normalize_yaml_key(key) for key in plot} - {
            "layout", "dpi", "tdbte_times", "path", "path_points", "path_segments"
        }
        if unsupported_plot:
            raise ValueError(
                "TD-BTE plot supports layout, dpi, tdbte-times and band-path settings; remove plot."
                + sorted(unsupported_plot)[0].replace("_", "-") + "."
            )
        if command not in {None, "run", "info"}:
            tdbte_key = keys["tdbte"]
            tdbte = dict(data[tdbte_key])
            if not {"kernel", "force_constants"} & {normalize_yaml_key(key) for key in tdbte}:
                output = data.get(keys.get("output"), {})
                tdbte["force-constants"] = next(
                    (value for key, value in output.items()
                     if normalize_yaml_key(key) == "result_dir"), "result"
                )
            # A direct command never prepares new force constants.
            return {**data, tdbte_key: tdbte, "workflow": {"steps": ["tdbte"]}}
        fc = data[keys["force_constant"]]
        forbidden = {"qha", "sscha", "four_phonon"} & {
            normalize_yaml_key(key) for key in fc
        }
        if forbidden:
            raise ValueError(
                "TD-BTE currently uses FC2/FC3 only; remove force-constant."
                + sorted(forbidden)[0].replace("_", "-") + "."
            )
        backend = next(
            (value for key, value in fc.items()
             if normalize_yaml_key(key) == "fc3_backend"), "phono3py"
        )
        if backend != "phono3py":
            raise ValueError("Six-section TD-BTE requires force-constant.fc3-backend: phono3py.")
        fc_format = next(
            (value for key, value in fc.items()
             if normalize_yaml_key(key) == "format"), "phono3py"
        )
        if fc_format not in {"phono3py", "both"}:
            raise ValueError("Six-section TD-BTE requires force-constant.format: phono3py or both.")
        tdbte_key = keys["tdbte"]
        tdbte = data[tdbte_key]
        if any(normalize_yaml_key(key) in {"force_constants", "kernel"} for key in tdbte):
            raise ValueError(
                "Six-section TD-BTE generates FC2/FC3; omit tdbte.force-constants "
                "and tdbte.kernel. Use the two-section input to reuse artifacts."
            )
        output = data[keys["output"]]
        result_dir = next(
            (value for key, value in output.items()
             if normalize_yaml_key(key) == "result_dir"), "result"
        )
        if not isinstance(result_dir, str) or not result_dir.strip():
            raise ValueError("output.result-dir must be a non-empty path.")
        structure = data[keys["structure"]]
        relaxation = next(
            (value for key, value in structure.items()
             if normalize_yaml_key(key) == "relaxation"), False
        )
        relax_requested = relaxation is True or (
            isinstance(relaxation, dict) and relaxation.get("enabled", True) is True
        )
        return {
            **data,
            tdbte_key: {**tdbte, "force-constants": result_dir},
            "workflow": {
                "steps": (["relax"] if relax_requested else [])
                + ["fc2fc3", "tdbte"]
            },
        }

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
    for obsolete in ("enabled", "qha_volumes"):
        if obsolete in kappa_values:
            raise ValueError(
                f"kappa.{obsolete.replace('_', '-')} is no longer accepted; "
                "set run-transport inside force-constant.qha or "
                "force-constant.sscha."
            )
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
        return {
            **data,
            keys["kappa"]: {**kappa, "enabled": True},
            "workflow": {"steps": (
                list(calculation_steps(command)) if command in STAGE_REQUIRED_SECTIONS
                and not (command == "plot" and explicit_kappa)
                else ["relax", "fc2fc3", "kappa"]
            )},
        }
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
    if "temps" in four_options and four_options["temps"] != kappa_values.get("temps"):
        raise ValueError(
            "force-constant.four-phonon.temps must match kappa.temps; "
            "set transport temperatures in kappa.temps."
        )
    selected_qha_transport = qha_options.get("run_transport")
    selected_sscha_transport = sscha_options.get("run_transport")
    for section, selected in (("qha", selected_qha_transport), ("sscha", selected_sscha_transport)):
        if selected is not None and not isinstance(selected, bool):
            raise ValueError(f"force-constant.{section}.run-transport must be true or false.")
    if four_on and qha_on and not sscha_on:
        raise ValueError("QHA plus four-phonon transport requires SSCHA.")
    if four_on and sscha_on and not qha_on:
        raise ValueError("Standalone SSCHA plus four-phonon transport is not supported.")
    if qha_on and sscha_on:
        if selected_qha_transport is not None and selected_sscha_transport is not None:
            if selected_qha_transport != selected_sscha_transport:
                raise ValueError(
                    "QHA+SSCHA transport switches disagree; QHA volumes and SSCHA FC2 "
                    "must be used together."
                )
        if four_on and (selected_qha_transport is False or selected_sscha_transport is False):
            raise ValueError("Four-phonon transport cannot disable QHA or SSCHA use.")
        qha_transport = sscha_transport = bool(
            four_on or selected_qha_transport or selected_sscha_transport
        )
    else:
        qha_transport = bool(selected_qha_transport)
        sscha_transport = bool(selected_sscha_transport)
    if qha_on:
        qha_options["run_transport"] = qha_transport
    kappa_enabled = four_on or qha_transport or sscha_transport or not (qha_on or sscha_on)
    transport_requested = command in {
        None, "run", "info", "kappa", "kappa4", "qha-kappa", "qha-sscha", "scph"
    }
    if transport_requested and kappa_enabled and (qha_on or sscha_on) and "temps" not in kappa_values:
        raise ValueError(
            "QHA/SSCHA transport requires explicit kappa.temps; "
            "choose temperatures within the QHA range and, when SSCHA is enabled, "
            "on its sampling grid."
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
    qha_volumes = qha_on and not sscha_on and qha_transport
    if qha_volumes and kappa_values.get("method", "rta") != "rta":
        raise ValueError("QHA-volume transport currently requires kappa.method: rta.")
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
    result[keys["kappa"]] = {
        key: value for key, value in kappa.items()
        if normalize_yaml_key(key) not in separate_methods
    }
    result[keys["kappa"]]["enabled"] = kappa_enabled
    result[keys["kappa"]]["qha-volumes"] = qha_volumes
    if qha_on:
        result["qha"] = {"enabled": True, **qha_options}
    if sscha_on:
        sscha_options["run_transport"] = sscha_transport and not four_on
        result["scph"] = {"enabled": True, **sscha_options}
    if four_on:
        result["fourphonon"] = {"enabled": True, **four_options}
    if qha_on and sscha_on:
        result["qha-sscha"] = {
            "three-phonon": sscha_options["run_transport"],
            "four-phonon": four_on and kappa_enabled,
        }
        steps = ["qha", "qha-sscha"]
    elif qha_on:
        steps = ["qha", "qha-kappa"] if kappa_enabled else ["qha"]
    elif sscha_on:
        needs_fc3 = sscha_transport or str2bool(sscha_options.get("bubble", False))
        force_step = "fc2fc3" if needs_fc3 else "fc2"
        steps = (["relax"] if relax_requested else []) + [force_step, "scph"]
    elif four_on:
        steps = ["relax", "fc2fc3", "fc4", "kappa4"]
    else:
        steps = ["relax", "fc2fc3", "kappa"]
    result["workflow"] = {"steps": steps}
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
    dynamic = "tdbte" in keys
    if fc_key is None or (kappa_key is None and not dynamic):
        raise ValueError("parallel requires force-constant and a calculation section.")
    fc = dict(result[fc_key])
    kappa = dict(result[kappa_key]) if kappa_key is not None else None
    for raw_target, settings in parallel.items():
        target = normalize_yaml_key(raw_target)
        if dynamic and target != "force_constant":
            raise ValueError("Six-section TD-BTE supports only parallel.force-constant.")
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
    if kappa_key is not None:
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



def resolve_workflow_plan(args, parser=None):
    """Validate the stages derived from static switches or the dynamic input."""
    steps = list(args.workflow_steps or [])
    error = None
    if not steps:
        error = "Could not derive run stages from the input."
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
