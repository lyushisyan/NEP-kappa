"""Compile a stage-first YAML plan into the existing scientific stage runner."""

from __future__ import annotations

from collections.abc import Mapping
import math


STAGE_FIELDS = {
    "structure": {"enabled", "method"},
    "forces": {"enabled", "method"},
    "force_constants": {"enabled", "method", "orders", "features"},
    "temperature": {"enabled", "method", "features"},
    "transport": {"enabled", "method", "features"},
    "analysis": {"enabled", "method", "features"},
}


def compile_stage_plan(stages, flat):
    """Resolve stage choices without changing backend-specific implementations."""
    if not isinstance(stages, Mapping) or not stages:
        raise ValueError("workflow.stages must be a non-empty mapping.")
    if "workflow_preset" in flat or "workflow_steps" in flat:
        raise ValueError(
            "workflow.stages cannot be combined with workflow.preset or workflow.steps."
        )
    normalized = {}
    for raw_name, raw_spec in stages.items():
        name = _key(raw_name)
        if name not in STAGE_FIELDS:
            raise ValueError(f"Unknown workflow stage '{raw_name}'.")
        if name in normalized:
            raise ValueError(f"Duplicate workflow stage '{raw_name}'.")
        if not isinstance(raw_spec, Mapping):
            raise ValueError(f"workflow.stages.{raw_name} must be a mapping.")
        spec = {}
        for raw_key, value in raw_spec.items():
            key = _key(raw_key)
            if key not in STAGE_FIELDS[name]:
                raise ValueError(
                    f"Unknown key '{raw_key}' in workflow.stages.{raw_name}."
                )
            if key in spec:
                raise ValueError(f"Duplicate key '{raw_key}' in workflow.stages.{raw_name}.")
            spec[key] = value
        normalized[name] = spec

    steps = []
    _structure(normalized.get("structure"), flat, steps)
    _forces(normalized.get("forces"), flat)
    _force_constants(normalized.get("force_constants"), flat, steps)
    temperature_method = _temperature(normalized.get("temperature"), flat, steps)
    transport_enabled = _transport(normalized.get("transport"), flat, steps)
    if temperature_method in {"qha", "sscha", "qha-sscha"} and transport_enabled:
        raise ValueError(
            "workflow.stages.temperature and transport cannot be chained as a "
            "temperature-corrected conductivity. For QHA-volume SSCHA, enable "
            "transport in temperature.features; otherwise use separate runs."
        )
    _analysis(normalized.get("analysis"), steps)
    if not steps:
        raise ValueError("workflow.stages does not enable any executable stage.")
    flat["workflow_preset"] = "custom"
    flat["workflow_steps"] = steps


def _structure(spec, flat, steps):
    if spec is None or not _enabled(spec, "structure"):
        return
    method = _method(spec, "structure", {"input", "relax"})
    if method == "relax":
        _set(flat, "do_relax", True, "structure.method")
        steps.append("relax")
    else:
        _set(flat, "do_relax", False, "structure.method")


def _forces(spec, flat):
    if spec is None or not _enabled(spec, "forces"):
        return
    method = spec.get("method")
    if not isinstance(method, str) or not method.strip():
        raise ValueError("workflow.stages.forces.method must name a calculator.")
    _set(flat, "calculator", method.strip(), "forces.method")


def _force_constants(spec, flat, steps):
    if spec is None or not _enabled(spec, "force-constants"):
        return
    method = _method(
        spec, "force-constants", {"finite-displacement", "hiphive", "thirdorder"}
    )
    orders = spec.get("orders", [2, 3])
    if not isinstance(orders, list) or orders not in ([2], [2, 3], [2, 3, 4]):
        raise ValueError(
            "workflow.stages.force-constants.orders must be [2], [2, 3], "
            "or [2, 3, 4]."
        )
    if method == "thirdorder" and orders == [2]:
        raise ValueError("Thirdorder requires order 3 in force-constants.orders.")
    _set(flat, "use_hiphive", method == "hiphive", "force-constants.method")
    _set(
        flat,
        "fc3_backend",
        "thirdorder" if method == "thirdorder" else "phono3py",
        "force-constants.method",
    )
    features = _features(spec, "force-constants", {"compact"})
    if "compact" in features:
        _set(
            flat,
            "compact_fc",
            _bool(features["compact"], "force-constants.features.compact"),
            "force-constants.features.compact",
        )
    steps.append("fc2" if orders == [2] else "fc2fc3")
    if 4 in orders:
        export_format = flat.get("fc_format", "both")
        if export_format not in {"both", "shengbte"}:
            raise ValueError(
                "force-constants.orders including 4 requires "
                "force-constant.format: both or shengbte."
            )
        flat["fc_format"] = export_format
        steps.append("fc4")


def _temperature(spec, flat, steps):
    if spec is None or not _enabled(spec, "temperature"):
        return None
    method = _method(spec, "temperature", {"qha", "sscha", "qha-sscha"})
    features = _features(
        spec, "temperature", {"bubble", "three_phonon", "four_phonon"}
    )
    if method == "qha" and features:
        raise ValueError("QHA does not support temperature.features in this stage plan.")
    if method == "sscha" and "four_phonon" in features:
        raise ValueError(
            "Standalone SSCHA does not support four-phonon transport; "
            "use temperature.method: qha-sscha."
        )
    if method in {"qha", "qha-sscha"}:
        if not flat.get("qha_enabled"):
            raise ValueError("temperature.method requires a qha section.")
        steps.append("qha")
    if method in {"sscha", "qha-sscha"}:
        if not flat.get("scph_enabled"):
            raise ValueError("temperature.method requires an scph section.")
        if "bubble" in features:
            _set(
                flat,
                "scph_bubble",
                _bool(features["bubble"], "temperature.features.bubble"),
                "temperature.features.bubble",
            )
        if "three_phonon" in features:
            key = (
                "scph_run_transport"
                if method == "sscha"
                else "qha_sscha_three_phonon"
            )
            _set(
                flat,
                key,
                _bool(
                    features["three_phonon"], "temperature.features.three-phonon"
                ),
                "temperature.features.three-phonon",
            )
        if method == "qha-sscha":
            flat["qha_sscha_enabled"] = True
            if "four_phonon" in features:
                _set(
                    flat,
                    "qha_sscha_four_phonon",
                    _bool(
                        features["four_phonon"],
                        "temperature.features.four-phonon",
                    ),
                    "temperature.features.four-phonon",
                )
            steps.append("qha-sscha")
        else:
            steps.append("scph")
    return method


def _transport(spec, flat, steps):
    if spec is None or not _enabled(spec, "transport"):
        return False
    method = _method(
        spec,
        "transport",
        {
            "three-phonon-rta",
            "three-phonon-lbte",
            "three-phonon-wigner",
            "four-phonon-rta",
            "four-phonon-wigner",
        },
    )
    if method.startswith("three-phonon"):
        features = _features(spec, "transport", {"isotope", "boundary_mfp"})
        _set(
            flat,
            "method",
            "lbte" if method == "three-phonon-lbte" else "rta",
            "transport.method",
        )
        _set(
            flat,
            "wigner",
            method == "three-phonon-wigner",
            "transport.method",
        )
        if "isotope" in features:
            _set(
                flat,
                "isotope",
                _bool(features["isotope"], "transport.features.isotope"),
                "transport.features.isotope",
            )
        if "boundary_mfp" in features:
            value = features["boundary_mfp"]
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
                or value <= 0
            ):
                raise ValueError("transport.features.boundary-mfp must be positive.")
            _set(flat, "bfmp", value, "transport.features.boundary-mfp")
        steps.append("kappa")
    else:
        features = _features(spec, "transport", {"isotope", "nonanalytic"})
        if not flat.get("fp_enabled"):
            raise ValueError(
                "Four-phonon transport requires a fourphonon section with "
                "the external executable settings."
            )
        _set(flat, "fp_solver", "rta", "transport.method")
        _set(
            flat,
            "fp_wigner",
            method == "four-phonon-wigner",
            "transport.method",
        )
        if "isotope" in features:
            _set(
                flat,
                "fp_isotopes",
                _bool(features["isotope"], "transport.features.isotope"),
                "transport.features.isotope",
            )
        if "nonanalytic" in features:
            _set(
                flat,
                "fp_nonanalytic",
                _bool(features["nonanalytic"], "transport.features.nonanalytic"),
                "transport.features.nonanalytic",
            )
        steps.append("kappa4")
    return True


def _analysis(spec, steps):
    if spec is None or not _enabled(spec, "analysis"):
        return
    selected = set()
    if "method" in spec:
        selected.add(_method(spec, "analysis", {"tdbte", "plot", "report"}))
    features = _features(spec, "analysis", {"tdbte", "plot", "report"})
    if selected and features.get(next(iter(selected))) is False:
        raise ValueError(
            "workflow.stages.analysis.method conflicts with a disabled "
            "analysis feature of the same name."
        )
    for feature in ("tdbte", "plot", "report"):
        if feature in features and _bool(
            features[feature], f"analysis.features.{feature}"
        ):
            selected.add(feature)
    for feature in ("tdbte", "plot", "report"):
        if feature in selected:
            steps.append(feature)


def _enabled(spec, stage):
    value = spec.get("enabled", True)
    enabled = _bool(value, f"{stage}.enabled")
    if not enabled and any(key != "enabled" for key in spec):
        raise ValueError(
            f"workflow.stages.{stage} is disabled but also configures a method or features."
        )
    return enabled


def _method(spec, stage, allowed):
    value = spec.get("method")
    if not isinstance(value, str):
        raise ValueError(f"workflow.stages.{stage}.method is required.")
    method = _key(value).replace("_", "-")
    if method not in allowed:
        raise ValueError(
            f"Unsupported workflow.stages.{stage}.method '{value}'. "
            f"Choose: {', '.join(sorted(allowed))}."
        )
    return method


def _features(spec, stage, allowed):
    raw = spec.get("features", {})
    if not isinstance(raw, Mapping):
        raise ValueError(f"workflow.stages.{stage}.features must be a mapping.")
    features = {}
    for raw_key, value in raw.items():
        key = _key(raw_key)
        if key not in allowed:
            raise ValueError(
                f"Unsupported feature '{raw_key}' in workflow.stages.{stage}."
            )
        if key in features:
            raise ValueError(
                f"Duplicate feature '{raw_key}' in workflow.stages.{stage}."
            )
        features[key] = value
    return features


def _bool(value, location):
    if not isinstance(value, bool):
        raise ValueError(f"workflow.stages.{location} must be true or false.")
    return value


def _set(flat, key, value, source):
    if key in flat and flat[key] != value:
        raise ValueError(
            f"workflow.stages.{source} conflicts with the detailed "
            f"configuration for '{key}'."
        )
    flat[key] = value


def _key(value):
    return str(value).strip().lower().replace("-", "_")
