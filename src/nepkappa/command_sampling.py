"""Read explicit sampling options without executing a custom phono3py command."""

from __future__ import annotations

import math
import shlex


def phono3py_sampling(command):
    """Return mesh and compact temperatures from common phono3py CLI forms.

    Unknown options are left to phono3py. Do not guess settings hidden in an
    external configuration file or script, or inherit phono3py's defaults.
    """
    aliases = {"--mesh": "mesh", "--ts": "ts", "--temperature": "ts",
               "--tmin": "tmin", "--tmax": "tmax", "--tstep": "tstep"}
    tokens = shlex.split(command)
    options = {}
    index = 0
    while index < len(tokens):
        token = tokens[index]
        index += 1
        if token == "--":
            break
        flag, separator, attached = token.partition("=")
        if flag not in aliases:
            continue
        name = aliases[flag]
        if name in options:
            raise ValueError(f"kappa.command repeats {flag}; specify it once.")
        values = attached.split() if separator else []
        if not separator:
            while index < len(tokens):
                parts = tokens[index].split()
                try:
                    [float(part) for part in parts]
                except ValueError:
                    break
                values.extend(parts)
                index += 1
                if name not in {"mesh", "ts"}:
                    break
        try:
            numbers = [float(value) for value in values]
        except ValueError:
            raise ValueError(f"kappa.command {flag} requires numeric values.") from None
        if not numbers or any(not math.isfinite(value) for value in numbers):
            raise ValueError(f"kappa.command {flag} requires finite numeric values.")
        if name not in {"mesh", "ts"} and len(numbers) != 1:
            raise ValueError(f"kappa.command {flag} requires one value.")
        options[name] = numbers

    result = {}
    if "mesh" in options:
        mesh = options["mesh"]
        if len(mesh) != 3 or any(value <= 0 or not value.is_integer() for value in mesh):
            raise ValueError("Automatic kappa.command mesh parsing requires three positive integers.")
        result["mesh"] = [int(value) for value in mesh]
    range_keys = {"tmin", "tmax", "tstep"} & options.keys()
    if "ts" in options:
        if range_keys:
            raise ValueError("kappa.command must use either --ts or --tmin/--tmax/--tstep.")
        points = options["ts"]
        if any(value <= 0 for value in points):
            raise ValueError("kappa.command --ts requires positive temperatures.")
        if len(points) == 1:
            result["temps"] = points
        else:
            step = points[1] - points[0]
            if step <= 0 or any(not math.isclose(b - a, step, rel_tol=1e-10, abs_tol=1e-10)
                                for a, b in zip(points, points[1:])):
                raise ValueError("Automatic kappa.command --ts parsing requires a single temperature "
                                 "or an increasing, equally spaced list.")
            result["temps"] = [points[0], points[-1], step]
    elif range_keys:
        if len(range_keys) != 3:
            raise ValueError("Specify all of --tmin, --tmax and --tstep in kappa.command.")
        result["temps"] = [options[name][0] for name in ("tmin", "tmax", "tstep")]
    return result
