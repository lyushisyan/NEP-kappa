"""NEP displacement-force batches evaluated by Calorine's GPU prediction tool."""

from pathlib import Path
import shutil

import numpy as np

from nepkappa.artifacts import ForceArtifactStore


def resolve_nep_executable(command):
    """Require one executable, so its exact binary can be fingerprinted."""
    executable = shutil.which(str(Path(command).expanduser()))
    if executable is None:
        raise FileNotFoundError(
            f"NEP GPU executable not found: {command!r}. Set calculator.nep-command "
            "to the nep executable from GPUMD (not gpumd)."
        )
    return str(Path(executable).resolve())


def predict_forces(structures, *, model, command):
    """Calculate one batch in an isolated temporary directory; never fall back to CPU."""
    from calorine import tools

    predict = getattr(tools, "batch_predict_properties", None)
    if predict is None:
        raise RuntimeError(
            "NEP GPU batches require Calorine >= 3.5 with "
            "calorine.tools.batch_predict_properties."
        )
    results = predict(
        structures, model=str(Path(model).resolve()), command=command
    )
    if len(results) != len(structures):
        raise RuntimeError("NEP GPU result count does not match input structures.")
    forces = []
    for index, (expected, actual) in enumerate(zip(structures, results), 1):
        if not ForceArtifactStore.structures_match(expected, actual) or not np.array_equal(
            expected.pbc, actual.pbc
        ):
            raise RuntimeError(f"NEP GPU result #{index} does not match its input structure.")
        forces.append(ForceArtifactStore.validate_forces(
            actual.get_forces(), len(expected), f"NEP GPU result #{index}"
        ))
    return forces
