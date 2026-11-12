"""Explicit force-cache dependencies, independent of fitting and transport.

Installed force providers must remain immutable for a process lifetime. Their
source and native binaries are hashed once per process; model inputs are always
hashed again. Opaque external backends calculate normally but do not reuse a
persistent cache that this dependency list cannot certify.
"""

import ast
from functools import lru_cache
from importlib import metadata
from pathlib import Path
import shutil
import shlex
import os
from uuid import uuid4

from ase.io import read

from nepkappa.identity import canonical_data, data_sha256, file_identity, file_sha256
from nepkappa.run_state import RunStateStore


def reuse_scope(core):
    """Limit uncertified backends to one calculation or Slurm pipeline."""
    scope = getattr(core, '_force_reuse_scope', None)
    if scope is None:
        state = {}
        if os.environ.get('NEPKAPPA_FORCE_COLLECT') == '1':
            state = RunStateStore(core.output_dir / 'force-slurm' / 'submission.yaml').read()
        scope = state.get('reuse_scope') or uuid4().hex
        core._force_reuse_scope = scope
    return scope


# Keep this explicit list in sync with force evaluation, serialization and parsing.
WORKFLOW_FORCE_METHODS = {
    '__init__', '_calculator_name', '_make_nep_calculator', '_make_external_backend',
    '_make_vasp_backend', '_make_force_store', '_calculate_forces',
    '_calculate_forces_cached', '_calculate_nep_batch_cached', '_force_job_dir',
    '_run_command', '_run_command_with_input', '_write_vasp_force_xml',
    '_run_nep_thirdorder_job', '_run_nep_fc4_job', '_run_external_order_job',
    '_run_vasp_fc4_job', '_run_vasp_thirdorder_job',
    '_pattern_digest', '_job_input_payload', '_job_input_fingerprint',
    '_job_matches_inputs', '_record_job_inputs', '_stage_force_job',
}


def source_digest(path, methods=None):
    """Hash executable syntax, excluding comments, docstrings and locations."""
    tree = ast.parse(Path(path).read_text(encoding='utf-8'))
    if methods is not None:
        context = next(node for node in tree.body if isinstance(node, ast.ClassDef)
                       and node.name == 'NEPPhononWorkflow')
        selected = [node for node in context.body if getattr(node, 'name', None) in methods]
        missing = set(methods) - {node.name for node in selected}
        if missing:
            raise ValueError(f'Missing force-identity methods: {sorted(missing)}')
        # Imports define the global names used by the selected methods.
        tree.body = selected + [node for node in tree.body
                                if isinstance(node, (ast.Import, ast.ImportFrom))]
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            if (node.body and isinstance(node.body[0], ast.Expr)
                    and isinstance(node.body[0].value, ast.Constant)
                    and isinstance(node.body[0].value.value, str)):
                node.body.pop(0)
    return data_sha256(ast.dump(tree, include_attributes=False))


@lru_cache(maxsize=None)
def force_implementation(calculator):
    root = Path(__file__).parent
    modules = ['identity.py', 'force_identity.py', 'artifacts.py', 'force_worker.py']
    if calculator == 'vasp':
        modules += ['adapters/vasp.py']
    elif calculator == 'nep':
        modules += ['nep_batch.py']
    else:
        modules += ['calculators.py']
    return {**{name: source_digest(root / name) for name in modules},
            'workflow.force_methods': source_digest(root / 'workflow.py', WORKFLOW_FORCE_METHODS)}


@lru_cache(maxsize=None)
def provider_identity(name):
    """Fingerprint installed Python/native provider files, including bundled libs."""
    try:
        distribution = metadata.distribution(name)
    except metadata.PackageNotFoundError:
        return {'version': None, 'implementation': None}
    files = {}
    for entry in distribution.files or ():
        if entry.suffix in {'.py', '.so', '.pyd', '.dll', '.dylib'} or '.so.' in entry.name:
            path = Path(distribution.locate_file(entry))
            if path.is_file():
                files[str(entry)] = file_sha256(path)
    return {'version': distribution.version,
            'implementation': data_sha256(files) if files else None}


def force_engine(calculator):
    names = ['numpy', 'ase']
    if calculator == 'nep':
        names.append('calorine')
    elif calculator not in {'vasp'}:
        names.append('scipy')
    return {name: provider_identity(name) for name in names}


def reusable_backend(core, payload):
    """Whether the explicit dependency contract supports cross-run reuse."""
    if not all(item['implementation'] for item in payload['software'].values()):
        return False
    calculator = core._calculator_name()
    if calculator == 'nep':
        model = payload['calculator_inputs']['model']
        return bool(model and model['exists'])
    if calculator == 'vasp':
        command = payload['calculator_inputs']['command']
        executable = command['executable']
        if not executable or len(command['files']) != 1:
            return False  # Launcher/wrapper dependency closure is not certified.
        with Path(executable['path']).open('rb') as handle:
            magic = handle.read(4)
        return magic in {b'\x7fELF', b'\xcf\xfa\xed\xfe', b'\xfe\xed\xfa\xcf', b'\xca\xfe\xba\xbe'}
    factory = getattr(core._make_external_backend(), 'factory_metadata', {})
    module = factory.get('implementation', {}).get('module')
    return module in {'ase.calculators.lj', 'ase.calculators.emt', 'ase.calculators.morse'}


def command_identity(command):
    """Describe a configured external command and its executable, if resolvable."""
    tokens = list(command)
    executable = shutil.which(tokens[0]) if tokens else None
    if tokens and Path(tokens[0]).expanduser().is_file():
        executable = Path(tokens[0]).expanduser()
    files = {}
    for token in tokens:
        candidate = Path(token).expanduser()
        resolved = candidate if candidate.is_file() else shutil.which(token)
        if resolved:
            files[str(token)] = file_identity(resolved)
    return {
        "tokens": tokens,
        "executable": file_identity(executable) if executable else None,
        "files": files,
    }


def potcar_inputs(core, atoms):
    """Return identities of the exact POTCAR sources selected for a VASP job."""
    configured = getattr(core.cfg, "potcar_path", None)
    symbols = list(dict.fromkeys(atoms.get_chemical_symbols()))
    sources = []
    if configured:
        path = Path(configured).expanduser()
        if path.is_file():
            sources = [path]
        elif path.is_dir():
            sources = [core._make_vasp_backend().find_potcar(path, symbol) for symbol in symbols]
    else:
        detected = core._make_vasp_backend().detect_potcar(symbols)
        if detected is not None:
            _, detected_sources = detected
            sources = list(detected_sources)
    return [
        file_identity(source) if source is not None else {"exists": False}
        for source in sources
    ] or [{"path": str(configured) if configured else None, "exists": False}]


def force_input_payload(core, pattern, backend):
    """Build all scientifically relevant inputs used to calculate one force job."""
    calculator = core._calculator_name()
    payload = {
        "schema_version": 3,
        "code": {
            "force_implementation": force_implementation(calculator),
        },
        "backend": backend,
        "displacement_sha256": core._pattern_digest(pattern),
        "calculator": calculator,
    }
    payload["software"] = force_engine(calculator)
    if backend == "thirdorder":
        payload["order_settings"] = {
            "command": command_identity(core._thirdorder_command()),
            "dimension": list(core.cfg.dim_fc3),
            "cutoff": float(core.cfg.cutoff_fc3),
        }
    elif backend == "fourthorder":
        payload["order_settings"] = {
            "command": command_identity(core._fourthorder_command()),
            "dimension": list(core.cfg.dim_fc4),
            "cutoff": float(core.cfg.cutoff_fc4),
        }

    if calculator == "nep":
        model = getattr(core.cfg, "nep_model", None)
        payload["calculator_inputs"] = {
            "model": file_identity(model) if model else None,
        }
        if getattr(core.cfg, "nep_force_backend", "cpu") == "gpu-batch":
            from nepkappa import nep_batch

            if core._nep_batch_identity is None:
                executable = nep_batch.resolve_nep_executable(core.cfg.nep_command)
                core._nep_batch_identity = {
                    "executable": file_identity(executable),
                }
            payload["calculator_inputs"].update({
                "force_backend": "gpu-batch",
                "batch_size": core.cfg.nep_batch_size,
                **core._nep_batch_identity,
            })
    elif calculator == "vasp":
        command = getattr(core.cfg, "vasp_command", None)
        if command:
            command = shlex.split(str(command))
        elif getattr(core.cfg, "vasp_path", None):
            command = [str(core.cfg.vasp_path)]
        else:
            detected = core._make_vasp_backend().detect_path()
            command = [str(detected)] if detected else ["vasp_std"]
        atoms = read(str(pattern), format="vasp")
        payload["calculator_inputs"] = {
            "command": command_identity(command),
            "incar": canonical_data(core._make_vasp_backend().combined_kwargs()),
            "potcar_sources": potcar_inputs(core, atoms),
        }
    else:
        payload["calculator_inputs"] = (
            core._make_external_backend().cache_inputs()
        )
    payload["persistent_reuse"] = reusable_backend(core, payload)
    if not payload["persistent_reuse"]:
        payload["reuse_scope"] = reuse_scope(core)
    return canonical_data(payload)
