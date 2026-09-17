"""Command-line interface for NEP-kappa."""

from __future__ import annotations

import argparse
import contextlib
import os
from pathlib import Path
import sys
import time

from nepkappa import __version__
from nepkappa.command_registry import (
    COMMAND_SPECS,
    EXECUTION_COMMANDS,
    VALIDATION_TARGETS,
    canonical_command,
)
from nepkappa.config import (
    format_compare_config,
    format_config,
    parse_compare_args,
    parse_workflow_args,
)
from nepkappa.provenance import ProvenanceRecorder


class Tee:
    """Write text to terminal and run.log, compacting dynamic terminal updates."""

    def __init__(self, terminal_stream, log_stream):
        self.terminal_stream = terminal_stream
        self.log_stream = log_stream
        self._rewrite_buffer = ""
        self._in_rewrite = False

    def write(self, data):
        self.terminal_stream.write(data)
        self.terminal_stream.flush()
        self._write_log(data)
        self.log_stream.flush()

    def flush(self):
        self.terminal_stream.flush()
        self.log_stream.flush()

    def isatty(self):
        return self.terminal_stream.isatty()

    def _write_log(self, data):
        for char in data:
            if char == "\r":
                self._rewrite_buffer = ""
                self._in_rewrite = True
                continue
            if self._in_rewrite:
                if char == "\n":
                    line = self._rewrite_buffer.rstrip()
                    if line:
                        self.log_stream.write(line + "\n")
                    self._rewrite_buffer = ""
                    self._in_rewrite = False
                else:
                    self._rewrite_buffer += char
                continue
            self.log_stream.write(char)


def main(argv: list[str] | None = None) -> int:
    """Run the NEP-kappa command-line interface."""
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command in EXECUTION_COMMANDS:
        return run_command(args.command, args.config, debug=args.debug)
    if args.command == "compare":
        return compare_command(args.config, debug=args.debug)
    if args.command == "converge":
        return convergence_command(args.config, debug=args.debug)
    if args.command == "status":
        return status_command(args.config)
    if args.command == "report":
        return report_command(args.config)
    if args.command == "init":
        return init_command(args)
    if args.command == "info":
        return info_command(args.config, command=args.workflow)
    if args.command == "validate":
        return validate_command(args.config, command=args.workflow)

    parser.error("missing command")
    return 2


def build_parser() -> argparse.ArgumentParser:
    """Build the top-level parser."""
    parser = argparse.ArgumentParser(
        prog="nepkappa",
        description="NEP-assisted lattice thermal conductivity workflow.",
    )
    parser.add_argument("--version", action="version", version=f"nepkappa {__version__}")
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Show Python tracebacks when a command fails.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    for spec in COMMAND_SPECS:
        command_parser = subparsers.add_parser(spec.name, help=spec.help)
        if spec.name == "init":
            command_parser.add_argument(
                "config",
                nargs="?",
                default="nepkappa.yaml",
                help="YAML file to create (default: nepkappa.yaml)",
            )
            command_parser.add_argument(
                "--preset",
                choices=(
                    "three-phonon",
                    "four-phonon",
                    "qha",
                    "scph",
                    "qha-sscha",
                ),
            )
            command_parser.add_argument(
                "--calculator", choices=("nep", "mace", "vasp")
            )
            command_parser.add_argument("--structure", help="POSCAR/structure path")
            command_parser.add_argument(
                "--model", help="NEP/MACE model or MACE foundation"
            )
            command_parser.add_argument("--result-dir")
            command_parser.add_argument(
                "--dim", nargs=3, type=int, metavar=("NX", "NY", "NZ")
            )
            command_parser.add_argument(
                "--mesh", nargs=3, type=int, metavar=("QX", "QY", "QZ")
            )
            command_parser.add_argument("--vasp-command")
            command_parser.add_argument("--potcar-path")
            command_parser.add_argument("--fourphonon-command")
            command_parser.add_argument("--slurm", action="store_true")
            command_parser.add_argument("--non-interactive", action="store_true")
            command_parser.add_argument("--force", action="store_true")
            continue
        config_help = {
            "compare": "YAML comparison input file",
            "convergence": "YAML convergence-study input file",
            "report": "Result directory or workflow YAML input file",
        }.get(spec.config_kind, "YAML input file")
        command_parser.add_argument("config", help=config_help)
        if spec.name in {"info", "validate"}:
            command_parser.add_argument(
                "--for",
                dest="workflow",
                choices=VALIDATION_TARGETS,
                default=None if spec.name == "info" else "run",
                help=(
                    "Show only settings for this workflow (default: all)."
                    if spec.name == "info"
                    else "Target workflow whose settings are checked (default: run)."
                ),
            )

    return parser


def run_command(command, config_path, *, debug=False) -> int:
    """Run one workflow command."""
    start_time = time.time()
    requested_command = command
    command = canonical_command(command)
    args = parse_workflow_args(config_path, command=command)
    args.config_path = str(Path(config_path).resolve())
    args.invoked_command = command
    os.makedirs(args.result_dir, exist_ok=True)
    log_path = os.path.join(args.result_dir, "run.log")
    log_mode = "w" if command == "run" else "a"
    recorder = ProvenanceRecorder(args.result_dir, requested_command, config_path, args)

    with open(log_path, log_mode, encoding="utf-8") as log_file:
        stdout = Tee(sys.stdout, log_file)
        stderr = Tee(sys.stderr, log_file)
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            if log_mode == "a" and os.path.getsize(log_path) > 0:
                print("\n" + "=" * 60)
            print(f"[main] Command: {requested_command}")
            if requested_command == "fc":
                print("[main] `fc` is deprecated; use `fc2fc3` instead.")
            if config_path is not None:
                print(f"[main] Reading arguments from {config_path}...")
            print(f"[main] Logging output to {log_path}")
            print("-" * 60)
            print(format_config(args, command=command))
            print("-" * 60)

            exit_code = 0
            execution = None
            failure = None
            try:
                manifest_path = recorder.start()
                print(f"[main] Provenance manifest: {manifest_path}")
                from nepkappa.application import execute_workflow_command

                execution = execute_workflow_command(command, args)
            except Exception as exc:
                exit_code = 1
                failure = exc
                print(f"\n[Error] Workflow execution failed: {exc}")
                if debug:
                    import traceback

                    traceback.print_exc()
                else:
                    print("[Hint] Re-run with 'nepkappa --debug ...' for a traceback.")
            finally:
                try:
                    recorder.finish(
                        "complete" if exit_code == 0 else "failed",
                        error=failure,
                        stage_timings=getattr(execution, "stage_timings", None),
                    )
                except Exception as exc:
                    exit_code = 1
                    print(f"\n[Error] Could not finalize provenance manifest: {exc}")
                print("-" * 60)
                print(f"Total Execution Time: {format_duration(time.time() - start_time)}")
                print("-" * 60)

    print("-" * 60)
    print(f"Run log saved to: {log_path}")
    return exit_code


def info_command(config_path, command=None) -> int:
    """Print parsed config values without running the workflow."""
    args = parse_workflow_args(config_path, command=command)
    print(format_config(args, command=command))
    return 0


def validate_command(config_path, command="run") -> int:
    """Validate one workflow configuration without creating result files."""
    parse_workflow_args(config_path, command=command)
    print(f"Configuration is valid for '{command}': {Path(config_path).resolve()}")
    return 0


def status_command(config_path) -> int:
    """Show persistent run state and live Slurm queue information."""
    args = parse_workflow_args(config_path, command="status")
    from nepkappa.run_state import discover_run_states, format_run_states
    from nepkappa.scheduler import SlurmScheduler

    records = discover_run_states(
        args.result_dir,
        scheduler=SlurmScheduler(),
    )
    print(format_run_states(records))
    return 0


def report_command(source) -> int:
    """Generate human- and machine-readable summaries of existing results."""
    from nepkappa.report import generate_report, resolve_result_directory

    try:
        result_dir = resolve_result_directory(source)
        outputs = generate_report(result_dir)
    except (FileNotFoundError, ValueError) as exc:
        print(f"[Error] Could not generate report: {exc}")
        return 1
    print(f"Report generated for: {result_dir}")
    for output in outputs:
        print(f"  - {output}")
    return 0


def init_command(args) -> int:
    """Create a validated starter YAML through the input initializer."""
    from nepkappa.initializer import initialize_input

    try:
        path, answers = initialize_input(
            args.config,
            preset=args.preset,
            calculator=args.calculator,
            structure=args.structure,
            model=args.model,
            result_dir=args.result_dir,
            dimension=args.dim,
            mesh=args.mesh,
            use_slurm=args.slurm,
            vasp_command=args.vasp_command,
            potcar_path=args.potcar_path,
            fourphonon_command=args.fourphonon_command,
            interactive=not args.non_interactive,
            overwrite=args.force,
        )
    except (FileExistsError, ValueError) as exc:
        print(f"[Error] Could not create input: {exc}")
        return 1
    print(f"Created validated {answers.preset} input: {path.resolve()}")
    print(f"Next: nepkappa run {path}")
    return 0


def compare_command(config_path, *, debug=False) -> int:
    """Run multi-model comparison plotting."""
    start_time = time.time()
    args = parse_compare_args(config_path)
    os.makedirs(args.compare_dir, exist_ok=True)
    log_path = os.path.join(args.compare_dir, "compare.log")
    recorder = ProvenanceRecorder(args.compare_dir, "compare", config_path, args)

    with open(log_path, "w", encoding="utf-8") as log_file:
        stdout = Tee(sys.stdout, log_file)
        stderr = Tee(sys.stderr, log_file)
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            print("[main] Command: compare")
            print(f"[main] Reading arguments from {config_path}...")
            print(f"[main] Logging output to {log_path}")
            print("-" * 60)
            print(format_compare_config(args))
            print("-" * 60)

            exit_code = 0
            failure = None
            try:
                manifest_path = recorder.start()
                print(f"[main] Provenance manifest: {manifest_path}")
                from nepkappa.plot import compare_results

                compare_results(args)
            except Exception as exc:
                exit_code = 1
                failure = exc
                print(f"\n[Error] Comparison failed: {exc}")
                if debug:
                    import traceback

                    traceback.print_exc()
                else:
                    print("[Hint] Re-run with 'nepkappa --debug ...' for a traceback.")
            finally:
                try:
                    recorder.finish(
                        "complete" if exit_code == 0 else "failed", error=failure
                    )
                except Exception as exc:
                    exit_code = 1
                    print(f"\n[Error] Could not finalize provenance manifest: {exc}")
                print("-" * 60)
                print(f"Total Execution Time: {format_duration(time.time() - start_time)}")
                print("-" * 60)

    print("-" * 60)
    print(f"Compare log saved to: {log_path}")
    return exit_code


def convergence_command(config_path, *, debug=False) -> int:
    """Prepare, optionally execute, and analyze a convergence study."""
    from nepkappa.convergence import parse_convergence_args, run_convergence_study

    try:
        config = parse_convergence_args(config_path)
        config.directory.mkdir(parents=True, exist_ok=True)
        summary = run_convergence_study(
            config,
            run_case=lambda path: run_command("run", path, debug=debug),
        )
        return 1 if summary["failed_cases"] else 0
    except Exception as exc:
        print(f"[Error] Convergence study failed: {exc}")
        if debug:
            import traceback

            traceback.print_exc()
        else:
            print("[Hint] Re-run with 'nepkappa --debug ...' for a traceback.")
        return 1


def format_duration(seconds):
    """Format elapsed seconds as a compact human-readable duration."""
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    secs = seconds % 60
    if hours:
        return f"{hours}h {minutes}m {secs:.2f}s"
    if minutes:
        return f"{minutes}m {secs:.2f}s"
    return f"{secs:.2f}s"
