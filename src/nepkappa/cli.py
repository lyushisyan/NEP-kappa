"""Command-line interface for NEP-kappa."""

from __future__ import annotations

import argparse
import contextlib
import os
from pathlib import Path
import sys
import time
import warnings

from nepkappa import __version__
from nepkappa.command_registry import (
    COMMAND_SPECS,
    EXECUTION_COMMANDS,
    PUBLIC_COMMANDS,
    VALIDATION_TARGETS,
    canonical_command,
    resolve_calculation_command,
)
from nepkappa.config import (
    format_config,
    parse_workflow_args,
)
from nepkappa.provenance import ProvenanceRecorder
from nepkappa.output import print_detail


class Tee:
    """Write text to terminal and run.log, compacting dynamic terminal updates."""

    def __init__(self, terminal_stream, log_stream, *, verbose=True):
        self.terminal_stream = terminal_stream
        self.log_stream = log_stream
        self.detail_stream = self if verbose else log_stream
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
        return run_command(
            args.command, args.config, debug=args.debug, verbose=args.verbose
        )
    if args.command == "report":
        return report_command(args.config)
    if args.command == "info":
        return info_command(args.config, command=args.workflow)

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
        "-v", "--verbose", action="store_true",
        help="Show detailed output on screen as well as in run.log.",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Show Python tracebacks when a command fails.",
    )
    subparsers = parser.add_subparsers(
        dest="command", required=True,
        metavar="{" + ",".join(PUBLIC_COMMANDS) + "}",
    )

    for spec in COMMAND_SPECS:
        command_parser = subparsers.add_parser(spec.name, help=spec.help)
        config_help = (
            "Result directory or workflow YAML input file"
            if spec.name == "report" else "YAML input file"
        )
        command_parser.add_argument("config", help=config_help)
        command_parser.add_argument(
            "-v", "--verbose", action="store_true", default=argparse.SUPPRESS,
            help="Show detailed output on screen as well as in run.log.",
        )
        if spec.name == "info":
            command_parser.add_argument(
                "--for",
                dest="workflow",
                choices=VALIDATION_TARGETS,
                default="run",
                help="Check and show this workflow's settings (default: run).",
            )

    return parser


@contextlib.contextmanager
def warning_output(stderr, *, verbose):
    """Keep full warning diagnostics in the log and readable messages on screen."""
    with warnings.catch_warnings():
        original_showwarning = warnings.showwarning

        def showwarning(message, category, filename, lineno, file=None, line=None):
            if file is not None:
                return original_showwarning(
                    message, category, filename, lineno, file=file, line=line
                )
            detail = warnings.formatwarning(message, category, filename, lineno, line)
            # Match category names without importing version-specific Phonopy
            # classes: older supported releases do not define both warnings.
            if category.__module__.startswith("phonopy.") and category.__name__ in {
                "PrimitiveMatrixAutoDefaultWarning", "GammaAcousticWarning",
            }:
                stderr.log_stream.write(detail)
                stderr.log_stream.flush()
                return
            stderr.detail_stream.write(detail)
            stderr.detail_stream.flush()
            if not verbose:
                stderr.terminal_stream.write(
                    f"[Warning] {category.__name__}: {' '.join(str(message).split())}\n"
                )
                stderr.terminal_stream.flush()

        warnings.showwarning = showwarning
        yield


def run_command(command, config_path, *, debug=False, verbose=False) -> int:
    """Run one workflow command."""
    start_time = time.time()
    requested_command = canonical_command(command)
    command, args = parse_selected_command(config_path, command)
    args.config_path = str(Path(config_path).resolve())
    args.invoked_command = command
    os.makedirs(args.result_dir, exist_ok=True)
    log_path = os.path.join(args.result_dir, "run.log")
    log_mode = "w" if command == "run" else "a"
    recorder = ProvenanceRecorder(args.result_dir, requested_command, config_path, args)

    with open(log_path, log_mode, encoding="utf-8") as log_file:
        stdout = Tee(sys.stdout, log_file, verbose=verbose)
        stderr = Tee(sys.stderr, log_file, verbose=verbose)
        with contextlib.ExitStack() as stack:
            stack.enter_context(contextlib.redirect_stdout(stdout))
            stack.enter_context(contextlib.redirect_stderr(stderr))
            stack.enter_context(warning_output(stderr, verbose=verbose))
            if log_mode == "a" and os.path.getsize(log_path) > 0:
                print_detail("\n" + "=" * 60)
            print(f"[main] Command: {requested_command}")
            if command != requested_command:
                print(f"[main] Input-selected operation: {command}")
            if config_path is not None:
                print_detail(f"[main] Reading arguments from {config_path}...")
            print_detail(f"[main] Logging output to {log_path}")
            print_detail("-" * 60)
            print_detail(format_config(args, command=command))
            print_detail("-" * 60)

            exit_code = 0
            execution = None
            failure = None
            try:
                manifest_path = recorder.start()
                print_detail(f"[main] Provenance manifest: {manifest_path}")
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
                print_detail("-" * 60)
                print(f"[main] {'Completed' if exit_code == 0 else 'Failed'}: {requested_command}")
                print(f"Total Execution Time: {format_duration(time.time() - start_time)}")
                print_detail("-" * 60)
                print(f"Run log saved to: {log_path}")

    return exit_code


def parse_selected_command(config_path, command):
    """Resolve input-selected operations before command-specific validation."""
    if command in {"fc2", "kappa"}:
        plan = parse_workflow_args(config_path, command="plot", input_command=command)
        try:
            selected = resolve_calculation_command(command, plan)
        except ValueError as exc:
            raise SystemExit(f"nepkappa {command}: error: {exc}") from None
    else:
        selected = command
    return selected, parse_workflow_args(config_path, command=selected)


def info_command(config_path, command="run") -> int:
    """Validate and print a workflow's settings without running it."""
    selected, args = parse_selected_command(config_path, command)
    if selected != command:
        print(f"Input-selected operation: {selected}")
    print(format_config(args, command=selected))
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
