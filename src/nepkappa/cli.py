"""Command-line interface for NEP-kappa."""

from __future__ import annotations

import argparse
import contextlib
import os
import sys
import time
from pathlib import Path

from nepkappa import __version__
from nepkappa.config import format_config, parse_workflow_args


class Tee:
    """Write text to multiple streams, e.g. terminal and run.log."""

    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for stream in self.streams:
            stream.write(data)
            stream.flush()

    def flush(self):
        for stream in self.streams:
            stream.flush()


def main(argv: list[str] | None = None) -> int:
    """Run the NEP-kappa command-line interface."""
    argv = normalize_argv(list(sys.argv[1:] if argv is None else argv))
    parser = build_parser()
    args, extra_args = parser.parse_known_args(argv)

    if args.command in {"run", "relax", "fc", "kappa"}:
        return run_command(args.command, args.config, extra_args)
    if args.command == "info":
        return info_command(args.config, extra_args)

    parser.error("missing command")
    return 2


def normalize_argv(argv: list[str]) -> list[str]:
    """Preserve the concise form ``nepkappa input.yaml``."""
    commands = {"run", "relax", "fc", "kappa", "info", "-h", "--help", "--version"}
    if not argv:
        return ["run"]
    if argv[0] in commands:
        return argv
    if Path(argv[0]).suffix.lower() in {".yaml", ".yml", ".txt"}:
        return ["run"] + argv
    return argv


def build_parser() -> argparse.ArgumentParser:
    """Build the top-level software-style CLI parser."""
    parser = argparse.ArgumentParser(
        prog="nepkappa",
        description="NEP-assisted lattice thermal conductivity workflow.",
    )
    parser.add_argument("--version", action="version", version=f"nepkappa {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="Run a NEP-kappa workflow.")
    run_parser.add_argument("config", nargs="?", help="YAML or legacy text input file")

    relax_parser = subparsers.add_parser(
        "relax", help="Relax the input structure only."
    )
    relax_parser.add_argument("config", nargs="?", help="YAML or legacy text input file")

    fc_parser = subparsers.add_parser(
        "fc", help="Generate force constants only."
    )
    fc_parser.add_argument("config", nargs="?", help="YAML or legacy text input file")

    kappa_parser = subparsers.add_parser(
        "kappa", help="Compute thermal conductivity from existing force constants."
    )
    kappa_parser.add_argument("config", nargs="?", help="YAML or legacy text input file")

    info_parser = subparsers.add_parser(
        "info", help="Print parsed workflow settings without running."
    )
    info_parser.add_argument("config", help="YAML or legacy text input file")

    return parser


def run_command(command, config_path=None, extra_args=None) -> int:
    """Run one workflow command."""
    start_time = time.time()
    args = parse_workflow_args(config_path, extra_args)
    os.makedirs(args.result_dir, exist_ok=True)
    log_path = os.path.join(args.result_dir, "run.log")
    log_mode = "w" if command == "run" else "a"

    with open(log_path, log_mode, encoding="utf-8") as log_file:
        stdout = Tee(sys.stdout, log_file)
        stderr = Tee(sys.stderr, log_file)
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            if log_mode == "a" and os.path.getsize(log_path) > 0:
                print("\n" + "=" * 60)
            print(f"[main] Command: {command}")
            if config_path is not None:
                print(f"[main] Reading arguments from {config_path}...")
            print(f"[main] Logging output to {log_path}")
            print("-" * 60)
            print(format_config(args))
            print("-" * 60)

            exit_code = 0
            try:
                from nepkappa.workflow import NEPPhononWorkflow

                workflow = NEPPhononWorkflow(args)
                if command == "relax":
                    workflow.run_relax()
                elif command == "fc":
                    workflow.run_force_constants()
                elif command == "kappa":
                    workflow.run_kappa()
                else:
                    workflow.run()
            except Exception as exc:
                exit_code = 1
                print(f"\n[Error] Workflow execution failed: {exc}")
                import traceback

                traceback.print_exc()
            finally:
                print("-" * 60)
                print(f"Total Execution Time: {format_duration(time.time() - start_time)}")
                print("-" * 60)

    print("-" * 60)
    print(f"Run log saved to: {log_path}")
    return exit_code


def info_command(config_path, extra_args=None) -> int:
    """Print parsed config values without running the workflow."""
    args = parse_workflow_args(config_path, extra_args)
    print(format_config(args))
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
