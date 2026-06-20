"""Command-line interface for NEP-kappa."""

from __future__ import annotations

import argparse
import contextlib
import os
import sys
import time

from nepkappa import __version__
from nepkappa.config import (
    format_compare_config,
    format_config,
    parse_compare_args,
    parse_workflow_args,
)


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

    if args.command in {"run", "relax", "fc", "kappa", "plot"}:
        return run_command(args.command, args.config)
    if args.command == "compare":
        return compare_command(args.config)
    if args.command == "info":
        return info_command(args.config)

    parser.error("missing command")
    return 2


def build_parser() -> argparse.ArgumentParser:
    """Build the top-level parser."""
    parser = argparse.ArgumentParser(
        prog="nepkappa",
        description="NEP-assisted lattice thermal conductivity workflow.",
    )
    parser.add_argument("--version", action="version", version=f"nepkappa {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser(
        "run", help="Run relax, fc, and kappa in sequence."
    )
    run_parser.add_argument("config", help="YAML input file")

    relax_parser = subparsers.add_parser(
        "relax", help="Relax the input structure only."
    )
    relax_parser.add_argument("config", help="YAML input file")

    fc_parser = subparsers.add_parser(
        "fc", help="Generate fc2.hdf5, fc3.hdf5, and phono3py_disp.yaml."
    )
    fc_parser.add_argument("config", help="YAML input file")

    kappa_parser = subparsers.add_parser(
        "kappa", help="Compute thermal conductivity from existing force constants."
    )
    kappa_parser.add_argument("config", help="YAML input file")

    plot_parser = subparsers.add_parser(
        "plot", help="Plot phonon and thermal-transport results."
    )
    plot_parser.add_argument("config", help="YAML input file")

    compare_parser = subparsers.add_parser(
        "compare", help="Compare DFT and NEP phonon/thermal-transport results."
    )
    compare_parser.add_argument("config", help="YAML comparison input file")

    info_parser = subparsers.add_parser(
        "info", help="Print parsed workflow settings without running."
    )
    info_parser.add_argument("config", help="YAML input file")

    return parser


def run_command(command, config_path) -> int:
    """Run one workflow command."""
    start_time = time.time()
    args = parse_workflow_args(config_path)
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
                elif command == "plot":
                    from nepkappa.plot import plot_results

                    plot_results(args)
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


def info_command(config_path) -> int:
    """Print parsed config values without running the workflow."""
    args = parse_workflow_args(config_path)
    print(format_config(args))
    return 0


def compare_command(config_path) -> int:
    """Run DFT-vs-NEP comparison plotting."""
    start_time = time.time()
    args = parse_compare_args(config_path)
    os.makedirs(args.compare_dir, exist_ok=True)
    log_path = os.path.join(args.compare_dir, "compare.log")

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
            try:
                from nepkappa.plot import compare_results

                compare_results(args)
            except Exception as exc:
                exit_code = 1
                print(f"\n[Error] Comparison failed: {exc}")
                import traceback

                traceback.print_exc()
            finally:
                print("-" * 60)
                print(f"Total Execution Time: {format_duration(time.time() - start_time)}")
                print("-" * 60)

    print("-" * 60)
    print(f"Compare log saved to: {log_path}")
    return exit_code


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
