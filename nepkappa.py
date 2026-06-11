#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ==============================================================================
# Author: Shixian Liu, Fei Yin
# Date: 2025
# Description: Main execution script.
# ==============================================================================

import argparse
import contextlib
import shlex
import os
import sys
import time

from workflow import NEPPhononWorkflow

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

def str2bool(v):
    if isinstance(v, bool): return v
    v = str(v).strip().lower()
    if v in ("true", "1", "yes", "y", "t"): return True
    if v in ("false", "0", "no", "n", "f"): return False
    raise argparse.ArgumentTypeError(f"Invalid boolean value: '{v}'")

def initialise_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="NEP + HiPhive + Phono3py Workflow")
    
    parser.add_argument("--poscar", default="POSCAR", help="Input structure file")
    parser.add_argument("--nep_model", required=True, help="Path to NEP model file")
    parser.add_argument("--do_relax", type=str2bool, nargs="?", const=True, default=False, help="Relax structure")
    parser.add_argument("--dim", type=int, nargs=3, default=[4, 4, 1], help="Supercell dimension")
    
    # HiPhive settings
    parser.add_argument("--use_hiphive", type=str2bool, nargs="?", const=True, default=False, help="Use HiPhive")
    parser.add_argument("--n_structures", type=int, default=50, help="[HiPhive] Number of structures")
    parser.add_argument("--rattle_std", type=float, default=0.03, help="[HiPhive] Rattle std")
    parser.add_argument("--cutoffs", type=float, nargs="+", default=[5.0], help="[HiPhive] Cutoffs")
    parser.add_argument("--min_dist", type=float, default=2.0, help="[HiPhive] Minimum atomic distance (d_min) for rattling")

    # Phono3py settings
    parser.add_argument("--mesh", type=int, nargs=3, default=[21, 21, 1], help="Phono3py mesh")
    parser.add_argument("--temps", type=float, nargs="+", default=[100, 1000, 100], help="Temperatures")
    parser.add_argument("--fc2fc3", type=str2bool, nargs="?", const=True, default=False, help="Fitting force constants")    
    parser.add_argument("--method", type=str, choices=["lbte", "rta"], default="lbte", help="Method")
    parser.add_argument("--wigner", type=str2bool, nargs="?", const=True, default=False, help="Wave-like contribution")
    parser.add_argument("--progress", type=str2bool, nargs="?", const=True, default=True, help="Show progress bars and timing summaries")
    parser.add_argument("--result_dir", default="result", help="Directory for generated outputs and run.log")
    parser.add_argument("--output_name", default=None, help="Optional output name for the final kappa HDF5 file")

    return parser

def parse_input_file(filename):
    if not os.path.exists(filename): return []
    with open(filename, "r", encoding="utf-8") as f:
        lines = f.readlines()
    clean_text = ""
    for line in lines:
        content = line.split('#')[0].strip()
        if content: clean_text += " " + content
    return shlex.split(clean_text)

def validate_temps(temps):
    if not isinstance(temps, (list, tuple)):
        return "--temps must be a list of floats."
    n = len(temps)
    if n not in (1, 3):
        return (
            f"--temps expects 1 value (single temperature) or 3 values "
            f"(tmin tmax tstep), got {n}: {temps}"
        )
    return None

def main():
    start_time = time.time()
    
    parser = initialise_parser()
    target_file = "input.txt"
    if len(sys.argv) > 1 and not sys.argv[1].startswith("-"):
        target_file = sys.argv[1]

    if os.path.exists(target_file):
        args = parser.parse_args(parse_input_file(target_file))
        input_source = target_file
        missing_input_warning = None
    else:
        input_source = None
        missing_input_warning = None
        if target_file != "input.txt":
            missing_input_warning = f"[Warning] '{target_file}' not found. Falling back to CLI."
        args = parser.parse_args()

    temps_error = validate_temps(args.temps)
    if temps_error:
        parser.error(temps_error)

    os.makedirs(args.result_dir, exist_ok=True)
    log_path = os.path.join(args.result_dir, "run.log")

    with open(log_path, "w", encoding="utf-8") as log_file:
        stdout = Tee(sys.stdout, log_file)
        stderr = Tee(sys.stderr, log_file)
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            if missing_input_warning is not None:
                print(missing_input_warning)
            if input_source is not None:
                print(f"[main] Reading arguments from {input_source}...")
            print(f"[main] Logging output to {log_path}")
            print("-" * 60)
            print("Running Workflow with configuration:")
            for arg, value in vars(args).items():
                print(f"  {arg:<15} : {value}")
            print("-" * 60)

            try:
                wf = NEPPhononWorkflow(args)
                wf.run()
            except Exception as e:
                print(f"\n[Error] Workflow execution failed: {e}")
                import traceback
                traceback.print_exc()
            finally:
                end_time = time.time()
                elapsed_seconds = end_time - start_time

                hours = int(elapsed_seconds // 3600)
                minutes = int((elapsed_seconds % 3600) // 60)
                seconds = elapsed_seconds % 60

                print("-" * 60)
                print(f"Total Execution Time: {hours}h {minutes}m {seconds:.2f}s")
                print("-" * 60)

    if os.path.exists(log_path):
        print("-" * 60)
        print(f"Run log saved to: {log_path}")

if __name__ == "__main__":
    main()
