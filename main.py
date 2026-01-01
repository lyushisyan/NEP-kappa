#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ==============================================================================
# Author: Shixian Liu, Fei Yin
# Date: 2025
# Description: Main execution script.
# ==============================================================================

import argparse
import shlex
import os
import sys
import time

try:
    from workflow import NEPPhononWorkflow
except ImportError:
    print("Error: Could not import 'workflow.py'.")
    sys.exit(1)

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
    parser.add_argument("--temps", type=float, nargs=3, default=[100, 1000, 100], help="Temperatures")
    parser.add_argument("--method", type=str, choices=["lbte", "rta"], default="lbte", help="Method")
    parser.add_argument("--wigner", type=str2bool, nargs="?", const=True, default=False, help="Wave-like contribution")

    return parser

def parse_input_file(filename):
    args_list = []
    if not os.path.exists(filename): return []
    print(f"[main] Reading arguments from {filename}...")
    with open(filename, "r", encoding="utf-8") as f:
        lines = f.readlines()
    clean_text = ""
    for line in lines:
        content = line.split('#')[0].strip()
        if content: clean_text += " " + content
    return shlex.split(clean_text)

def main():
    start_time = time.time()
    
    parser = initialise_parser()
    target_file = "input.txt"
    if len(sys.argv) > 1 and not sys.argv[1].startswith("-"):
        target_file = sys.argv[1]

    if os.path.exists(target_file):
        args = parser.parse_args(parse_input_file(target_file))
    else:
        if target_file != "input.txt":
            print(f"[Warning] '{target_file}' not found. Falling back to CLI.")
        args = parser.parse_args()

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

if __name__ == "__main__":
    main()