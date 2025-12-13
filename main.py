#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import shlex

from nep_phonon import PhononWorkflow


def str2bool(v: str) -> bool:
    """
    把字符串解析成 bool:
      true/1/yes/y/t -> True
      false/0/no/n/f -> False
    """
    if isinstance(v, bool):
        return v
    v = v.strip().lower()
    if v in ("true", "1", "yes", "y", "t"):
        return True
    if v in ("false", "0", "no", "n", "f"):
        return False
    raise argparse.ArgumentTypeError(f"Invalid boolean value: '{v}'")


def initialise_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="NEP + HiPhive + phono3py phonon/kappa workflow"
    )

    # ===== 基本输入 =====
    parser.add_argument("--poscar", default="POSCAR",
                        help="input POSCAR file for primitive cell")
    parser.add_argument("--nep_model", default="Si_2025_XuKe.txt",
                        help="NEP potential file")

    # ===== 统一超胞尺寸：训练 + phonopy 同一个 dim =====
    parser.add_argument(
        "--dim",
        type=int,
        nargs=3,
        default=[4, 4, 1],
        help="supercell repeat (nx ny nz), used for both training supercell and phonopy/phono3py"
    )

    # ===== 训练结构参数 =====
    parser.add_argument("--n_structures", type=int, default=500,
                        help="number of rattled structures for training")
    parser.add_argument("--rattle_std", type=float, default=0.02,
                        help="std dev for random displacements (Å)")
    parser.add_argument("--min_dist", type=float, default=2.0,
                        help="minimum interatomic distance in rattling (Å)")
    parser.add_argument("--seed", type=int, default=2025,
                        help="random seed")

    # 这里改成接收一个布尔值，而不是 action="store_true"
    parser.add_argument(
        "--do_relax",
        type=str2bool,
        nargs="?",          # 允许 "--do_relax"（默认 True）或 "--do_relax true/false"
        const=True,
        default=False,
        help="whether to relax primitive cell with NEP before training (true/false)"
    )

    # ===== hiphive 截断 =====
    parser.add_argument(
        "--cutoffs",
        type=float,
        nargs="+",
        default=[5.0, 4.0, 0.0],
        help="cluster cutoffs for HiPhive, e.g. --cutoffs 5.0 4.0 0.0"
    )

    # ===== phono3py / κ 计算参数：统一 mesh =====
    parser.add_argument(
        "--mesh",
        type=int,
        nargs=3,
        default=[15, 15, 1],
        help="q-point mesh (mx my mz) for phono3py"
    )

    parser.add_argument(
        "--temps",
        type=int,
        nargs="+",
        default=[100, 200, 300, 400],
        help="temperatures for kappa(T), e.g. --temps 100 200 300 400"
    )

    return parser


def main():
    parser = initialise_parser()

    # 从 input.txt 读取“命令行风格”的参数
    with open("input.txt", "r", encoding="utf-8") as infile:
        text = infile.read()
        # 如果是 Python 3.11，可以用 comments=True 支持 # 注释
        # args_list = shlex.split(text, comments=True)
        args_list = shlex.split(text)

    args = parser.parse_args(args_list)
    print("[main] Parsed args from input.txt:")
    print(args)

    wf = PhononWorkflow(args)
    wf.run_all()


if __name__ == "__main__":
    main()
