#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import argparse
import json
import os
import pickle
import random

from ta_gen.tool.crem_runner import CremRunner
from ta_gen.tool.post_processing import PostProcessor
from ta_gen.utils.common_utils import read_scaffolds
from ta_gen.utils.const import (MAXINUM_NUM_OF_MOLS_TO_GROW,
                                MAXINUM_NUM_OF_OUTPUT_MOLS)
from ta_gen.utils.logger import LOGGER
from ta_gen.utils.prop_filter import Filter


def schema_parser():
    parser = argparse.ArgumentParser(conflict_handler="resolve", add_help=False)
    group = parser.add_argument_group("Required arguments")
    group.add_argument("-w", "--work_dir", required=True, help="working directory.")
    group.add_argument(
        "--result_out",
        type=str,
        required=True,
        help="result file path",
        default="result.csv",
    )
    group.add_argument(
        "--scaffolds_file", type=str, required=True, help="scaffold file"
    )
    group.add_argument("--grow_mol_dir", type=str, required=True)

    parser.add_argument("--result_file", default="", help="result file")

    parser.add_argument("--run_ta_prop", action="store_true")
    parser.add_argument("--num_cpus", type=int, default=2)
    parser.add_argument("--cleanup", action="store_true")
    parser.add_argument("--max_mols_gen", type=int, default=MAXINUM_NUM_OF_OUTPUT_MOLS)
    parser.add_argument("--last_iteration", action="store_true")

    Filter.add_to_argparse(parser)
    CremRunner.add_to_argparse(parser)

    return parser


def parse_args():
    """Parse arguments and set default values."""
    parser = schema_parser()
    args, _ = parser.parse_known_args()
    return args


def get_all_produced_mols(grow_mol_dir):
    all_produced_mols = set()
    for file in os.listdir(grow_mol_dir):
        if file.endswith(".json"):
            with open(os.path.join(grow_mol_dir, file), "r") as f:
                all_produced_mols |= set(json.load(f))
    return list(all_produced_mols)


def run_post_grow(args):
    d_args = vars(args)
    all_produced_mols = get_all_produced_mols(args.grow_mol_dir)
    scaffold = read_scaffolds(args.scaffolds_file)[0]

    crem_runner = CremRunner(
        scaffold=scaffold,
        db_config=args.db_config,
        db_engine=args.db_engine,
        radius=args.radius,
        min_atoms=args.min_atoms,
        max_atoms=args.max_atoms,
        result_out=args.result_out,
        num_cpus=args.num_cpus,
        work_dir=args.work_dir,
        max_mols_gen=args.max_mols_gen,
    )

    if args.last_iteration or len(all_produced_mols) == 0:
        PostProcessor.post_process(
            result_data={"smiles": all_produced_mols},
            run_ta_prop=args.run_ta_prop,
            tool_name="crem",
            scaffold=crem_runner.scaffold,
            attachment_points_mapping=crem_runner.attachment_points_mapping,
            num_cpus=args.num_cpus,
            max_mols_gen=args.max_mols_gen,
            result_file=args.result_out,
            ta_prop_args=d_args,
        )

        result_dict = {
            "finished": True,
        }

    else:
        if len(all_produced_mols) > MAXINUM_NUM_OF_MOLS_TO_GROW:
            LOGGER.info(
                f"Sample from out_smis: {len(all_produced_mols)} -> {MAXINUM_NUM_OF_MOLS_TO_GROW}"
            )
            all_produced_mols = random.sample(
                all_produced_mols, MAXINUM_NUM_OF_MOLS_TO_GROW
            )

        result_dict = {
            "finished": False,
            "out_smis": all_produced_mols,
        }

        LOGGER.info(f"Total output: {len(all_produced_mols)}")

    with open(args.result_file, "wb") as f:
        pickle.dump(result_dict, f)
        f.flush()

    return result_dict


def main():  # pragma: no cover
    args = parse_args()
    run_post_grow(args)


if __name__ == "__main__":  # pragma: no cover
    main()
