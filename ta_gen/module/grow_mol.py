#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import json
import os


from ta_base.common.util import ArgumentParser
from ta_base.exception.ta_module_exception import TAModuleException

from ta_gen.utils.logger import LOGGER
from ta_gen.tool.crem_runner import CremRunner
from ta_gen.utils.common_utils import load_db_query, read_scaffolds


def schema_parser():
    parser = ArgumentParser(conflict_handler="resolve", add_help=False)
    group = parser.add_argument_group("Required arguments")
    group.add_argument("-w", "--work_dir", required=True, help="working directory.")
    group.add_argument("--task_index", type=int, required=True, help="index number of task")
    group.add_argument("--scaffolds_file", type=str, required=True, help="scaffold file")
    group.add_argument("--smarts_file", required=True, type=str)
    group.add_argument("--protect_id", required=True, type=int)
    group.add_argument("--max_replacements", required=True, type=int)
    group.add_argument("--input_smi", required=True, type=str)

    parser.add_argument("--result_file", default="", help="result file")

    CremRunner.add_to_argparse(parser)

    return parser


def parse_args():
    """Parse arguments and set default values."""
    parser = schema_parser()
    args, _ = parser.parse_known_args()

    try:
        args.db_query = load_db_query(args.db_query)
    except Exception as e:
        print("db_query: ", args.db_query)
        raise RuntimeError(f"Failed to load 'db_query' conditions: {e}")

    return args


def run_grow_mol(args):
    work_dir = args.work_dir
    if not os.path.exists(work_dir):
        os.makedirs(work_dir, exist_ok=True)

    scaffold = read_scaffolds(args.scaffolds_file)[0]

    crem_runner = CremRunner(
        scaffold=scaffold,
        db_config=args.db_config,
        db_engine=args.db_engine,
        radius=args.radius,
        min_atoms=args.min_atoms,
        max_atoms=args.max_atoms,
        num_cpus=args.num_cpus,
        work_dir=work_dir,
        db_query=args.db_query,
    )

    with open(args.input_smi, "r") as f:
        input_smi = f.readlines()[0].strip()

    with open(args.smarts_file, "r") as f:
        smarts = f.readlines()[0].strip()

    LOGGER.info("Start grow")
    LOGGER.info(f"input_smi: {input_smi}")
    LOGGER.info(f"smarts: {smarts}")
    LOGGER.info(f"protect_id: {args.protect_id}")
    LOGGER.info(f"max_replacements: {args.max_replacements}")
    LOGGER.info(f"db_query: {args.db_query}")
    LOGGER.info(f"db_config: {args.db_config}")
    LOGGER.info(f"db_engine: {args.db_engine}")

    try:
        max_replacements = (
            args.max_replacements if args.max_replacements > 0 else None
        )  # -1 = no limit
        out_smis = crem_runner.additional_grow(
            input_smi, smarts, args.protect_id, max_replacements
        )
    except Exception as e:
        print(e)
        raise RuntimeError(f"Failed to grow molecules: {e}")

    LOGGER.info(f"Total number of generated molecules: {len(out_smis)}")

    result_file = os.path.join(work_dir, f"task_{args.task_index}.json")
    with open(result_file, "wt") as f:
        json.dump(list(out_smis), f, indent=4)

    return out_smis


def main():  # pragma: no cover
    args = parse_args()
    run_grow_mol(args)


if __name__ == "__main__":  # pragma: no cover
    main()