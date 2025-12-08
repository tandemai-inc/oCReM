#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import copy
import json
import os
import pickle

from rdkit import Chem
from ta_base.common.util import ArgumentParser
from ta_base.exception.ta_module_exception import TAModuleException

from ta_gen.utils.logger import LOGGER
from ta_gen.tool.crem_runner import CremRunner
from ta_gen.tool.post_processing import PostProcessor
from ta_gen.utils.common_utils import (
    get_smiles_list_from_csv,
    load_db_query,
    read_scaffolds,
)
from ta_gen.utils.const import MAXINUM_NUM_OF_MOLS_TO_GROW, MAXINUM_NUM_OF_OUTPUT_MOLS
from ta_gen.utils.prop_filter import Filter


def schema_parser():
    parser = ArgumentParser(conflict_handler="resolve", add_help=False)
    group = parser.add_argument_group("Required arguments")
    group.add_argument("-w", "--work_dir", required=True, help="working directory.")
    group.add_argument(
        "--result_out", type=str, required=True, help="result file path", default="result.csv"
    )
    group.add_argument("--scaffolds_file", type=str, required=True, help="scaffold file")

    parser.add_argument("--result_file", default="", help="result file")

    parser.add_argument("--run_ta_prop", action="store_true")
    parser.add_argument("--num_cpus", type=int, default=2)
    parser.add_argument("--cleanup", action="store_true")
    parser.add_argument("--max_mols_gen", type=int, default=MAXINUM_NUM_OF_OUTPUT_MOLS)
    parser.add_argument("--user_provided_output_file", type=str, default="")

    Filter.add_to_argparse(parser)
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


def run_init_grow(args):
    d_args = vars(args)
    os.makedirs(args.work_dir, exist_ok=True)
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
        db_query=args.db_query,
    )

    user_provided_output_file = args.user_provided_output_file
    if user_provided_output_file:
        LOGGER.info("Read user provided output file...")
        out_smis = get_smiles_list_from_csv(user_provided_output_file)
        finished = True
    else:
        try:
            mol_h, protect_id_list, attach_id_list = crem_runner.prepare()

            finished = len(protect_id_list) == 1

            # if final round,
            # limit is max_mols_gen else MAXINUM_NUM_OF_MOLS_TO_GROW(for next round)
            max_replacements = args.max_mols_gen if finished else MAXINUM_NUM_OF_MOLS_TO_GROW
            if max_replacements <= 0:
                max_replacements = None
            out_smis = crem_runner.grow(
                mol_h,
                protect_id_list,
                attach_id_list,
                max_replacements=max_replacements,
            )

            finished |= len(out_smis) == 0

        except Exception as e:
            print(e)
            raise RuntimeError(f"Failed to initialize molecules: {e}")

    if finished:
        PostProcessor.post_process(
            result_data={"smiles": out_smis},
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
            "finished": finished,
        }
    else:
        smarts_list = [""]

        pre_mol = mol_h
        for i in range(1, len(protect_id_list)):
            new_mol = copy.deepcopy(pre_mol)
            pre_attach_id = attach_id_list[i - 1]
            # convert previous attachment to *
            # get smarts of scaffold
            atom = new_mol.GetAtomWithIdx(pre_attach_id)
            atom.SetAtomicNum(0)
            smarts = Chem.MolToSmiles(new_mol, canonical=False)
            smarts_list.append(smarts)
            pre_mol = new_mol

        result_dict = {
            "finished": finished,
            "out_smis": out_smis,
            "protect_id_list": protect_id_list,
            "attach_id_list": attach_id_list,
            "smarts_list": smarts_list,
        }

    with open(args.result_file, "wb") as f:
        pickle.dump(result_dict, f)
        f.flush()

    return result_dict


def main():  # pragma: no cover
    args = parse_args()
    run_init_grow(args)


if __name__ == "__main__":  # pragma: no cover
    main()