#!/usr/bin/env python
# -*- coding:utf-8 -*-

import argparse
import os
import pathlib
import pickle

import yaml
from src.utils.cmd import cmd
from src.utils.common_utils import dict_to_cmdline
from src.utils.const import (MAXINUM_NUM_OF_MOLS_TO_GROW,
                             MAXINUM_NUM_OF_OUTPUT_MOLS)

MODULE_PATH = pathlib.Path(__file__).parents[1] / "module"


def parse_args():  # pragma: no cover
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--in_config", help="config file")
    return parser.parse_args()


def read_config(config_file):
    """
    read config yaml file and parser it.

    Args:
        None

    Returns:
        DotDict: most like dict, support dot call

    """
    try:
        with open(config_file, "r") as f:
            return yaml.load(f, Loader=yaml.FullLoader)
    except IOError:
        raise Exception("Cannot read options from config file %s." % config_file)


def init_grow(paras):
    work_dir = "./scr/init_grow"
    result_file = f"./scr/init_grow/init_grow.pkl"

    run_ta_prop = paras["parameter"]["master"]["run_ta_prop"]
    paras = {
        "work_dir": work_dir,
        "result_file": result_file,
        "result_out": paras["parameter"]["master"]["result_out"],
        "scaffolds_file": paras["parameter"]["master"]["scaffolds_file"],
        "run_ta_prop": run_ta_prop,
        "num_cpus": paras["parameter"]["master"]["num_cpus"],
        "max_mols_gen": paras["parameter"]["master"]["max_mols_gen"],
        "user_provided_output_file": paras["parameter"]["master"][
            "user_provided_output_file"
        ],
    }
    paras.update(paras["parameter"]["crem"])
    if run_ta_prop:
        paras.update(paras["parameter"]["ta_prop"])

    return_code, stdout, stderr = cmd(
        f"python {MODULE_PATH}/init_grow.py {dict_to_cmdline(paras)}"
    )
    if return_code != 0:
        raise Exception(
            f"init_grow failed, return code: {return_code}, stdout: {stdout}, stderr: {stderr}"
        )

    with open(result_file, "rb") as f:
        result_dict = pickle.load(f)

    crem_rg_finished = result_dict["finished"]
    crem_rg_next_iter = not crem_rg_finished
    return crem_rg_finished, crem_rg_next_iter, result_dict


def grow_mol(pre_grow_result, index, paras):
    work_dir = f"./scr/grow_mol/{index}/outputs"
    result_file = f"./scr/grow_mol/{index}/grow_mol.pkl"

    smarts = pre_grow_result["smarts_list"][index]
    smarts_file = os.path.join(work_dir, "in.smarts")
    with open(smarts_file, "w") as f:
        f.write(smarts)

    protect_id = pre_grow_result["protect_id"][index]

    if index != len(protect_id) - 1:  # not final grow mol
        max_replacements = MAXINUM_NUM_OF_MOLS_TO_GROW  # intermediate grow
    else:
        max_replacements = (
            paras["parameter"]["master"].get(
                "max_mols_gen", MAXINUM_NUM_OF_OUTPUT_MOLS  # final grow
            ),
        )

    pre_out_smiles = pre_grow_result["out_smis"]
    for i, smi in enumerate(pre_out_smiles):
        task_id = f"grow_mol_{index}_{i}"
        input_smi = os.path.join(work_dir, f"{task_id}.smi")
        with open(input_smi, "w") as f:
            f.write(smi)

        paras = {
            "work_dir": work_dir,
            "result_file": result_file,
            "task_index": i,
            "input_smi": input_smi,
            "scaffolds_file": paras["parameter"]["master"]["scaffolds_file"],
            "smarts_file": smarts_file,
            "protect_id": protect_id,
            "max_replacements": max_replacements,
        }
        paras.update(paras["parameter"]["crem"])

    return_code, stdout, stderr = cmd(
        f"python {MODULE_PATH}/init_grow.py {dict_to_cmdline(paras)}"
    )
    if return_code != 0:
        raise Exception(
            f"grow mol failed, return code: {return_code}, stdout: {stdout}, stderr: {stderr}"
        )

    with open(result_file, "rb") as f:
        result_dict = pickle.load(f)

    return result_dict


def post_grow_mol(grow_result, index, paras, protect_id):
    grow_mol_dir = f"./scr/grow_mol/{index}/outputs"
    work_dir = f"./scr/grow_mol/{index}"
    result_file = f"./scr/init_grow/init_grow.pkl"

    run_ta_prop = paras["parameter"]["master"]["run_ta_prop"]

    max_mols_gen = paras["parameter"]["master"].get(
        "max_mols_gen", MAXINUM_NUM_OF_OUTPUT_MOLS
    )

    paras = {
        "work_dir": work_dir,
        "result_file": result_file,
        "result_out": paras["parameter"]["master"]["result_out"],
        "scaffolds_file": paras["parameter"]["master"]["scaffolds_file"],
        "grow_mol_dir": grow_mol_dir,
        "run_ta_prop": run_ta_prop,
        "num_cpus": paras["parameter"]["master"]["num_cpus"],
        "max_mols_gen": max_mols_gen,
        "last_iteration": index == len(protect_id) - 1,
    }
    paras.update(paras["parameter"]["crem"])
    if run_ta_prop:
        paras.update(paras["parameter"]["ta_prop"])
    return_code, stdout, stderr = cmd(
        f"python {MODULE_PATH}/init_grow.py {dict_to_cmdline(paras)}"
    )
    if return_code != 0:
        raise Exception(
            f"init_grow failed, return code: {return_code}, stdout: {stdout}, stderr: {stderr}"
        )

    with open(result_file, "rb") as f:
        result_dict = pickle.load(f)

    crem_rg_finished = result_dict["finished"]
    crem_rg_next_iter = not crem_rg_finished
    return crem_rg_finished, crem_rg_next_iter, result_dict


def main(paras):
    print("Create Crem RGroup Enumeration tasks")
    crem_rg_finished, crem_rg_next_iter, pre_grow_result = init_grow(paras)
    index = 1
    while not crem_rg_finished:
        if crem_rg_next_iter:
            protect_id = pre_grow_result["protect_id"]
            grow_result = grow_mol(pre_grow_result, index, paras)
            crem_rg_finished, crem_rg_next_iter, pre_grow_result = post_grow_mol(
                grow_result, index, paras, protect_id
            )
            index += 1
        yield False

    yield True


if __name__ == "__main__":
    args = parse_args()
    paras = read_config(args.in_config)
    main(paras)
