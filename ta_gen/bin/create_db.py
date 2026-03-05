#!/usr/bin/env python
# -*- coding:utf-8 -*-

import argparse
import csv
import os
import subprocess
import sys
from collections import Counter
from functools import partial
from itertools import permutations
from multiprocessing import Pool
from queue import Queue
from threading import Thread

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMMPA
from tqdm import tqdm

from ta_gen.db import create_db_manager
from ta_gen.utils.mol_context import (combine_core_env_to_rxn_smarts,
                                      get_std_context_core_permutations)


def schema_parser():
    parser = argparse.ArgumentParser(description="Create database")
    parser.add_argument(
        "-i", "--input_file", type=str, required=True, help="Input compounds"
    )
    parser.add_argument(
        "-o",
        "--out",
        metavar="output.txt",
        required=False,
        default="output.txt",
        help="output file of fragmented molecules.",
    )
    parser.add_argument(
        "--debug", action="store_true", default=False, help="print debug message."
    )
    parser.add_argument("--chunk_size", type=int, default=100, help="Chunk size")
    parser.add_argument(
        "-c",
        "--ncpu",
        metavar="NUMBER",
        required=False,
        default=1,
        type=int,
        help="number of cpus used for computation. Default: 1.",
    )
    parser.add_argument(
        "-s",
        "--sep",
        metavar="STRING",
        required=False,
        default=",",
        help="separator in input file. Default: Tab.",
    )
    parser.add_argument(
        "-d",
        "--sep_out",
        metavar="STRING",
        required=False,
        default=",",
        help="separator in the output file. Default: comma",
    )
    parser.add_argument(
        "-m",
        "--mode",
        metavar="INTEGER",
        required=False,
        default=0,
        choices=[0, 1, 2],
        type=int,
        help="fragmentation mode: 0 - all atoms constitute a fragment, 1- heavy atoms only, "
        "2 - hydrogen atoms only. Default: 0.",
    )

    group = parser.add_argument_group("Fragment to Environment Parameters")
    group.add_argument(
        "--max_heavy_atoms",
        metavar="NUMBER",
        required=False,
        default=20,
        help="maximum number of heavy atoms in cores. If the number of atoms exceeds the limit "
        "fragment will be discarded. Default: 20.",
    )
    parser.add_argument(
        "--radius",
        metavar="NUMBER",
        required=False,
        type=int,
        default=1,
        help="radius of molecular context (in bonds) which will be taken into account. Default: 1.",
    )
    parser.add_argument(
        "--keep_stereo",
        action="store_true",
        default=False,
        help="set this flag if you want to keep stereo in context and core parts.",
    )
    group = parser.add_argument_group("Database Parameters")
    group.add_argument(
        "--db_type",
        metavar="STRING",
        required=False,
        default="sqlite",
        choices=["sqlite", "postgres"],
        help="database type. Default: sqlite.",
    )
    group.add_argument(
        "--db_path",
        metavar="STRING",
        required=False,
        default="ocrem.db",
        help="database path. Default: :memory: (in-memory database).",
    )
    group.add_argument(
        "--ini_file",
        metavar="STRING",
        required=False,
        default=None,
        help="postgresql ini file. Default: None.",
    )
    group.add_argument(
        "--reset_db",
        action="store_true",
        default=False,
        help="set this flag if you want to reset the database.",
    )
    return parser


def parse_args():
    parser = schema_parser()
    args, _ = parser.parse_known_args()
    return args


def preprocess_input_file(input_file):
    """Preprocess input file"""
    # remove duplicated
    name, ext = os.path.splitext(input_file)
    output_file = f"{name}_deduped{ext}"
    if ext == ".csv":
        os.system(f'( head -n 1 "{input_file}"; tail -n +2 "{input_file}" | sort -u ) > "{output_file}"')
    else:
        os.system(f"sort -u {input_file} -o {output_file}")
    return output_file


def count_total_rows(input_file):
    """Count total number of rows in input file using wc -l"""
    result = subprocess.run(["wc", "-l", input_file], capture_output=True, text=True)
    total_rows = int(result.stdout.split()[0])
    _, ext = os.path.splitext(input_file)
    if ext == ".csv":
        total_rows -= 1
    return total_rows


def read_chunks(input_file, chunk_size, sep):
    _, ext = os.path.splitext(input_file)
    if ext == ".csv":
        chunks = pd.read_csv(input_file, sep=sep, chunksize=chunk_size)
    elif ext == ".smi":
        chunks = pd.read_csv(
            input_file,
            header=None,
            names=["smiles", "smi_id"],
            sep=sep,
            na_filter=False,
            chunksize=chunk_size,
        )
    else:
        chunks = pd.read_csv(
            input_file,
            header=None,
            names=["smiles", "smi_id"],
            sep=sep,
            na_filter=False,
            chunksize=chunk_size,
        )
    return chunks


def calc_mp(env, core):
    sma = combine_core_env_to_rxn_smarts(core, env, False)
    if core.count("*") == 2:
        mol = Chem.MolFromSmiles(core, sanitize=False)
        mat = Chem.GetDistanceMatrix(mol)
        ids = []
        for a in mol.GetAtoms():
            if not a.GetAtomicNum():
                ids.append(a.GetIdx())
        dist2 = mat[ids[0], ids[1]]
    else:
        dist2 = 0
    return sma, int(dist2)


def frag_to_env(smi, core, contexts, max_heavy_atoms, radius, keep_stereo):
    results = []
    if not core and not contexts:
        return results

    if not core:  # one split
        residues = contexts.split(".")
        if len(residues) == 2:
            for context, core in permutations(residues, 2):
                if context == "[H][*:1]":  # ignore such cases
                    continue

                mm = Chem.MolFromSmiles(core, sanitize=False)
                num_heavy_atoms = mm.GetNumHeavyAtoms() if mm else float("inf")
                if num_heavy_atoms <= max_heavy_atoms:
                    env, cores = get_std_context_core_permutations(
                        context, core, radius, keep_stereo
                    )
                    if env and cores:
                        sma, dist2 = calc_mp(env, cores[0])
                        results.append((env, cores[0], num_heavy_atoms, sma, dist2))
        else:
            sys.stderr.write(
                f"more than two fragments in context ({contexts}) where core is empty for smiles: {smi}\n"
            )
            sys.stderr.flush()
    else:  # two or more splits
        mm = Chem.MolFromSmiles(core, sanitize=False)
        num_heavy_atoms = mm.GetNumHeavyAtoms() if mm else float("inf")
        if num_heavy_atoms <= max_heavy_atoms:
            env, cores = get_std_context_core_permutations(
                contexts, core, radius, keep_stereo
            )
            if env and cores:
                for c in cores:
                    sma, dist2 = calc_mp(env, c)
                    results.append((env, c, num_heavy_atoms, sma, dist2))

    return results


def __fragment_mol_heavy_atoms(df, max_heavy_atoms, radius, keep_stereo):
    results = []
    for smi, smi_id in df[["smiles", "smi_id"]].values:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            sys.stderr.write(f"Can't generate mol for: {smi}\n")
            continue

        # heavy atoms
        frags = rdMMPA.FragmentMol(
            mol,
            pattern="[!#1]!@!=!#[!#1]",
            maxCuts=4,
            resultsAsMols=False,
            maxCutBonds=30,
        )
        frags += rdMMPA.FragmentMol(
            mol,
            pattern="[!#1]!@!=!#[!#1]",
            maxCuts=3,
            resultsAsMols=False,
            maxCutBonds=30,
        )
        frags = set(frags)
        for core, chains in frags:
            env_results = frag_to_env(
                smi, core, chains, max_heavy_atoms, radius, keep_stereo
            )
            for env, cores, num_heavy_atoms, sma, dist2 in env_results:
                results.append(
                    (smi, smi_id, core, chains, env, cores, num_heavy_atoms, sma, dist2)
                )
    return results


def __fragment_mol_hydrogen(df, max_heavy_atoms, radius, keep_stereo):
    results = []
    for smi, smi_id in df[["smiles", "smi_id"]].values:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            sys.stderr.write(f"Can't generate mol for: {smi}\n")
            continue

        # hydrogen splitting
        mol = Chem.AddHs(mol)
        n = mol.GetNumAtoms() - mol.GetNumHeavyAtoms()
        if n < 60:
            frags = rdMMPA.FragmentMol(
                mol,
                pattern="[#1]!@!=!#[!#1]",
                maxCuts=1,
                resultsAsMols=False,
                maxCutBonds=100,
            )
            for core, chains in frags:
                env_results = frag_to_env(
                    smi, core, chains, max_heavy_atoms, radius, keep_stereo
                )
                for env, cores, num_heavy_atoms, sma, dist2 in env_results:
                    results.append(
                        (
                            smi,
                            smi_id,
                            core,
                            chains,
                            env,
                            cores,
                            num_heavy_atoms,
                            sma,
                            dist2,
                        )
                    )
    return results


def batch_insert_db(data, db_manager, radius):
    envs = set()
    fragments = {}
    env_fragment_combo = {}
    for row in data:
        smi, smi_id, core, chains, env, core_smi, num_heavy_atoms, core_sma, dist2 = row
        envs.add(env)
        fragments.update({core_smi: num_heavy_atoms})
        if (env, core_smi) in env_fragment_combo:
            env_fragment_combo[(env, core_smi)]["freq"] += 1
        else:
            env_fragment_combo[(env, core_smi)] = {
                "core_sma": core_sma,
                "dist2": dist2,
                "freq": 1,
            }

    db_manager.insert(list(envs), fragments, env_fragment_combo, radius)


def upload_to_db(q, db_manager, radius, total_chunks):
    with tqdm(total=total_chunks, desc="Uploading to database") as pbar:
        while True:
            data = q.get()
            if data is None:
                break
            batch_insert_db(data, db_manager, radius)
            pbar.update(1)


def fragment_mols(args):
    # remove duplicated smiles
    args.input_file = preprocess_input_file(args.input_file)
    # create db manager
    db_manager = create_db_manager(
        args.db_type, args.db_path, args.ini_file, args.reset_db
    )
    total_rows = count_total_rows(args.input_file)
    print(f"Start import {total_rows} rows to {args.db_type}")
    total_chunks = total_rows // args.chunk_size + 1
    if args.mode == 1:
        total_chunks *= 2
    q = Queue()
    t = Thread(target=upload_to_db, args=(q, db_manager, args.radius, total_chunks))
    t.start()
    if args.debug:
        with open(args.out, "w") as f_output:
            # create csv writer
            csv_writer = csv.writer(f_output, delimiter=args.sep_out)
            # write header
            csv_writer.writerow(
                [
                    "smiles",
                    "smi_id",
                    "core",
                    "chains",
                    "env",
                    "cores",
                    "num_heavy_atoms",
                    "sma",
                    "dist2",
                ]
            )
            if args.mode in [0, 1]:
                print(f"Start fragment heavy atoms")
                chunks = read_chunks(args.input_file, args.chunk_size, args.sep)
                with Pool(args.ncpu) as p:
                    __fragment_mol = partial(
                        __fragment_mol_heavy_atoms,
                        max_heavy_atoms=args.max_heavy_atoms,
                        radius=args.radius,
                        keep_stereo=args.keep_stereo,
                    )
                    frags_heavy_atoms = p.imap_unordered(__fragment_mol, chunks)
                    for frag in frags_heavy_atoms:
                        if frag is None:
                            continue
                        q.put(frag)
                        csv_writer.writerows(frag)
            if args.mode in [1, 2]:
                print(f"Start fragment hydrogen atoms")
                chunks = read_chunks(args.input_file, args.chunk_size, args.sep)
                with Pool(args.ncpu) as p:
                    __fragment_mol = partial(
                        __fragment_mol_hydrogen,
                        max_heavy_atoms=args.max_heavy_atoms,
                        radius=args.radius,
                        keep_stereo=args.keep_stereo,
                    )
                    frags_hydrogen = p.imap_unordered(__fragment_mol, chunks)
                    for frag in frags_hydrogen:
                        if frag is None:
                            continue
                        q.put(frag)
                        csv_writer.writerows(frag)
    else:
        if args.mode in [0, 1]:
            print(f"Start fragment heavy atoms")
            chunks = read_chunks(args.input_file, args.chunk_size, args.sep)
            __fragment_mol = partial(
                __fragment_mol_heavy_atoms,
                max_heavy_atoms=args.max_heavy_atoms,
                radius=args.radius,
                keep_stereo=args.keep_stereo,
            )
            with Pool(args.ncpu) as p:
                frags_heavy_atoms = p.imap_unordered(__fragment_mol, chunks)
                for frag in frags_heavy_atoms:
                    if frag is None:
                        continue
                    q.put(frag)
        if args.mode in [1, 2]:
            print(f"Start fragment hydrogen atoms")
            chunks = read_chunks(args.input_file, args.chunk_size, args.sep)
            __fragment_mol = partial(
                __fragment_mol_hydrogen,
                max_heavy_atoms=args.max_heavy_atoms,
                radius=args.radius,
                keep_stereo=args.keep_stereo,
            )
            with Pool(args.ncpu) as p:
                frags_hydrogen = p.imap_unordered(__fragment_mol, chunks)
                for frag in frags_hydrogen:
                    if frag is None:
                        continue
                    q.put(frag)

    q.put(None)
    t.join()
    print(f"finished uploading {args.input_file} to {args.db_type} database")


def create_db(args):
    fragment_mols(args)


if __name__ == "__main__":
    args = parse_args()
    create_db(args)
