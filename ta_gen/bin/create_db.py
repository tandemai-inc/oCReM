#!/usr/bin/env python
# -*- coding:utf-8 -*-

import argparse
import csv
import os
import sys
from functools import partial
from itertools import permutations
from multiprocessing import Pool, cpu_count

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMMPA

from ta_gen.utils.mol_context import get_std_context_core_permutations


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
        default=1,
        help="radius of molecular context (in bonds) which will be taken into account. Default: 1.",
    )
    parser.add_argument(
        "--keep_stereo",
        action="store_true",
        default=False,
        help="set this flag if you want to keep stereo in context and core parts.",
    )
    return parser


def parse_args():
    parser = schema_parser()
    args, _ = parser.parse_known_args()
    return args


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


def frag_to_env(smi, core, contexts, max_heavy_atoms, radius, keep_stereo):
    if not core and not contexts:
        return None

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
                        return env, cores[0], num_heavy_atoms
        else:
            sys.stderr.write(
                f"more than two fragments in context ({contexts}) where core is empty for smiles: {smi}"
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
                return env, cores[0], num_heavy_atoms

    return None


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
            if env_results is not None:
                env, cores, num_heavy_atoms = env_results
                results.append((smi, smi_id, core, chains, env, cores, num_heavy_atoms))
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
                if env_results is not None:
                    env, cores, num_heavy_atoms = env_results
                    results.append(
                        (smi, smi_id, core, chains, env, cores, num_heavy_atoms)
                    )
    return results


def fragment_mols(args):
    if args.debug:
        with open(args.out, "w") as f_output:
            # create csv writer
            csv_writer = csv.writer(f_output, delimiter=args.sep_out)
            # write header
            csv_writer.writerow(["smiles", "smi_id", "core", "chains"])
            if args.mode in [0, 1]:
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
                        csv_writer.writerows(frag)
            if args.mode in [1, 2]:
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
                        csv_writer.writerows(frag)
    else:
        if args.mode in [0, 1]:
            chunks = read_chunks(args.input_file, args.chunk_size, args.sep)
            with Pool(args.ncpu) as p:
                frags_heavy_atoms = p.imap_unordered(__fragment_mol_heavy_atoms, chunks)
        if args.mode in [1, 2]:
            chunks = read_chunks(args.input_file, args.chunk_size, args.sep)
            with Pool(args.ncpu) as p:
                frags_hydrogen = p.imap_unordered(__fragment_mol_hydrogen, chunks)


def create_db(args):
    fragment_mols(args)


if __name__ == "__main__":
    args = parse_args()
    create_db(args)
