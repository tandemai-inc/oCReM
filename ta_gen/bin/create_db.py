#!/usr/bin/env python
# -*- coding:utf-8 -*-

import argparse
import pandas as pd
import os
from multiprocessing import Pool, cpu_count
from rdkit import Chem
import sys
from rdkit.Chem import rdMMPA
import csv

def schema_parser():
    parser = argparse.ArgumentParser(description="Create database")
    parser.add_argument("-i", "--input_file", type=str, required=True, help="Input compounds")
    parser.add_argument('-o', '--out', metavar='output.txt', required=False, default='output.txt', help='output file of fragmented molecules.')
    parser.add_argument( '--debug', action='store_true', default=False, help='print debug message.')
    parser.add_argument("--chunk_size", type=int, default=100, help="Chunk size")
    parser.add_argument('-c', '--ncpu', metavar='NUMBER', required=False, default=1, help='number of cpus used for computation. Default: 1.')
    parser.add_argument('-s', '--sep', metavar='STRING', required=False, default=",",
                        help='separator in input file. Default: Tab.')
    parser.add_argument('-d', '--sep_out', metavar='STRING', required=False, default=',',
                        help='separator in the output file. Default: comma')
    parser.add_argument('-m', '--mode', metavar='INTEGER', required=False, default=0,
                        choices=[0, 1, 2], type=int,
                        help='fragmentation mode: 0 - all atoms constitute a fragment, 1- heavy atoms only, '
                             '2 - hydrogen atoms only. Default: 0.')
    return parser


def parse_args():
    parser = schema_parser()
    args, _ = parser.parse_known_args()
    return args


def read_chunks(input_file, chunk_size, sep):
    input_file = args.input_file
    _, ext = os.path.splitext(input_file)
    if ext == ".csv":
        chunks = pd.read_csv(input_file, sep=args.sep, chunksize=args.chunk_size)
    elif ext == ".smi":
        chunks = pd.read_csv(input_file, header=None, names=["smiles", "smi_id"], sep=args.sep,
                             chunksize=args.chunk_size)
    else:
        chunks = pd.read_csv(input_file, header=None, names=["smiles", "smi_id"], sep=args.sep,
                             chunksize=args.chunk_size)
    return chunks


def __fragment_mol_heavy_atoms(df):
    results = []
    for smi, smi_id in df[["smiles", "smi_id"]].values:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            sys.stderr.write(f"Can't generate mol for: {smi}\n")
            continue

        # heavy atoms
        frags = rdMMPA.FragmentMol(mol, pattern="[!#1]!@!=!#[!#1]", maxCuts=4,
                                   resultsAsMols=False, maxCutBonds=30)
        frags += rdMMPA.FragmentMol(mol, pattern="[!#1]!@!=!#[!#1]", maxCuts=3,
                                    resultsAsMols=False, maxCutBonds=30)
        frags = set(frags)
        results.extend([(smi, smi_id, core, chains) for core, chains in frags])
    return results


def __fragment_mol_hydrogen(df):
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
            frags = rdMMPA.FragmentMol(mol, pattern="[#1]!@!=!#[!#1]", maxCuts=1,
                                       resultsAsMols=False, maxCutBonds=100)
            results.extend([(smi, smi_id, core, chains) for core, chains in frags])
    return results

def fragment_mols(args):
    if args.debug:
        with open(args.out, 'w') as f_output:
            # create csv writer
            csv_writer = csv.writer(f_output, delimiter=args.sep_out)
            # write header
            csv_writer.writerow(["smiles", "smi_id", "core", "chains"])
            if args.mode in [0, 1]:
                chunks = read_chunks(args.input_file, args.chunk_size, args.sep)
                with Pool(args.ncpu) as p:
                    frags_heavy_atoms = p.imap_unordered(__fragment_mol_heavy_atoms, chunks)
                    for frag in frags_heavy_atoms:
                        if frag is None:
                            continue
                        csv_writer.writerows(frag)
            if args.mode in [1, 2]:
                chunks = read_chunks(args.input_file, args.chunk_size, args.sep)
                with Pool(args.ncpu) as p:
                    frags_hydrogen = p.imap_unordered(__fragment_mol_hydrogen, chunks)
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
