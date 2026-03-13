#!/usr/bin/env python
# -*- coding:utf-8 -*-

import argparse
import os
import subprocess
from multiprocessing import Process, Queue, cpu_count

import pandas as pd
from tqdm import tqdm

from ta_gen.db import create_db_manager


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


def prepare_args(args):
    if args.db_type == "postgres":
        assert os.path.exists(
            args.ini_file
        ), f"ini_file {args.ini_file} does not exist."


def parse_args():
    parser = schema_parser()
    args, _ = parser.parse_known_args()
    prepare_args(args)
    return args


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
        chunks = pd.read_csv(
            input_file, sep=sep, chunksize=chunk_size, on_bad_lines="skip"
        )
    else:
        chunks = pd.read_csv(
            input_file,
            sep=sep,
            chunksize=chunk_size,
            on_bad_lines="skip",
            names=[
                "smi",
                "smi_id",
                "core",
                "chains",
                "env",
                "core_smi",
                "num_heavy_atoms",
                "core_sma",
                "dist2",
            ],
        )

    return chunks


def batch_insert_db(data, db_manager, radius):
    envs = set()
    fragments = {}
    env_fragment_combo = {}
    for frag in data.values:
        smi, smi_id, core, chains, env, core_smi, num_heavy_atoms, core_sma, dist2 = (
            frag
        )
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


def upload_to_db(q, db_type, db_path, ini_file, reset_db, radius, total_chunks):
    db_manager = create_db_manager(db_type, db_path, ini_file, reset_db)
    with tqdm(total=total_chunks, desc="Uploading to database") as pbar:
        while True:
            data = q.get()
            if data is None:
                break
            batch_insert_db(data, db_manager, radius)
            pbar.update(1)
        db_manager.close()


class DBManager(object):

    def __init__(self, use_db, args):
        self.use_db = use_db
        self.args = args

        # create thread to upload
        self.q = Queue()
        self.upload_thread = Process(
            target=upload_to_db,
            args=(
                self.q,
                args.db_type,
                args.db_path,
                args.ini_file,
                args.reset_db,
                args.radius,
                args.total_chunks,
            ),
        )
        self.upload_thread.start()

        self.update_queue = self.__update_queue
        self.join = self.upload_thread.join

    def __update_queue(self, data):
        self.q.put(data)


def result_to_db(args):
    total_rows = count_total_rows(args.input_file)
    print(f"Start import {total_rows} rows from {args.input_file} to {args.db_type}")
    print(f"Applied cpus: {args.ncpu}. Total cpus: {cpu_count()}")
    args.ncpu = min(args.ncpu, cpu_count())
    args.total_chunks = total_rows // args.chunk_size + 1
    # create db manager
    db_manager = DBManager(True, args)
    print(f"Start save results to db")
    chunks = read_chunks(args.input_file, args.chunk_size, args.sep)
    for chunk in chunks:
        db_manager.update_queue(chunk)

    db_manager.update_queue(None)
    db_manager.join()
    print(f"finished uploading {args.input_file} to {args.db_type} database")


if __name__ == "__main__":
    args = parse_args()
    result_to_db(args)
