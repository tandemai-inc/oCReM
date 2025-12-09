#!/usr/bin/env python
# -*- coding:utf-8 -*-

import argparse
import os
import sqlite3

import pandas as pd
from tqdm import tqdm

EXPORT_DIR = "/nfs/workspace/yinhuiyi/codedir/CrEM-Database/tmp/export_20251205_143245"
SQLITE_DB = "/nfs/workspace/yinhuiyi/codedir/CrEM-Database/tmp/gen_chembl_tmp.db"


def create_sqlite_db(db_file):
    """create sqlite db"""
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()

    # create env table
    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS env (
            id INTEGER,
            name TEXT,
            chembl INTEGER,
            enamine INTEGER,
            pubchem INTEGER,
            surechem INTEGER,
            zinc INTEGER
        )
    """
    )

    # create fragment table
    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS fragment (
            id INTEGER,
            core_smi TEXT,
            core_num_atoms INTEGER,
            core_sma TEXT,
            dist2 INTEGER
        )
    """
    )

    # create env_fragment table
    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS env_fragment (
            env_id INTEGER,
            fragment_id INTEGER
        )
    """
    )

    conn.commit()
    conn.close()
    print(f"SQLite database created: {db_file}")


def import_to_sqlite(db_file, data_folder):
    """import csv data to sqlite db"""
    conn = sqlite3.connect(db_file)

    csv_files = [f for f in os.listdir(data_folder) if f.endswith(".csv")]

    for csv_file in csv_files:
        table_name = csv_file.split(".")[0]
        csv_path = os.path.join(data_folder, csv_file)

        # get total rows in csv file
        with open(csv_path, "r", encoding="utf-8") as f:
            total_rows = sum(1 for _ in f) - 1

        print(f"Importing {csv_file} to {table_name} ({total_rows:,} rows)...")

        # read csv file in chunks to avoid memory issues
        chunk_size = 100000
        chunks = pd.read_csv(csv_path, chunksize=chunk_size)

        with tqdm(
            total=total_rows, desc=f"Importing {table_name}", unit="rows"
        ) as pbar:
            for i, chunk in enumerate(chunks):
                # create table if not exists, append if exists
                if_exists = "append" if i > 0 else "replace"
                chunk.to_sql(table_name, conn, if_exists=if_exists, index=False)
                pbar.update(len(chunk))

        print(f"Imported {csv_file}")

    conn.close()
    print("All files imported to SQLite")


def create_sqlite_indexes(db_file):
    """create indexes in sqlite db"""
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()

    index_definitions = [
        ("CREATE UNIQUE INDEX IF NOT EXISTS env_pkey ON env (id)"),
        ("CREATE UNIQUE INDEX IF NOT EXISTS env_name_key ON env (name)"),
        ("CREATE UNIQUE INDEX IF NOT EXISTS fragment_pkey ON fragment (id)"),
        (
            "CREATE UNIQUE INDEX IF NOT EXISTS fragment_core_smi_key ON fragment (core_smi)"
        ),
        (
            "CREATE UNIQUE INDEX IF NOT EXISTS env_fragment_pkey ON env_fragment (env_id, fragment_id)"
        ),
    ]

    for sql in index_definitions:
        try:
            cursor.execute(sql)
            index_name = sql.split("INDEX IF NOT EXISTS ")[1].split(" ON")[0]
            print(f"Created index: {index_name}")
        except Exception as e:
            print(f"Error creating index: {e}")

    conn.commit()
    conn.close()
    print("SQLite indexes created")


def optimize_sqlite(db_file):
    """optimize sqlite db performance"""
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()

    # enable foreign key constraints
    cursor.execute("PRAGMA foreign_keys = ON")

    # set journal mode to WAL (improves concurrent performance)
    cursor.execute("PRAGMA journal_mode = WAL")

    # set synchronous mode to NORMAL (better performance)
    cursor.execute("PRAGMA synchronous = NORMAL")

    # set cache size
    cursor.execute("PRAGMA cache_size = -10000")  # 10MB cache

    # execute VACUUM to optimize database file
    cursor.execute("VACUUM")

    conn.commit()
    conn.close()
    print("SQLite database optimized")


def verify_import(db_file):
    """verify import results"""
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()

    tables = ["env", "fragment", "env_fragment"]

    print("\nVerifying import results:")
    print("-" * 40)

    for table in tables:
        cursor.execute(f"SELECT COUNT(*) FROM {table}")
        count = cursor.fetchone()[0]
        print(f"{table}: {count:,} rows")

    # check indexes
    print("\nChecking indexes:")
    cursor.execute("SELECT name FROM sqlite_master WHERE type='index'")
    indexes = cursor.fetchall()

    for idx in indexes:
        if idx[0].startswith("sqlite_autoindex"):  # skip auto-created indexes
            continue
        print(f"  - {idx[0]}")

    conn.close()



def schema_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-D", "--db_file", type=str, default="gen_chembl.db")
    parser.add_argument("-d", "--data_folder", type=str)
    return parser.parse_args()


def parse_args():
    parser = schema_parser()
    return parser.parse_args()


if __name__ == "__main__":
    print("=" * 60)
    print("IMPORTING TO SQLITE DATABASE")
    print("=" * 60)

    args = parse_args()
    db_file = args.db_file
    data_folder = args.data_folder

    # step 1: create sqlite database and tables
    create_sqlite_db(db_file)

    # step 2: import data
    import_to_sqlite(db_file, data_folder)

    # step 3: create indexes
    create_sqlite_indexes(db_file)

    # step 4: optimize database
    optimize_sqlite(db_file)

    # step 5: verify results
    verify_import(db_file)

    print("\n" + "=" * 60)
    print("SQLITE IMPORT COMPLETED SUCCESSFULLY!")
    print(f"Database file: {db_file}")
    print("=" * 60)
