import argparse
import copy
import os

import psycopg2
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT
from tqdm import tqdm

import configparser

def load_ini(ini_file):
    config = configparser.ConfigParser()
    config.read(ini_file)
    return dict(config["database"])


def create_db(conn_params):
    # check if new db already exists
    try:
        default_db_config = copy.deepcopy(conn_params)
        default_db_config["user"] = "postgres"
        conn = psycopg2.connect(**default_db_config)
        conn.set_isolation_level(ISOLATION_LEVEL_AUTOCOMMIT)
        cursor = conn.cursor()

        # 检查数据库是否存在
        cursor.execute("SELECT 1 FROM pg_database WHERE datname = 'gen_chembl_tmp'")
        exists = cursor.fetchone()

        if not exists:
            cursor.execute("CREATE DATABASE gen_chembl_tmp")
            print("Database created")
        else:
            print("Database already exists")

        cursor.close()
        conn.close()
    except Exception as e:
        print(f"Error creating database: {e}")
        raise Exception(f"Error creating database: {e}")

    # create tables
    try:
        conn = psycopg2.connect(**conn_params)
        cursor = conn.cursor()

        # Drop tables if they exist
        cursor.execute("DROP TABLE IF EXISTS env_fragment")
        cursor.execute("DROP TABLE IF EXISTS env")
        cursor.execute("DROP TABLE IF EXISTS fragment")

        # Create env table with correct structure
        cursor.execute(
            """
            CREATE TABLE env (
                id BIGINT,
                name TEXT,
                chembl BIGINT,
                enamine BIGINT,
                pubchem BIGINT,
                surechem BIGINT,
                zinc BIGINT
            )
        """
        )

        # Create fragment table with correct structure
        cursor.execute(
            """
            CREATE TABLE fragment (
                id BIGINT,
                core_smi TEXT,
                core_num_atoms BIGINT,
                core_sma TEXT,
                dist2 BIGINT
            )
        """
        )

        # Create env_fragment table with correct structure
        cursor.execute(
            """
            CREATE TABLE env_fragment (
                env_id BIGINT,
                fragment_id BIGINT
            )
        """
        )

        conn.commit()
        cursor.close()
        conn.close()
    except Exception as e:
        print(f"Error creating tables: {e}")
        raise Exception(f"Error creating tables: {e}")


def import_data(conn_params, data_folder, delimiter=",", null_string=""):
    for f in os.listdir(data_folder):
        if f.endswith(".csv"):
            csv_path = os.path.join(data_folder, f)
            table_name = f.split(".")[0]
            print(f"importing {table_name} from {csv_path}")
            conn = psycopg2.connect(**conn_params)
            cursor = conn.cursor()

            # Get total rows for progress bar (excluding header)
            with open(csv_path, "r", encoding="utf-8") as f:
                total_rows = sum(1 for _ in f) - 1

            print(f"Importing {table_name} ({total_rows:,} rows)...")

            # Use COPY command with progress tracking
            with open(csv_path, "r", encoding="utf-8") as f:
                # Skip header
                next(f)

                # Create a temporary table-like structure for copy
                # Use copy_expert for maximum flexibility
                copy_sql = f"""
                    COPY {table_name} FROM STDIN 
                    WITH (
                        FORMAT CSV,
                        DELIMITER '{delimiter}',
                        NULL '{null_string}',
                        ENCODING 'UTF8'
                    )
                """

                # Execute copy with progress callback
                with tqdm(
                    total=total_rows, desc=f"COPY {table_name}", unit="rows"
                ) as pbar:
                    # Use a wrapper to track progress
                    class ProgressReader:
                        def __init__(self, file_obj, progress_bar):
                            self.file_obj = file_obj
                            self.progress_bar = progress_bar
                            self.line_count = 0

                        def read(self, size=-1):
                            data = self.file_obj.read(size)
                            if size == -1 or "\n" in data:
                                self.line_count += data.count("\n")
                                self.progress_bar.update(data.count("\n"))
                            return data

                    progress_reader = ProgressReader(f, pbar)
                    cursor.copy_expert(copy_sql, progress_reader)

            conn.commit()
            cursor.close()
            conn.close()


def restore_indexes(conn_params):
    INDEX_SQL = [
        "CREATE UNIQUE INDEX IF NOT EXISTS env_pkey ON env USING btree (id)",
        "CREATE UNIQUE INDEX IF NOT EXISTS env_name_key ON env USING btree (name)",
        "CREATE UNIQUE INDEX IF NOT EXISTS fragment_pkey ON fragment USING btree (id)",
        "CREATE UNIQUE INDEX IF NOT EXISTS fragment_core_smi_key ON fragment USING btree (core_smi)",
        "CREATE UNIQUE INDEX IF NOT EXISTS env_fragment_pkey ON env_fragment USING btree (env_id, fragment_id)",
    ]

    conn = psycopg2.connect(**conn_params)
    cursor = conn.cursor()

    for sql in INDEX_SQL:
        try:
            cursor.execute(sql)
            print(f"Created index: {sql[:50]}...")
        except Exception as e:
            print(f"Error: {e}")

    conn.commit()
    cursor.close()
    conn.close()
    print("Index restoration completed")

    print("\nIndex restoration completed successfully!")


def schema_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--ini_file", type=str)
    parser.add_argument("-d", "--data_folder", type=str)
    return parser.parse_args()


def parse_args():
    parser = schema_parser()
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    conn_params = load_ini(args.ini_file)
    data_folder = args.data_folder
    create_db(conn_params)
    import_data(conn_params, data_folder)
    restore_indexes(conn_params)
