#!/usr/bin/env python
# -*- coding:utf-8 -*-


import configparser
import copy
import traceback

import psycopg2
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT

from ta_gen.db.db_manager import DBManager


def load_ini(ini_file):
    config = configparser.ConfigParser()
    config.read(ini_file)
    return dict(config["database"])


class PostGresManager(DBManager):

    def __init__(self, ini_file, reset_db=False):
        super().__init__()
        self.conn_params = load_ini(ini_file)
        db_exists = self.db_exist()
        if db_exists and reset_db:
            self.clear_db()
        self.create_db()

    def db_exist(self):
        # check if db already exists
        try:
            default_db_config = copy.deepcopy(self.conn_params)
            db_name = self.conn_params["database"]
            default_db_config["user"] = "postgres"
            default_db_config["database"] = "postgres"

            conn = psycopg2.connect(**default_db_config)
            conn.set_isolation_level(ISOLATION_LEVEL_AUTOCOMMIT)
            cursor = conn.cursor()

            # check if db already exists
            cursor.execute(f"SELECT 1 FROM pg_database WHERE datname = '{db_name}'")
            exists = cursor.fetchone()
            cursor.close()
            conn.close()

            if not exists:
                return False
            else:
                return True
        except Exception as e:
            raise Exception(f"Error checking database existence: {e}")

    def create_db(self):
        # create database
        try:
            default_db_config = copy.deepcopy(self.conn_params)
            db_name = self.conn_params["database"]
            default_db_config["user"] = "postgres"
            default_db_config["database"] = "postgres"

            conn = psycopg2.connect(**default_db_config)
            conn.set_isolation_level(ISOLATION_LEVEL_AUTOCOMMIT)
            cursor = conn.cursor()

            # check if db already exists
            cursor.execute(f"SELECT 1 FROM pg_database WHERE datname = '{db_name}'")
            exists = cursor.fetchone()

            if not exists:
                cursor.execute(f'CREATE DATABASE "{db_name}"')
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
            conn = psycopg2.connect(**self.conn_params)
            cursor = conn.cursor()

            # Create env table with correct structure
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS env (
                    id BIGSERIAL PRIMARY KEY,
                    name TEXT UNIQUE NOT NULL,
                    radis SMALLINT
                )
            """)

            # Create fragment table with correct structure
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS fragment (
                    id BIGSERIAL PRIMARY KEY,
                    core_smi TEXT UNIQUE NOT NULL,
                    core_num_atoms INTEGER,
                    core_sma TEXT,
                    dist2 BIGINT
                )
            """)

            # Create env_fragment table with correct structure
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS env_fragment (
                    env_id BIGINT REFERENCES env(id),
                    fragment_id BIGINT REFERENCES fragment(id),
                    frequency BIGINT,
                    PRIMARY KEY (env_id, fragment_id)
                )
            """)

            cursor.execute(
                "CREATE INDEX IF NOT EXISTS idx_env_fragment_env_id ON env_fragment(env_id)"
            )
            cursor.execute(
                "CREATE INDEX IF NOT EXISTS idx_env_fragment_fragment_id ON env_fragment(fragment_id)"
            )

            conn.commit()
            cursor.close()
            conn.close()
        except Exception as e:
            print(f"Error creating tables: {e}")
            raise Exception(f"Error creating tables: {e}")

    def clear_db(self):
        print(f"clearing database {self.conn_params['database']}")
        try:
            conn = psycopg2.connect(**self.conn_params)
            cursor = conn.cursor()

            cursor.execute("TRUNCATE TABLE env_fragment CASCADE")
            cursor.execute("TRUNCATE TABLE fragment CASCADE")
            cursor.execute("TRUNCATE TABLE env CASCADE")

            conn.commit()
            cursor.close()
            conn.close()
        except Exception as e:
            print(f"Error clearing tables: {e}")
            raise Exception(f"Error clearing tables: {e}")

    def connect_db(self):
        self.conn = psycopg2.connect(**self.conn_params)
        self.cursor = self.conn.cursor()

    def insert_new_env(self, envs, radius):
        placeholders = ",".join(["%s"] * len(envs))
        self.cursor.execute(
            f"SELECT name, id FROM env WHERE name IN ({placeholders})", envs
        )
        env_map = {row[0]: row[1] for row in self.cursor.fetchall()}
        missing = [name for name in envs if name not in env_map]

        if missing:
            self.cursor.executemany(
                "INSERT INTO env (name, radis) VALUES (%s, %s)",
                [(name, radius) for name in missing],
            )
            self.cursor.execute(
                f"SELECT name, id FROM env WHERE name IN ({','.join(['%s'] * len(missing))})",
                missing,
            )
            for row in self.cursor.fetchall():
                env_map[row[0]] = row[1]

        return env_map

    def insert_new_fragment(self, fragments):
        core_smis = list(fragments.keys())
        placeholders = ",".join(["%s"] * len(core_smis))
        self.cursor.execute(
            f"SELECT core_smi, id FROM fragment WHERE core_smi IN ({placeholders})",
            core_smis,
        )
        fragment_map = {row[0]: row[1] for row in self.cursor.fetchall()}
        missing = [name for name in core_smis if name not in fragment_map]
        if missing:
            data = [(name, fragments[name][0], fragments[name][1]) for name in missing]
            self.cursor.executemany(
                "INSERT INTO fragment (core_smi, core_num_atoms, core_sma, dist2) VALUES (%s, %s, %s, %s)",
                data,
            )
            self.cursor.execute(
                f"SELECT core_smi, id FROM fragment WHERE core_smi IN ({','.join(['%s'] * len(missing))})",
                missing,
            )
            for row in self.cursor.fetchall():
                fragment_map[row[0]] = row[1]

        return fragment_map

    def insert_env_fragment(self, env_fragment_counter, fragment_ids, env_ids):
        upsert_data = [
            (env_ids[env], fragment_ids[core_smi], count)
            for (env, core_smi), count in env_fragment_counter.items()
        ]

        upsert_sql = """
            INSERT INTO env_fragment (env_id, fragment_id, frequency)
            VALUES (%s, %s, %s)
            ON CONFLICT (env_id, fragment_id) DO UPDATE SET
            frequency = env_fragment.frequency + EXCLUDED.frequency
        """
        self.cursor.executemany(upsert_sql, upsert_data)

    def insert(self, envs, fragments, env_fragment_counter, radius):
        self.connect_db()
        try:
            with self.conn:
                env_ids = self.insert_new_env(envs, radius)
                fragment_ids = self.insert_new_fragment(fragments)
                self.insert_env_fragment(env_fragment_counter, fragment_ids, env_ids)
                self.conn.commit()
        except Exception as e:
            traceback.print_exc()
        finally:
            self.close()

    def close(self):
        """close cursor and connection"""
        self.cursor.close()
        self.conn.close()
