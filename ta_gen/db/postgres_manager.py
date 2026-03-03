#!/usr/bin/env python
# -*- coding:utf-8 -*-


import configparser
import copy

import psycopg2
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT

from ta_gen.db.db_manager import DBManager


def load_ini(ini_file):
    config = configparser.ConfigParser()
    config.read(ini_file)
    return dict(config["database"])


class PostGresManager(DBManager):

    def __init__(self, ini_file):
        super().__init__()
        self.conn_params = load_ini(ini_file)

    def create_db(self, ini_file):
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
                    id BIGINT PRIMARY KEY,
                    name TEXT,
                    radis SMALLINT
                )
            """)

            # Create fragment table with correct structure
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS fragment (
                    id BIGINT PRIMARY KEY,
                    core_smi TEXT,
                    core_num_atoms BIGINT,
                    dist2 BIGINT
                )
            """)

            # Create env_fragment table with correct structure
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS env_fragment (
                    env_id BIGINT REFERENCES env(id),
                    fragment_id BIGINT REFERENCES fragment(id),
                    frequency BIGINT
                )
            """)

            conn.commit()
            cursor.close()
            conn.close()
        except Exception as e:
            print(f"Error creating tables: {e}")
            raise Exception(f"Error creating tables: {e}")


    def connect_db(self):
        self.conn = psycopg2.connect(**self.conn_params)
        self.cursor = self.conn.cursor()

    def insert_new_fragment(self, fragments):
        core_smis = list(fragments.keys())

