#!/usr/bin/env python
# -*- coding:utf-8 -*-


import sqlite3
import traceback

from ta_gen.db.db_manager import DBManager


class SqliteManager(DBManager):

    def __init__(self, db_path="ocrem.db", reset_db=False):
        super().__init__()
        self.db_path = db_path
        if reset_db:
            self.clear_db(db_path)
        self.create_db(db_path)
        # create connection
        self.connect_db()

    def create_db(self, db_path):
        """create database"""
        # connect to database
        conn = sqlite3.connect(db_path)
        # support foreign key
        conn.execute("PRAGMA foreign_keys = ON")
        cursor = conn.cursor()

        try:
            cursor.execute("""
                    CREATE TABLE IF NOT EXISTS env (
                        id INTEGER PRIMARY KEY,  
                        name TEXT UNIQUE,
                        radius INTEGER           
                    )
                """)

            cursor.execute("""
                    CREATE TABLE IF NOT EXISTS fragment (
                        id INTEGER PRIMARY KEY,
                        core_smi TEXT UNIQUE,
                        core_num_atoms INTEGER  
                    )
                """)

            cursor.execute("""
                    CREATE TABLE IF NOT EXISTS env_fragment (
                        env_id INTEGER,
                        fragment_id INTEGER,
                        core_sma TEXT,
                        dist2 INTEGER,
                        frequency INTEGER,
                        PRIMARY KEY (env_id, fragment_id),
                        FOREIGN KEY (env_id) REFERENCES env(id) ON DELETE CASCADE,
                        FOREIGN KEY (fragment_id) REFERENCES fragment(id) ON DELETE CASCADE
                    )
                """)

            cursor.execute(
                "CREATE INDEX IF NOT EXISTS idx_env_fragment_env_id ON env_fragment(env_id)"
            )
            cursor.execute(
                "CREATE INDEX IF NOT EXISTS idx_env_fragment_fragment_id ON env_fragment(fragment_id)"
            )

            conn.commit()
            print("all tables created successfully (or already exist)")

        except sqlite3.Error as e:
            print(f"database error: {e}")
            conn.rollback()
        finally:
            conn.close()

    def clear_db(self, db_path):
        """clear database"""
        print(f"clearing database {db_path}")
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        try:
            cursor.execute("DROP TABLE IF EXISTS env_fragment")
            cursor.execute("DROP TABLE IF EXISTS fragment")
            cursor.execute("DROP TABLE IF EXISTS env")

            conn.commit()
            print("all tables cleared successfully")
        except sqlite3.Error as e:
            print(f"database error: {e}")
            conn.rollback()
        finally:
            conn.close()

    def connect_db(self):
        """connect database"""
        self.conn = sqlite3.connect(self.db_path, timeout=10.0)
        # support foreign key
        self.conn.execute("PRAGMA foreign_keys = ON")
        # turn on synchronous mode
        self.conn.execute("PRAGMA journal_mode=WAL")
        self.conn.execute("PRAGMA synchronous=NORMAL")
        self.cursor = self.conn.cursor()

    def insert_new_env(self, envs, radius):
        placeholders = ",".join(["?"] * len(envs))
        self.cursor.execute(
            f"SELECT name, id FROM env WHERE name IN ({placeholders})", envs
        )
        env_map = {row[0]: row[1] for row in self.cursor.fetchall()}
        missing_envs = [name for name in envs if name not in env_map]

        if missing_envs:
            # insert new env
            self.cursor.executemany(
                "INSERT OR IGNORE INTO env (name, radius) VALUES (?, ?)",
                [(name, radius) for name in missing_envs],
            )
            # get new ids
            self.cursor.execute(
                f"SELECT name, id FROM env WHERE name IN ({','.join(['?'] * len(missing_envs))})",
                missing_envs,
            )
            for row in self.cursor.fetchall():
                env_map[row[0]] = row[1]

        return env_map

    def insert_new_fragment(self, fragments):
        core_smis = list(fragments.keys())
        placeholders = ",".join(["?"] * len(core_smis))
        self.cursor.execute(
            f"SELECT core_smi, id FROM fragment WHERE core_smi IN ({placeholders})",
            core_smis,
        )
        fragment_map = {row[0]: row[1] for row in self.cursor.fetchall()}
        missing_fragments = [name for name in core_smis if name not in fragment_map]

        if missing_fragments:
            # insert new fragment
            self.cursor.executemany(
                "INSERT OR IGNORE INTO fragment (core_smi, core_num_atoms) VALUES (?, ?)",
                [(name, fragments.get(name)) for name in missing_fragments],
            )
            # get new ids
            self.cursor.execute(
                f"SELECT core_smi, id FROM fragment WHERE core_smi IN ({','.join(['?'] * len(missing_fragments))})",
                missing_fragments,
            )
            for row in self.cursor.fetchall():
                fragment_map[row[0]] = row[1]

        return fragment_map

    def insert_env_fragment(self, env_fragment_combo, fragment_ids, env_ids):
        upsert_data = [
            (
                env_ids[env],
                fragment_ids[core_smi],
                attr["core_sma"],
                attr["dist2"],
                attr["freq"],
            )
            for (env, core_smi), attr in env_fragment_combo.items()
        ]
        upsert_sql = """
            INSERT INTO env_fragment (env_id, fragment_id, core_sma, dist2, frequency)
            VALUES (?, ?, ?, ?, ?)
            ON CONFLICT(env_id, fragment_id) DO UPDATE SET
            frequency = frequency + excluded.frequency
        """
        self.cursor.executemany(upsert_sql, upsert_data)

    def insert(self, envs, fragments, env_fragment_combo, radius):
        try:
            env_ids = self.insert_new_env(envs, radius)
            fragment_ids = self.insert_new_fragment(fragments)
            self.insert_env_fragment(env_fragment_combo, fragment_ids, env_ids)
            self.conn.commit()
        except Exception as e:
            traceback.print_exc()
            self.conn.rollback()  # rollback

    def execute(self, sql):
        self.connect_db()
        try:
            self.cursor.execute(sql)
            return self.cursor.fetchall() or []
        except Exception as e:
            traceback.print_exc()
            self.conn.rollback()  # rollback
            return []

    def close(self):
        self.cursor.close()
        self.conn.close()
