#!/usr/bin/env python
# -*- coding:utf-8 -*-

import json
import os

import psycopg2
import psycopg2.extras


class PostgresImporter:
    def __init__(self, host, database, user, password, port=5432, export_dir=None):
        self.connection_params = {
            "host": host,
            "database": database,
            "user": user,
            "password": password,
            "port": port,
        }
        self.export_dir = export_dir or self.find_latest_export()

    def find_latest_export(self):
        """find latest export directory"""
        export_dirs = [
            d for d in os.listdir(".") if os.path.isdir(d) and d.startswith("export_")
        ]
        if not export_dirs:
            raise Exception("no export directory found")
        return max(export_dirs)  # sort by name, latest one is first

    def connect(self):
        """connect to database"""
        try:
            self.conn = psycopg2.connect(**self.connection_params)
            self.cur = self.conn.cursor()
            print("database connection successful")
        except Exception as e:
            print(f"connection failed: {e}")
            raise

    def execute_sql_file(self, sql_file):
        """execute sql file"""
        if not os.path.exists(sql_file):
            print(f"sql file not found: {sql_file}")
            return False

        with open(sql_file, "r", encoding="utf-8") as f:
            sql_content = f.read()

        try:
            # split sql statements by semicolon
            statements = [
                stmt.strip() for stmt in sql_content.split(";") if stmt.strip()
            ]

            for stmt in statements:
                if stmt:
                    self.cur.execute(stmt)

            self.conn.commit()
            print(f"successfully executed: {sql_file}")
            return True
        except Exception as e:
            self.conn.rollback()
            print(f"failed to execute {sql_file}: {e}")
            return False

    def import_data_from_csv(self, table_name):
        """import data from csv file"""
        csv_file = os.path.join(self.export_dir, "data", f"{table_name}.csv")

        if not os.path.exists(csv_file):
            print(f"{table_name}: csv file not found, skip")
            return True

        try:
            # disable all triggers on the table (to speed up import)
            self.cur.execute(f"ALTER TABLE {table_name} DISABLE TRIGGER ALL")

            # use COPY command to import (fastest method)
            with open(csv_file, "r", encoding="utf-8") as f:
                # skip header row
                next(f)
                self.cur.copy_expert(f"COPY {table_name} FROM STDIN WITH CSV", f)

            # enable all triggers on the table
            self.cur.execute(f"ALTER TABLE {table_name} ENABLE TRIGGER ALL")

            self.conn.commit()
            print(f"{table_name}: data import successful")
            return True
        except Exception as e:
            self.conn.rollback()
            print(f"{table_name}: copy failed, use insert instead: {e}")
            return self.import_with_insert(table_name, csv_file)

    def execute_sql_file(self, sql_file):
        """execute sql file"""
        if not os.path.exists(sql_file):
            print(f"sql file not found: {sql_file}")
            return False

        with open(sql_file, "r", encoding="utf-8") as f:
            sql_content = f.read()

        try:
            # 按分号分割SQL语句
            statements = [
                stmt.strip() for stmt in sql_content.split(";") if stmt.strip()
            ]

            for stmt in statements:
                if stmt:
                    self.cur.execute(stmt)

            self.conn.commit()
            print(f"{sql_file}: successfully executed")
            return True
        except Exception as e:
            self.conn.rollback()
            print(f"{sql_file}: failed to execute: {e}")
            return False

    def import_all_tables(self):
        """import all tables"""
        print(f"importing from directory: {self.export_dir}")

        # read metadata
        metadata_file = os.path.join(self.export_dir, "metadata.json")
        if os.path.exists(metadata_file):
            with open(metadata_file, "r") as f:
                metadata = json.load(f)
            print(f"database: {metadata.get('database')}")
            print(f"export time: {metadata.get('export_time')}")

        # 1. create table structure
        print("\n1. creating table structure...")
        self.execute_sql_file(os.path.join(self.export_dir, "schema_ddl.sql"))

        # 2. import data
        print("\n2. importing data...")

        # get table list
        self.cur.execute(
            """
            SELECT table_name 
            FROM information_schema.tables 
            WHERE table_schema = 'public' 
            AND table_type = 'BASE TABLE'
            ORDER BY table_name;
        """
        )

        tables = [row[0] for row in self.cur.fetchall()]

        success_count = 0
        for table in tables:
            if self.import_data_from_csv(table):
                success_count += 1

        # 3. 创建索引（数据导入后创建，提高导入性能）
        print("\n3. creating indexes...")
        self.execute_sql_file(os.path.join(self.export_dir, "indexes.sql"))

        # 4. add foreign key constraints
        print("\n4. adding foreign key constraints...")
        self.execute_sql_file(os.path.join(self.export_dir, "foreign_keys.sql"))

        # 5. create sequences
        print("\n5. creating sequences...")
        self.execute_sql_file(os.path.join(self.export_dir, "sequences.sql"))

        print(f"\n{'=' * 50}")
        print(f"import completed successfully!")
        print(f"successfully imported: {success_count}/{len(tables)} tables of data")
        print(f"{'=' * 50}")

    def verify_import(self):
        """verify imported data"""
        print("\nverifying imported data:")

        schema = self.read_schema()
        for table_name in schema.keys():
            self.cur.execute(f"SELECT COUNT(*) FROM {table_name}")
            count = self.cur.fetchone()[0]

            csv_file = os.path.join(self.export_dir, f"{table_name}.csv")
            if os.path.exists(csv_file):
                with open(csv_file, "r", encoding="utf-8") as f:
                    csv_count = sum(1 for _ in f) - 1  # subtract header row
                print(
                    f"table {table_name}: imported {count} rows, csv has {csv_count} rows"
                )

    def close(self):
        """close database connection"""
        if hasattr(self, "cur"):
            self.cur.close()
        if hasattr(self, "conn"):
            self.conn.close()


# example
if __name__ == "__main__":
    # parameters for database connection
    DB_CONFIG = {
        "host": "192.168.100.2",
        "port": "5432",
        "user": "postgres",
        "password": "T1@2n3d4emai",
        "database": "gen-chembl_tmp",
    }

    # specify export directory, if not specified, use the latest one
    EXPORT_DIR = "/nfs/workspace/yinhuiyi/codedir/CrEM-Database/tmp/export_20251205_143245"  # modify to your export directory

    # create importer and execute import
    importer = PostgresImporter(export_dir=EXPORT_DIR, **DB_CONFIG)

    try:
        importer.connect()
        importer.import_all_tables()
        importer.verify_import()
    except Exception as e:
        print(f"importing process failed: {e}")
    finally:
        importer.close()
