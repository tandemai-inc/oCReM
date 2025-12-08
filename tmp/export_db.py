#!/usr/bin/env python
# -*- coding:utf-8 -*-

import csv
import json
import os
from datetime import datetime

import psycopg2
import psycopg2.extras


class PostgresExporter:
    def __init__(self, host, database, user, password, port=5432):
        self.connection_params = {
            "host": host,
            "database": database,
            "user": user,
            "password": password,
            "port": port,
        }
        self.export_dir = f"export_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

    def connect(self):
        """create db connection"""
        try:
            self.conn = psycopg2.connect(**self.connection_params)
            self.cur = self.conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
            print("successfully connected to db")
        except Exception as e:
            print(f"failed to connect db: {e}")
            raise

    def export_schema_ddl(self):
        """export complete ddl (table schema, indexes, constraints)"""
        ddl_file = os.path.join(self.export_dir, "schema_ddl.sql")

        # 1. export all table ddl (including primary keys, foreign keys, constraints)
        self.cur.execute(
            """
            SELECT 
                'CREATE TABLE ' || schemaname || '.' || tablename || ' (' ||
                string_agg(column_definition, ', ') || 
                COALESCE(', ' || constraint_definition, '') ||
                ');' as create_table_ddl
            FROM (
                SELECT 
                    c.table_schema as schemaname,
                    c.table_name as tablename,
                    c.column_name || ' ' || c.data_type || 
                    COALESCE('(' || c.character_maximum_length || ')', '') || 
                    CASE WHEN c.is_nullable = 'NO' THEN ' NOT NULL' ELSE '' END ||
                    COALESCE(' DEFAULT ' || c.column_default, '') as column_definition,
                    con.constraint_definition
                FROM information_schema.columns c
                LEFT JOIN (
                    SELECT 
                        tc.table_schema, tc.table_name,
                        'PRIMARY KEY (' || string_agg(kcu.column_name, ', ') || ')' as constraint_definition
                    FROM information_schema.table_constraints tc
                    JOIN information_schema.key_column_usage kcu 
                        ON tc.constraint_name = kcu.constraint_name
                    WHERE tc.constraint_type = 'PRIMARY KEY'
                    GROUP BY tc.table_schema, tc.table_name, tc.constraint_name
                ) con ON c.table_schema = con.table_schema AND c.table_name = con.table_name
                WHERE c.table_schema = 'public'
                ORDER BY c.table_schema, c.table_name, c.ordinal_position
            ) t
            GROUP BY schemaname, tablename, constraint_definition;
        """
        )

        with open(ddl_file, "w", encoding="utf-8") as f:
            # write table ddl
            for row in self.cur.fetchall():
                f.write(row["create_table_ddl"] + "\n\n")

        print(f"table schema ddl exported to: {ddl_file}")

    def export_table_data(self, table_name):
        """export table data to csv"""
        csv_file = os.path.join(self.export_dir, "data", f"{table_name}.csv")

        # ensure data directory exists
        os.makedirs(os.path.dirname(csv_file), exist_ok=True)

        try:
            # use COPY command to export (most efficient)
            with open(csv_file, "w", encoding="utf-8") as f:
                self.cur.copy_expert(f"COPY {table_name} TO STDOUT WITH CSV HEADER", f)
            print(f"{table_name}: data exported successfully")
            return True
        except Exception as e:
            # use SELECT if COPY fails
            print(f"{table_name}: COPY failed, using SELECT: {e}")
            try:
                self.cur.execute(f"SELECT * FROM {table_name}")
                rows = self.cur.fetchall()

                if not rows:
                    print(f"{table_name}: is empty")
                    return True

                column_names = [desc[0] for desc in self.cur.description]

                with open(csv_file, "w", newline="", encoding="utf-8") as f:
                    writer = csv.writer(f)
                    writer.writerow(column_names)
                    for row in rows:
                        writer.writerow(row)

                print(f"{table_name}: data exported successfully ({len(rows)} rows)")
                return True
            except Exception as e2:
                print(f"{table_name}: data export failed: {e2}")
                return False

    def export_foreign_keys(self):
        """export foreign key constraints"""
        fk_file = os.path.join(self.export_dir, "foreign_keys.sql")

        self.cur.execute(
            """
            SELECT
                tc.table_schema,
                tc.table_name,
                tc.constraint_name,
                kcu.column_name,
                ccu.table_schema AS foreign_table_schema,
                ccu.table_name AS foreign_table_name,
                ccu.column_name AS foreign_column_name
            FROM information_schema.table_constraints AS tc
            JOIN information_schema.key_column_usage AS kcu
                ON tc.constraint_name = kcu.constraint_name
            JOIN information_schema.constraint_column_usage AS ccu
                ON ccu.constraint_name = tc.constraint_name
            WHERE tc.constraint_type = 'FOREIGN KEY' AND tc.table_schema = 'public';
        """
        )

        fks = self.cur.fetchall()

        with open(fk_file, "w", encoding="utf-8") as f:
            for fk in fks:
                sql = f"ALTER TABLE {fk['table_schema']}.{fk['table_name']} "
                sql += f"ADD CONSTRAINT {fk['constraint_name']} "
                sql += f"FOREIGN KEY ({fk['column_name']}) "
                sql += f"REFERENCES {fk['foreign_table_schema']}.{fk['foreign_table_name']} "
                sql += f"({fk['foreign_column_name']});\n"
                f.write(sql)

        print(f"foreign key constraints exported to: {fk_file}")

    def export_sequences(self):
        """export sequences"""
        seq_file = os.path.join(self.export_dir, "sequences.sql")

        self.cur.execute(
            """
            SELECT 
                schemaname, sequencename, 
                start_value, minimum_value, maximum_value,
                increment, cycle_option
            FROM information_schema.sequences
            WHERE schemaname = 'public';
        """
        )

        sequences = self.cur.fetchall()

        with open(seq_file, "w", encoding="utf-8") as f:
            for seq in sequences:
                sql = f"CREATE SEQUENCE {seq['schemaname']}.{seq['sequencename']} "
                sql += f"START WITH {seq['start_value']} "
                sql += f"INCREMENT BY {seq['increment']} "
                sql += f"MINVALUE {seq['minimum_value']} "
                sql += f"MAXVALUE {seq['maximum_value']} "
                sql += f"{'CYCLE' if seq['cycle_option'] == 'YES' else 'NO CYCLE'};\n"
                f.write(sql)

        print(f"sequences exported to: {seq_file}")

    def export_metadata(self, tables, success_count):
        """export metadata"""
        metadata = {
            "export_time": datetime.now().isoformat(),
            "database": self.connection_params["database"],
            "host": self.connection_params["host"],
            "total_tables": len(tables),
            "exported_tables": success_count,
            "tables_list": tables,
            "postgres_version": self.get_postgres_version(),
        }

        with open(os.path.join(self.export_dir, "metadata.json"), "w") as f:
            json.dump(metadata, f, indent=2, default=str)
        print(f"metadata exported to: {os.path.join(self.export_dir, 'metadata.json')}")

    def export_all_tables(self):
        os.makedirs(os.path.join(self.export_dir, "data"), exist_ok=True)
        print(f"starting exporting all tables to: {self.export_dir}")

        # export table schema
        self.export_schema_ddl()

        # 2. export all table data
        self.cur.execute(
            """
            SELECT table_name 
            FROM information_schema.tables 
            WHERE table_schema = 'public' 
            AND table_type = 'BASE TABLE'
            ORDER BY table_name;
        """
        )

        tables = [row["table_name"] for row in self.cur.fetchall()]
        print(f"found {len(tables)} tables")

        # 3. export all table data
        print("exporting table data...")
        success_count = 0
        for table in tables:
            if self.export_table_data(table):
                success_count += 1

        # 4. export indexes
        print("\nexporting indexes...")
        self.export_indexes()

        # 5. export foreign key constraints
        print("exporting foreign key constraints...")
        self.export_foreign_keys()

        # 6. export sequences
        print("exporting sequences...")
        self.export_sequences()

        # 7. export metadata
        print("exporting metadata...")
        self.export_metadata(tables, success_count)

        print(f"\n{'=' * 50}")
        print(f"export completed successfully!")
        print(f"export directory: {self.export_dir}")
        print(f"successfully exported: {success_count}/{len(tables)} tables")
        print(f"{'=' * 50}")

    def close(self):
        """close db connection"""
        if hasattr(self, "cur"):
            self.cur.close()
        if hasattr(self, "conn"):
            self.conn.close()


# export db example
if __name__ == "__main__":
    DB_CONFIG = {
        "host": "192.168.100.2",
        "port": "5432",
        "user": "postgres",
        "password": "T1@2n3d4emai",
        "database": "gen-chembl",
    }

    exporter = PostgresExporter(**DB_CONFIG)

    try:
        exporter.connect()
        exporter.export_all_tables()
    except Exception as e:
        print(f"export process failed: {e}")
    finally:
        exporter.close()
