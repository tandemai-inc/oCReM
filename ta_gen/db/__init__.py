#!/usr/bin/env python
# -*- coding:utf-8 -*-


def create_db_manager(db_type, db_path=None, ini_file=None):
    if db_type == "sqlite":
        from ta_gen.db.sqlite3_manager import SqliteManager

        return SqliteManager(db_path)
    elif db_type == "postgres":
        from ta_gen.db.postgres_manager import PostGresManager

        return PostGresManager(ini_file)
