#!/usr/bin/env python
# -*- coding:utf-8 -*-


def create_db_manager(db_type, db_path=None, ini_file=None, reset_db=False):
    if db_type == "sqlite":
        from ta_gen.db.sqlite3_manager import SqliteManager

        return SqliteManager(db_path, reset_db)
    elif db_type == "postgres":
        from ta_gen.db.postgres_manager import PostGresManager

        return PostGresManager(ini_file, reset_db)
