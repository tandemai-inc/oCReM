#!/usr/bin/env python
# -*- coding:utf-8 -*-


from ta_gen.db.postgres_manager import PostGresManager
from ta_gen.db.sqlite3_manager import SqliteManager


def create_db_manager(db_type):
    if db_type == "sqlite":
        return SqliteManager()
    elif db_type == "postgres":
        return PostGresManager()
