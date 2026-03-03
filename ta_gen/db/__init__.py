#!/usr/bin/env python
# -*- coding:utf-8 -*-




def create_db_manager(db_type):
    if db_type == "sqlite":
        from ta_gen.db.sqlite3_manager import SqliteManager

        return SqliteManager()
    elif db_type == "postgres":
        from ta_gen.db.postgres_manager import PostGresManager

        return PostGresManager()
