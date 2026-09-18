"""Alembic environment: the url comes from .env, the metadata from every models module in the repo.

The numbered source directories are not importable packages, so their '<source>_models.py' files
are loaded by path. A new directory needs no edit here - it only has to follow the naming.
"""

import glob
import importlib.util
import os
import sys

from alembic import context
from logging.config import fileConfig
from sqlalchemy import engine_from_config, pool
from sqlmodel import SQLModel

from integrator_utils.python.db import REPO_ROOT, db_url
from integrator_utils.python import models_core  # noqa: F401  - registers the core tables

config = context.config
if config.config_file_name is not None:
    fileConfig(config.config_file_name)

config.set_main_option("sqlalchemy.url", db_url().replace("%", "%%"))


#########################################
def import_directory_models() -> None:
    """Import every <nn>_<source>/<source>_models.py, so autogenerate sees all tables at once."""
    for path in sorted(glob.glob(os.path.join(REPO_ROOT, "[0-9][0-9]_*", "*_models.py"))):
        module_name = os.path.splitext(os.path.basename(path))[0]
        if module_name in sys.modules:
            continue
        spec = importlib.util.spec_from_file_location(module_name, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = module
        spec.loader.exec_module(module)


import_directory_models()
target_metadata = SQLModel.metadata


#########################################
def run_migrations_offline() -> None:
    url = config.get_main_option("sqlalchemy.url")
    context.configure(url=url, target_metadata=target_metadata, literal_binds=True,
                      dialect_opts={"paramstyle": "named"}, compare_type=True)
    with context.begin_transaction():
        context.run_migrations()


def run_migrations_online() -> None:
    section = config.get_section(config.config_ini_section, {})
    connectable = engine_from_config(section, prefix="sqlalchemy.", poolclass=pool.NullPool)
    with connectable.connect() as connection:
        context.configure(connection=connection, target_metadata=target_metadata, compare_type=True)
        with context.begin_transaction():
            context.run_migrations()


#########################################
if context.is_offline_mode():
    run_migrations_offline()
else:
    run_migrations_online()
