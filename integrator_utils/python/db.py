#! /usr/bin/python3
"""Engine, session and bulk helpers shared by all loaders, plus the gene resolution rules.

Two rules are enforced here, and they are the reason this module exists:

  * the schema belongs to alembic - no loader calls SQLModel.metadata.create_all();
  * 01_hgnc is the only writer of 'genes' - every other loader resolves symbols through
    gene_id()/gene_ids(), and exits if the table or the symbol is missing.
"""

import os
import sys
from contextlib import contextmanager
from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence, Set, Type

from dotenv import load_dotenv
from sqlalchemy import Engine, inspect
from sqlalchemy.dialects.postgresql import insert
from sqlmodel import Session, SQLModel, create_engine, select

from integrator_utils.python.models_core import Gene, GeneAlias, HgncWithdrawn

# realpath, not abspath: several directories still hold an 'integrator_utils' symlink from the
# days before the package was installed, and the repo root has to come out the same either way
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))


#########################################
# engine and session
def db_url() -> str:
    load_dotenv(os.path.join(REPO_ROOT, ".env"))
    password = os.getenv("POSTGRES_PASSWD")
    if not password:
        errmsg = "POSTGRES_PASSWD not found in the environment (.env)"
        sys.exit(errmsg)
    user = os.getenv("POSTGRES_USER", "postgres").strip('"')
    host = os.getenv("POSTGRES_HOST", "127.0.0.1").strip('"')
    port = os.getenv("POSTGRES_PORT", "5432").strip('"')
    name = os.getenv("POSTGRES_DB", "integrator").strip('"')
    return f"postgresql+psycopg://{user}:{password}@{host}:{port}/{name}"


def db_engine(echo: bool = False) -> Engine:
    return create_engine(db_url(), echo=echo)


@contextmanager
def db_session(engine: Optional[Engine] = None) -> Iterator[Session]:
    """Commit on a clean exit, roll back on any exception - a half-loaded table helps nobody."""
    with Session(engine if engine is not None else db_engine()) as session:
        try:
            yield session
            session.commit()
        except Exception:
            session.rollback()
            raise


#########################################
# write helpers - every loader is idempotent, so every write is an upsert
def upsert(session: Session, model: Type[SQLModel], values: Dict[str, Any], conflict_on: Sequence[str]) -> int:
    """Insert or update a single row; return its id."""
    updatable = {field: value for field, value in values.items() if field not in conflict_on}
    statement = insert(model).values(values)
    if updatable:
        statement = statement.on_conflict_do_update(index_elements=list(conflict_on), set_=updatable)
    else:
        statement = statement.on_conflict_do_nothing(index_elements=list(conflict_on))
    return session.exec(statement.returning(model.__table__.c.id)).one()[0]


def batches(rows: List[Dict[str, Any]]) -> Iterator[List[Dict[str, Any]]]:
    """Split so that no statement exceeds the 65535 bound parameters postgres accepts."""
    per_statement = max(1, 65000 // max(1, len(rows[0])))
    for start in range(0, len(rows), per_statement):
        yield rows[start:start + per_statement]


def upsert_many(session: Session, model: Type[SQLModel], rows: List[Dict[str, Any]],
                conflict_on: Sequence[str]) -> None:
    """Insert or update. The caller must not pass the same key twice within one batch."""
    if not rows:
        return
    for batch in batches(rows):
        statement = insert(model).values(batch)
        updatable = {column: statement.excluded[column] for column in batch[0] if column not in conflict_on}
        if updatable:
            statement = statement.on_conflict_do_update(index_elements=list(conflict_on), set_=updatable)
        else:
            statement = statement.on_conflict_do_nothing(index_elements=list(conflict_on))
        session.exec(statement)


def insert_ignore(session: Session, model: Type[SQLModel], rows: List[Dict[str, Any]],
                  conflict_on: Sequence[str]) -> None:
    if not rows:
        return
    for batch in batches(rows):
        statement = insert(model).values(batch).on_conflict_do_nothing(index_elements=list(conflict_on))
        session.exec(statement)


def copy_from_tsv(engine: Engine, model: Type[SQLModel], path: str, columns: Sequence[str]) -> int:
    """COPY a headerless tsv into a table - the only sane way to load gnomad-scale data."""
    column_list = ", ".join(columns)
    statement = f"COPY {model.__tablename__} ({column_list}) FROM STDIN WITH (FORMAT text, NULL '\\N')"
    raw = engine.raw_connection()
    try:
        with raw.cursor() as cursor:
            with cursor.copy(statement) as copy, open(path, "rb") as inf:
                for block in iter(lambda: inf.read(1 << 20), b""):
                    copy.write(block)
            rowcount = cursor.rowcount
        raw.commit()
    finally:
        raw.close()
    return rowcount


#########################################
# the gene hub - resolution only, never creation
def require_genes(engine: Engine) -> None:
    """Called once by every loader except 01_hgnc, before any parsing is done."""
    missing = [table for table in ("genes", "gene_aliases") if not inspect(engine).has_table(table)]
    if missing:
        tables = " and ".join(f"'{table}'" for table in missing)
        errmsg = f"the {tables} table does not exist - run 'alembic upgrade head', then 01_hgnc"
        sys.exit(errmsg)
    with Session(engine) as session:
        if session.exec(select(Gene).limit(1)).first() is None:
            errmsg = "the 'genes' table is empty - run 01_hgnc/01_hgnc_genes.py first"
            sys.exit(errmsg)


def resolve(session: Session, symbol: str) -> Optional[int]:
    """Approved symbol first, then aliases and previous symbols; None if unresolvable."""
    gene = session.exec(select(Gene).where(Gene.symbol == symbol)).first()
    if gene is not None:
        return gene.id
    gene_ids_found = {alias.gene_id for alias in session.exec(select(GeneAlias).where(GeneAlias.alias == symbol))}
    if len(gene_ids_found) == 1:
        return gene_ids_found.pop()
    return None  # unknown, or an alias shared by several genes - the caller decides how loudly to die


def withdrawn_hgnc_ids(session: Session) -> Set[str]:
    """The HGNC ids HGNC itself has retired - see HgncWithdrawn for why a loader wants these."""
    return {row.hgnc_id for row in session.exec(select(HgncWithdrawn))}


def gene_id(session: Session, symbol: str) -> int:
    resolved = resolve(session, symbol)
    if resolved is None:
        errmsg = f"gene symbol '{symbol}' does not resolve to a single row of 'genes' - not an HGNC symbol?"
        sys.exit(errmsg)
    return resolved


def gene_id_map(session: Session) -> Dict[str, int]:
    """Every resolvable name - approved symbol, then unambiguous alias - lowercased, in one dict.

    Worth the memory (45k genes, 59k aliases) as soon as a loader has more than a few hundred
    symbols to resolve: the alternative is a query per symbol.
    """
    by_name: Dict[str, int] = {}
    ambiguous: Set[str] = set()
    for alias in session.exec(select(GeneAlias)):
        name = alias.alias.lower()
        if name in by_name and by_name[name] != alias.gene_id:
            ambiguous.add(name)
        by_name[name] = alias.gene_id
    for name in ambiguous:
        del by_name[name]
    for gene in session.exec(select(Gene)):  # an approved symbol always wins over an alias
        by_name[gene.symbol.lower()] = gene.id
    return by_name


def gene_ids(session: Session, symbols: Iterable[str]) -> Dict[str, int]:
    """Resolve a whole symbol set up front: a bulk load should fail in seconds, not halfway."""
    by_name = gene_id_map(session)
    resolved: Dict[str, int] = {}
    unresolved: List[str] = []
    for symbol in dict.fromkeys(symbols):  # de-duplicated, order kept
        found = by_name.get(symbol.lower())
        if found is None:
            unresolved.append(symbol)
        else:
            resolved[symbol] = found
    if unresolved:
        offenders = ", ".join(unresolved)
        errmsg = f"{len(unresolved)} gene symbol(s) do not resolve to a single row of 'genes': {offenders}"
        sys.exit(errmsg)
    return resolved
