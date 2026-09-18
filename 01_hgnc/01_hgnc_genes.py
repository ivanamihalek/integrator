#! /usr/bin/python3
"""Fill 'genes' and 'gene_aliases' from the HGNC complete set - the first loader to run.

usage:  01_hgnc_genes.py [hgnc_complete_set.txt [withdrawn.txt]]

HGNC is the authority for gene identity in this database: this is the only script that creates
gene rows. Every other loader resolves symbols against what this one stores, and exits if it
cannot. Re-running replaces the alias set of each gene it sees, so a fresher HGNC download
propagates without leaving stale aliases behind.

The withdrawn list goes into 'hgnc_withdrawn' in the same run. Later loaders need it to tell a
source record that names a retired gene - which is outside their selection and fine - from one
that names an id nobody has ever heard of, which is not fine and stops them.

The downloads, both from https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/ :
  hgnc_complete_set.txt
  withdrawn.txt
"""

import csv
import os
import sys
from typing import Any, Dict, Iterator, List, Optional

from sqlmodel import Session, delete, select

from integrator_utils.python.db import db_engine, upsert_many
from integrator_utils.python.models_core import Gene, GeneAlias, HgncWithdrawn

HGNC_COMPLETE_SET = "/storage/databases/hgnc/hgnc_complete_set.txt"
HGNC_WITHDRAWN = "/storage/databases/hgnc/withdrawn.txt"
CHUNK = 2000

# hgnc column -> model field, for the columns that go into 'genes' as they stand
SCALAR_FIELDS: Dict[str, str] = {
    "hgnc_id": "hgnc_id",
    "symbol": "symbol",
    "name": "name",
    "locus_group": "locus_group",
    "locus_type": "locus_type",
    "location": "location",
    "ensembl_gene_id": "ensembl_gene_id",
    "ucsc_id": "ucsc_id",
}

# hgnc column -> model field, for the pipe-joined columns that become text arrays
ARRAY_FIELDS: Dict[str, str] = {
    "uniprot_ids": "uniprot_ids",
    "refseq_accession": "refseq_ids",
    "ccds_id": "ccds_ids",
    "mgd_id": "mgd_ids",
}

# hgnc column -> gene_aliases.alias_type
ALIAS_FIELDS: Dict[str, str] = {
    "alias_symbol": "alias",
    "prev_symbol": "previous",
}


#########################################
def split_list(value: str) -> Optional[List[str]]:
    items = [item.strip() for item in value.split("|") if item.strip()]
    return items or None


def gene_row(record: Dict[str, str]) -> Dict[str, Any]:
    row: Dict[str, Any] = {field: (record.get(column) or None) for column, field in SCALAR_FIELDS.items()}
    row.update({field: split_list(record.get(column, "")) for column, field in ARRAY_FIELDS.items()})
    entrez_id = record.get("entrez_id", "").strip()
    row["ncbi_gene_id"] = int(entrez_id) if entrez_id.isdigit() else None
    return row


def alias_rows(record: Dict[str, str], gene_id: int) -> List[Dict[str, Any]]:
    """One row per (alias, type), de-duplicated: HGNC occasionally lists a name twice."""
    seen: Dict[str, Dict[str, Any]] = {}
    for column, alias_type in ALIAS_FIELDS.items():
        for alias in split_list(record.get(column, "")) or []:
            key = f"{alias.lower()}|{alias_type}"
            if key not in seen:
                seen[key] = {"gene_id": gene_id, "alias": alias, "alias_type": alias_type}
    return list(seen.values())


def approved_records(path: str) -> Iterator[Dict[str, str]]:
    """Approved symbols only - withdrawn entries have no place in a table used for resolution."""
    with open(path, newline="") as inf:
        for record in csv.DictReader(inf, delimiter="\t"):
            if record.get("status") != "Approved" or not record.get("symbol"):
                continue
            yield record


#########################################
def store(session: Session, records: List[Dict[str, str]]) -> int:
    upsert_many(session, Gene, [gene_row(record) for record in records], conflict_on=["hgnc_id"])
    session.commit()

    hgnc_ids = [record["hgnc_id"] for record in records]
    gene_id_of = {gene.hgnc_id: gene.id for gene in session.exec(select(Gene).where(Gene.hgnc_id.in_(hgnc_ids)))}
    session.exec(delete(GeneAlias).where(GeneAlias.gene_id.in_(list(gene_id_of.values()))))

    aliases: List[Dict[str, Any]] = []
    for record in records:
        aliases.extend(alias_rows(record, gene_id_of[record["hgnc_id"]]))
    upsert_many(session, GeneAlias, aliases, conflict_on=["gene_id", "alias", "alias_type"])
    session.commit()
    return len(aliases)


def load(engine, path: str) -> None:
    genes = 0
    aliases = 0
    chunk: List[Dict[str, str]] = []
    with Session(engine) as session:
        for record in approved_records(path):
            chunk.append(record)
            if len(chunk) < CHUNK:
                continue
            aliases += store(session, chunk)
            genes += len(chunk)
            chunk = []
            progress = f"\t {genes} genes"
            print(progress)
        if chunk:
            aliases += store(session, chunk)
            genes += len(chunk)
    report = f"{genes} genes, {aliases} aliases and previous symbols stored"
    print(report)


def load_withdrawn(engine, path: str) -> None:
    rows: List[Dict[str, Any]] = []
    with open(path, newline="") as inf:
        for record in csv.DictReader(inf, delimiter="\t"):
            hgnc_id = (record.get("HGNC_ID") or "").strip()
            if not hgnc_id:
                continue
            merged = record.get("MERGED_INTO_REPORT(S) (i.e HGNC_ID|SYMBOL|STATUS)") or ""
            # "HGNC:123|SYMBOL|Approved, HGNC:456|..." - keep the ids the entry was merged into
            merged_into = [item.split("|")[0].strip() for item in merged.split(",") if item.strip()]
            rows.append({"hgnc_id": hgnc_id, "status": (record.get("STATUS") or "").strip(),
                         "withdrawn_symbol": (record.get("WITHDRAWN_SYMBOL") or "").strip() or None,
                         "merged_into": merged_into or None})
    with Session(engine) as session:
        upsert_many(session, HgncWithdrawn, rows, conflict_on=["hgnc_id"])
        session.commit()
    report = f"{len(rows)} withdrawn hgnc entries stored"
    print(report)


#########################################
def main():
    if len(sys.argv) > 3:
        usage = f"usage: {sys.argv[0]} [hgnc_complete_set.txt [withdrawn.txt]]"
        sys.exit(usage)
    path = sys.argv[1] if len(sys.argv) > 1 else HGNC_COMPLETE_SET
    if not os.path.exists(path):
        errmsg = f"{path} not found - see the download command in the docstring of this script"
        sys.exit(errmsg)

    engine = db_engine()
    load(engine, path)
    withdrawn = sys.argv[2] if len(sys.argv) == 3 else HGNC_WITHDRAWN
    if not os.path.exists(withdrawn):
        errmsg = f"{withdrawn} not found - see the download commands in the docstring of this script"
        sys.exit(errmsg)
    load_withdrawn(engine, withdrawn)


#########################################
if __name__ == '__main__':
    main()
