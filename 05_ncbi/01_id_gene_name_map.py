#! /usr/bin/python3
"""Fill genes.ncbi_gene_id from the NCBI gene_info file, and report where NCBI and HGNC disagree.

usage:  01_id_gene_name_map.py [Homo_sapiens.gene_info.gz]

Selection: human records whose HGNC cross-reference names a gene currently in 'genes'. NCBI also
carries entries HGNC has withdrawn since this gene_info was cut; those are outside the selection,
are counted in the report, and create nothing - 01_hgnc remains the only writer of gene rows.

Matching is on the HGNC id and not on the symbol: NCBI and HGNC do not always agree on the
symbol, and that disagreement is what this script reports rather than something it works around.

The download:
  wget -O /storage/databases/ncbi/Homo_sapiens.gene_info.gz \\
    https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
"""

import gzip
import os
import re
import sys
from typing import Dict, Iterator, List, Optional, Tuple

from sqlmodel import Session, select

from integrator_utils.python.db import db_engine, require_genes
from integrator_utils.python.models_core import Gene

GENE_INFO = "/storage/databases/ncbi/Homo_sapiens.gene_info.gz"
HGNC_XREF_PATTERN = re.compile(r"HGNC:(HGNC:\d+)")
CHUNK = 5000


#########################################
def gene_info_records(path: str) -> Iterator[Tuple[str, int, str]]:
    """(hgnc id, ncbi gene id, ncbi symbol) for every human record cross-referenced to HGNC."""
    with gzip.open(path, "rt") as inf:
        for line in inf:
            if line.startswith("#"):
                continue
            field = line.rstrip("\n").split("\t")
            if len(field) < 11 or field[0] != "9606":
                continue
            xref_match = HGNC_XREF_PATTERN.search(field[5])
            if not xref_match or not field[1].isdigit():
                continue
            symbol = field[10] if field[10] != "-" else field[2]
            yield xref_match.group(1), int(field[1]), symbol


#########################################
def update(session: Session, records: List[Tuple[str, int, str]]) -> Tuple[int, int, List[str]]:
    hgnc_ids = [hgnc_id for hgnc_id, _, _ in records]
    genes = {gene.hgnc_id: gene for gene in session.exec(select(Gene).where(Gene.hgnc_id.in_(hgnc_ids)))}
    stored = 0
    unknown = 0
    disagreement: List[str] = []
    for hgnc_id, ncbi_gene_id, symbol in records:
        gene = genes.get(hgnc_id)
        if gene is None:  # an HGNC entry withdrawn since this gene_info was cut - not our row to make
            unknown += 1
            continue
        if gene.symbol.lower() != symbol.lower():
            note = f"{hgnc_id}: HGNC '{gene.symbol}' vs NCBI '{symbol}'"
            disagreement.append(note)
        if gene.ncbi_gene_id != ncbi_gene_id:
            gene.ncbi_gene_id = ncbi_gene_id
            session.add(gene)
            stored += 1
    session.commit()
    return stored, unknown, disagreement


def load(engine, path: str) -> None:
    seen = 0
    stored = 0
    unknown = 0
    disagreements: List[str] = []
    chunk: List[Tuple[str, int, str]] = []
    with Session(engine) as session:
        for record in gene_info_records(path):
            chunk.append(record)
            if len(chunk) < CHUNK:
                continue
            updated, missing, notes = update(session, chunk)
            stored += updated
            unknown += missing
            disagreements += notes
            seen += len(chunk)
            chunk = []
        if chunk:
            updated, missing, notes = update(session, chunk)
            stored += updated
            unknown += missing
            disagreements += notes
            seen += len(chunk)

    report = f"{seen} hgnc cross-referenced records read, {stored} ncbi gene ids stored"
    print(report)
    skipped = f"{unknown} record(s) refer to an hgnc id that is not in 'genes' (withdrawn since)"
    print(skipped)
    mismatch = f"{len(disagreements)} symbol disagreement(s) between HGNC and NCBI"
    print(mismatch)
    for note in disagreements[:20]:
        print(f"\t {note}")


#########################################
def main():
    if len(sys.argv) > 2:
        usage = f"usage: {sys.argv[0]} [Homo_sapiens.gene_info.gz]"
        sys.exit(usage)
    path = sys.argv[1] if len(sys.argv) == 2 else GENE_INFO
    if not os.path.exists(path):
        errmsg = f"{path} not found - see the download command in the docstring of this script"
        sys.exit(errmsg)

    engine = db_engine()
    require_genes(engine)
    load(engine, path)


#########################################
if __name__ == '__main__':
    main()
