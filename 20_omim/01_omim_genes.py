#! /usr/bin/python3
"""Fill 'omim_genes' from the OMIM mim2gene file.

usage:  01_omim_genes.py [mim2gene.txt]

Selection: every mim2gene record; those whose Entrez gene id is in 'genes' are attached to it.

Attachment is by Entrez id, not by the 'Approved Gene Symbol (HGNC)' column, and that is not a
detail: this file is from 2019, and 353 of its symbols no longer resolve. Some are genes renamed
since (QARS -> QARS1, SEPT2 -> SEPTIN2, TAZ -> TAFAZZIN), others are phenotype locus names HGNC
never approved as genes at all (DFNB5, MYP2, SCZD1). Every symbol-bearing record carries an
Entrez id, and that id survived all of those renamings.

Records whose Entrez id is not in 'genes' - phenotype loci, and genes withdrawn since - keep
gene_id null and are counted in the report. Nothing here creates a gene row.

The source: https://omim.org/downloads (mim2gene.txt is the one file OMIM serves without a
licence key).
"""

import os
import sys
from typing import Any, Dict, Iterator, List, Optional

from sqlmodel import Session, select

from integrator_utils.python.db import db_engine, require_genes, upsert_many
from integrator_utils.python.models_core import Gene, OmimGene

MIM2GENE = "/storage/databases/omim/mim2gene.txt"


#########################################
def records(path: str) -> Iterator[Dict[str, Any]]:
    with open(path) as inf:
        for line in inf:
            if line.startswith("#"):
                continue
            field = (line.rstrip("\n").split("\t") + [""] * 5)[:5]
            if not field[0].isdigit():
                continue
            yield {"mim_number": int(field[0]), "entry_type": field[1].strip(),
                   "symbol": field[3].strip() or None,
                   "ncbi_gene_id": int(field[2]) if field[2].strip().isdigit() else None,
                   "ensembl_gene_id": field[4].strip() or None}


#########################################
def load(engine, path: str) -> None:
    rows = list(records(path))
    with Session(engine) as session:
        gene_id_of = {gene.ncbi_gene_id: gene.id for gene in session.exec(select(Gene)) if gene.ncbi_gene_id}
        attached = 0
        unattached: List[str] = []
        for row in rows:
            symbol = row.pop("symbol")
            row["gene_id"] = gene_id_of.get(row["ncbi_gene_id"])
            if row["gene_id"] is not None:
                attached += 1
            elif symbol:
                unattached.append(f"{row['mim_number']} ({symbol}, {row['entry_type']})")
        upsert_many(session, OmimGene, rows, conflict_on=["mim_number"])
        session.commit()

    report = f"{len(rows)} omim entries stored, {attached} of them attached to a gene"
    print(report)
    orphans = f"{len(unattached)} record(s) name a symbol whose entrez id is not in 'genes'"
    print(orphans)
    for note in unattached[:20]:
        print(f"\t {note}")


#########################################
def main():
    if len(sys.argv) > 2:
        usage = f"usage: {sys.argv[0]} [mim2gene.txt]"
        sys.exit(usage)
    path = sys.argv[1] if len(sys.argv) == 2 else MIM2GENE
    if not os.path.exists(path):
        errmsg = f"{path} not found - see the source url in the docstring of this script"
        sys.exit(errmsg)

    engine = db_engine()
    require_genes(engine)
    load(engine, path)


#########################################
if __name__ == '__main__':
    main()
