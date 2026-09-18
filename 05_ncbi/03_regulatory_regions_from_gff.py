#! /usr/bin/python3
"""Store the RefSeq regulatory features of GRCh37 in 'ncbi_regulatory_regions'.

usage:  03_regulatory_regions_from_gff.py [--bed <output.bed>]

Selection: features carrying a regulatory_class attribute, on assembled molecules only - the
alt loci and patch scaffolds of the gff duplicate them in coordinates nothing else here uses.

No gene is attached: the Dbxref of a regulatory feature is that feature's own GeneID, not the
gene it acts on.

The downloads, both from
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/
  GCF_000001405.25_GRCh37.p13_genomic.gff.gz
  GCF_000001405.25_GRCh37.p13_assembly_report.txt
"""

import os
import sys
from typing import Any, Dict, Iterator, List, Optional, TextIO

from sqlalchemy.dialects.postgresql import Range
from sqlmodel import Session

from integrator_utils.python.db import db_engine, upsert_many

import gff
from ncbi_models import RegulatoryRegion

HG19 = "/storage/databases/ncbi/hg19"
GFF = os.path.join(HG19, "GCF_000001405.25_GRCh37.p13_genomic.gff.gz")
ASSEMBLY_REPORT = os.path.join(HG19, "GCF_000001405.25_GRCh37.p13_assembly_report.txt")
CHUNK = 5000


#########################################
def regulatory_rows(path: str, chrom_of: Dict[str, str], bed: Optional[TextIO]) -> Iterator[Dict[str, Any]]:
    for field in gff.lines(path):
        chrom = chrom_of.get(field[0])
        if chrom is None or "regulatory_class=" not in field[8]:
            continue
        attribute = gff.attributes(field[8])
        # gff is one based and closed, bed and int4range here are zero based and half open
        begin = int(field[3]) - 1
        end = int(field[4])
        strand = field[6] if field[6] in ("+", "-") else None
        regulatory_class = attribute["regulatory_class"]
        if bed is not None:
            print("\t".join([chrom, str(begin), str(end), strand or ".", regulatory_class]), file=bed)
        yield {"chrom": chrom, "region": Range(begin, end, bounds="[)"), "strand": strand,
               "regulatory_class": regulatory_class, "feature_id": attribute.get("ID", ""),
               "note": attribute.get("Note")}


def load(engine, path: str, chrom_of: Dict[str, str], bed: Optional[TextIO]) -> None:
    stored = 0
    chunk: List[Dict[str, Any]] = []
    with Session(engine) as session:
        for row in regulatory_rows(path, chrom_of, bed):
            chunk.append(row)
            if len(chunk) < CHUNK:
                continue
            upsert_many(session, RegulatoryRegion, chunk, conflict_on=["feature_id"])
            session.commit()
            stored += len(chunk)
            chunk = []
        if chunk:
            upsert_many(session, RegulatoryRegion, chunk, conflict_on=["feature_id"])
            session.commit()
            stored += len(chunk)
    report = f"{stored} regulatory regions stored"
    print(report)


#########################################
def main():
    if len(sys.argv) not in (1, 3) or (len(sys.argv) == 3 and sys.argv[1] != "--bed"):
        usage = f"usage: {sys.argv[0]} [--bed <output.bed>]"
        sys.exit(usage)
    for path in (GFF, ASSEMBLY_REPORT):
        if not os.path.exists(path):
            errmsg = f"{path} not found - see the download commands in the docstring of this script"
            sys.exit(errmsg)

    chrom_of = gff.chromosome_names(ASSEMBLY_REPORT)
    engine = db_engine()
    if len(sys.argv) == 3:
        with open(sys.argv[2], "w") as bed:
            load(engine, GFF, chrom_of, bed)
    else:
        load(engine, GFF, chrom_of, None)


#########################################
if __name__ == '__main__':
    main()
