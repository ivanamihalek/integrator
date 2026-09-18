#! /usr/bin/python3
"""Store the transcription start sites of GRCh37 in 'ncbi_tss'.

usage:  05_tss_from_gff.py [--bed <output.bed>]

The RefSeq functional elements documentation promises transcription_start_site features, but the
gff holds none - so the site is taken from the transcript itself: its start on the plus strand,
its end on the minus one.

Selection: transcript features (gff ID 'rna-...') on assembled molecules, carrying an HGNC
cross-reference that names a gene in 'genes'. The HGNC id, not the gene= attribute, is what the
row is attached by, so a selected transcript resolves by construction.

The downloads: see 03_regulatory_regions_from_gff.py.
"""

import os
import sys
from typing import Any, Dict, Iterator, List, Optional, TextIO

from sqlmodel import Session, select

from integrator_utils.python.db import db_engine, require_genes, upsert_many
from integrator_utils.python.models_core import Gene

import gff
from ncbi_models import TranscriptionStartSite

HG19 = "/storage/databases/ncbi/hg19"
GFF = os.path.join(HG19, "GCF_000001405.25_GRCh37.p13_genomic.gff.gz")
ASSEMBLY_REPORT = os.path.join(HG19, "GCF_000001405.25_GRCh37.p13_assembly_report.txt")
CHUNK = 5000


#########################################
def hgnc_id_from(dbxref: str) -> Optional[str]:
    for xref in dbxref.split(","):
        if xref.startswith("HGNC:HGNC:"):
            return xref[len("HGNC:"):]
    return None


def tss_rows(path: str, chrom_of: Dict[str, str], gene_id_of: Dict[str, int],
             bed: Optional[TextIO]) -> Iterator[Dict[str, Any]]:
    for field in gff.lines(path):
        chrom = chrom_of.get(field[0])
        if chrom is None or "transcript_id=" not in field[8]:
            continue
        attribute = gff.attributes(field[8])
        if not attribute.get("ID", "").startswith("rna-"):
            continue  # exons carry a transcript_id too, and they are not transcripts
        hgnc_id = hgnc_id_from(attribute.get("Dbxref", ""))
        if hgnc_id is None or hgnc_id not in gene_id_of:
            continue  # outside the selection: no HGNC cross-reference, or a gene HGNC has withdrawn
        strand = field[6]
        if strand not in ("+", "-"):
            continue
        # gff is one based and closed; bed and the stored position are zero based
        position = int(field[3]) - 1 if strand == "+" else int(field[4]) - 1
        transcript_id = attribute["transcript_id"]
        if bed is not None:
            print("\t".join([chrom, str(position), str(position + 1), transcript_id, strand]), file=bed)
        yield {"gene_id": gene_id_of[hgnc_id], "transcript_id": transcript_id, "transcript_type": field[2],
               "chrom": chrom, "position": position, "strand": strand}


def load(engine, path: str, chrom_of: Dict[str, str], bed: Optional[TextIO]) -> None:
    with Session(engine) as session:
        gene_id_of = {gene.hgnc_id: gene.id for gene in session.exec(select(Gene))}
        stored = 0
        chunk: List[Dict[str, Any]] = []
        conflict_on = ["transcript_id", "chrom", "position"]
        for row in tss_rows(path, chrom_of, gene_id_of, bed):
            chunk.append(row)
            if len(chunk) < CHUNK:
                continue
            upsert_many(session, TranscriptionStartSite, chunk, conflict_on=conflict_on)
            session.commit()
            stored += len(chunk)
            chunk = []
        if chunk:
            upsert_many(session, TranscriptionStartSite, chunk, conflict_on=conflict_on)
            session.commit()
            stored += len(chunk)
    report = f"{stored} transcription start sites stored"
    print(report)


#########################################
def main():
    if len(sys.argv) not in (1, 3) or (len(sys.argv) == 3 and sys.argv[1] != "--bed"):
        usage = f"usage: {sys.argv[0]} [--bed <output.bed>]"
        sys.exit(usage)
    for path in (GFF, ASSEMBLY_REPORT):
        if not os.path.exists(path):
            errmsg = f"{path} not found - see the download commands in 03_regulatory_regions_from_gff.py"
            sys.exit(errmsg)

    engine = db_engine()
    require_genes(engine)
    chrom_of = gff.chromosome_names(ASSEMBLY_REPORT)
    if len(sys.argv) == 3:
        with open(sys.argv[2], "w") as bed:
            load(engine, GFF, chrom_of, bed)
    else:
        load(engine, GFF, chrom_of, None)


#########################################
if __name__ == '__main__':
    main()
