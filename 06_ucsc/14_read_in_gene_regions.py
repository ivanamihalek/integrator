#! /usr/bin/python3
"""Store the UCSC refGene exon and intron regions in 'ucsc_refgene_regions'.

usage:  14_read_in_gene_regions.py [hg19]

Reads the per-chromosome bed files that 13_get_gene_regions_from_ucsc.py writes out of UCSC's
public mysql server, so this script needs no network of its own.

Selection, in two parts:

  * primary assembly only - the _hap, _random and chrUn files repeat the same transcripts in
    coordinates nothing downstream uses, and the MHC haplotypes would multiply them sevenfold;
  * transcripts that name a gene this database knows - resolved first by RefSeq accession
    against 'ncbi_tss', and failing that by the symbol in the fifth bed column.

UCSC's refGene carries some 1,600 names HGNC has never approved: LOC provisional names, clone
names (DKFZp..., FLJ...), readthrough fusions, and symbols HGNC has withdrawn. Those name no
gene here, so their transcripts cannot be stored. They are written to a report file and the
script exits non-zero after loading everything that did resolve - the load is idempotent, so
re-running once they are dealt with costs nothing.
"""

import os
import sys
from typing import Any, Dict, Iterator, List, Optional, Set, Tuple

from sqlalchemy import text
from sqlalchemy.dialects.postgresql import Range
from sqlmodel import Session

from integrator_utils.python.db import db_engine, gene_id_map, require_genes, upsert_many

from ucsc_models import RefgeneRegion

GENE_REGIONS = "/storage/databases/ucsc/gene_regions"
UNRESOLVED_REPORT = "/tmp/ucsc_unresolved_transcripts.tsv"
CHUNK = 5000
CONFLICT_ON = ["assembly", "transcript_id", "chrom", "feature", "region"]


#########################################
def primary_assembly_files(directory: str) -> List[str]:
    """chr1.bed ... chrY.bed, leaving out the alt haplotypes, the randoms and chrUn."""
    paths: List[str] = []
    for entry in sorted(os.scandir(directory), key=lambda item: item.name):
        if not entry.name.endswith(".bed"):
            continue
        chrom = entry.name[:-len(".bed")]
        if "_" in chrom:
            continue
        paths.append(entry.path)
    return paths


def bed_rows(path: str) -> Iterator[Tuple[str, int, int, str, str, str, str]]:
    with open(path) as inf:
        for line in inf:
            field = line.rstrip("\n").split("\t")
            if len(field) < 7:
                continue
            yield field[0], int(field[1]), int(field[2]), field[3], field[4], field[5], field[6]


#########################################
def transcript_gene_ids(session: Session) -> Dict[str, int]:
    """RefSeq accession (unversioned) -> gene id, from the tss table 05_ncbi fills.

    Read as sql rather than through the model: 05_ncbi is not an importable package, and the
    dependency here is on the table, not on that directory's code.
    """
    statement = text("select split_part(transcript_id, \'.\', 1), gene_id from ncbi_tss")
    return {transcript_id: gene_id for transcript_id, gene_id in session.exec(statement)}


def load(engine, assembly: str) -> None:
    directory = os.path.join(GENE_REGIONS, assembly)
    if not os.path.isdir(directory):
        errmsg = f"{directory} not found - run 13_get_gene_regions_from_ucsc.py first"
        sys.exit(errmsg)

    stored = 0
    unresolved: Dict[str, str] = {}
    chunk: List[Dict[str, Any]] = []
    with Session(engine) as session:
        gene_id_of_transcript = transcript_gene_ids(session)
        gene_id_of_symbol = gene_id_map(session)
        for path in primary_assembly_files(directory):
            for chrom, begin, end, transcript_id, symbol, strand, feature in bed_rows(path):
                gene_id = gene_id_of_transcript.get(transcript_id) or gene_id_of_symbol.get(symbol.lower())
                if gene_id is None:
                    unresolved[transcript_id] = symbol
                    continue
                chunk.append({"gene_id": gene_id, "assembly": assembly, "transcript_id": transcript_id,
                              "chrom": chrom, "region": Range(begin, end, bounds="[)"), "strand": strand,
                              "feature": feature})
                if len(chunk) < CHUNK:
                    continue
                upsert_many(session, RefgeneRegion, chunk, conflict_on=CONFLICT_ON)
                session.commit()
                stored += len(chunk)
                chunk = []
            progress = f"\t {os.path.basename(path)}: {stored} regions"
            print(progress)
        if chunk:
            upsert_many(session, RefgeneRegion, chunk, conflict_on=CONFLICT_ON)
            session.commit()
            stored += len(chunk)

    report = f"{stored} refgene regions stored for {assembly}"
    print(report)
    if not unresolved:
        return
    with open(UNRESOLVED_REPORT, "w") as outf:
        for transcript_id, symbol in sorted(unresolved.items()):
            print(f"{transcript_id}\t{symbol}", file=outf)
    errmsg = f"{len(unresolved)} transcript(s) name no gene in 'genes' - listed in {UNRESOLVED_REPORT}"
    sys.exit(errmsg)


#########################################
def main():
    if len(sys.argv) > 2:
        usage = f"usage: {sys.argv[0]} [hg19]"
        sys.exit(usage)
    assembly = sys.argv[1] if len(sys.argv) == 2 else "hg19"

    engine = db_engine()
    require_genes(engine)
    load(engine, assembly)


#########################################
if __name__ == '__main__':
    main()
