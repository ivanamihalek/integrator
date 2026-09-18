#! /usr/bin/python3
"""Fill 'uniprot_basic_infos' from the uniprot (SwissProt) flat file.

usage:  01_uniprot_basic_info.py [uniprot_sprot.dat]

Selection: reviewed entries carrying an HGNC cross-reference that names a gene in 'genes'. That
cross-reference, not the GN Name= line, is what the row is attached by - uniprot and HGNC do not
always agree on the symbol, and the flat file on disk is not of the same vintage as the HGNC
download. Entries whose HGNC id is no longer in 'genes' are counted in the report and skipped;
an entry selected but unresolvable stops the run.

This replaces 01_uniprot_parser.pl, which wrote a tsv for mysqlimport. Two differences worth
knowing: the canonical sequence is taken straight from the SQ block instead of from a blast
database, and the ';'-joined columns are stored as arrays.

The source: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/
            complete/uniprot_sprot.dat.gz
"""

import gzip
import os
import re
import sys
from typing import Any, Dict, Iterator, List, Optional

from sqlmodel import Session, select

from integrator_utils.python.db import db_engine, require_genes, upsert_many
from integrator_utils.python.models_core import Gene, UniprotBasicInfo

UNIPROT_SPROT = "/storage/databases/uniprot/uniprot_sprot.dat"
CHUNK = 2000

HGNC_PATTERN = re.compile(r"^DR   HGNC; (HGNC:\d+);")
ENSEMBL_GENE_PATTERN = re.compile(r"(ENSG\d{11})")
EC_PATTERN = re.compile(r"EC=([\d\.\-]+)")
AA_LENGTH_PATTERN = re.compile(r"^SQ   SEQUENCE\s+(\d+)\s+AA")
COFACTOR_NAME_PATTERN = re.compile(r"Name=(.+?);")
EVIDENCE_PATTERN = re.compile(r"\s*\{[^{}]*\}")


#########################################
def entries(path: str) -> Iterator[List[str]]:
    """The flat file as blocks of lines, one block per entry."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as inf:
        block: List[str] = []
        for line in inf:
            if line.startswith("//"):
                if block:
                    yield block
                block = []
            else:
                block.append(line.rstrip("\n"))
        if block:
            yield block


def strip_evidence(text: str) -> str:
    """Drop the '{ECO:...}' evidence tags uniprot scatters through its free text."""
    return re.sub(r"\s+", " ", EVIDENCE_PATTERN.sub("", text)).strip()


#########################################
def parse_entry(block: List[str]) -> Optional[Dict[str, Any]]:
    """Everything 'uniprot_basic_infos' wants, or None if the entry has no HGNC cross-reference."""
    accessions: List[str] = []
    hgnc_ids: List[str] = []
    ensembl_gene_ids: List[str] = []
    cofactors: List[str] = []
    tissues: List[str] = []
    full_name: Optional[str] = None
    ec_number: Optional[str] = None
    aa_length: Optional[int] = None
    function: List[str] = []
    location: List[str] = []
    sequence: List[str] = []
    reading: Optional[str] = None

    for line in block:
        if line.startswith("AC   "):
            accessions += [item.strip() for item in line[5:].split(";") if item.strip()]
        elif line.startswith("DE   RecName: Full="):
            if full_name is None:
                full_name = strip_evidence(line[len("DE   RecName: Full="):]).rstrip(";")
        elif line.startswith("DE ") and "EC=" in line and ec_number is None:
            ec_match = EC_PATTERN.search(line)
            ec_number = ec_match.group(1) if ec_match else None
        elif line.startswith("DR   HGNC; "):
            hgnc_match = HGNC_PATTERN.match(line)
            if hgnc_match:
                # 80 entries name two genes - readthrough transcripts and the like. The row is
                # attached to the first, which is the one uniprot names the entry after
                hgnc_ids.append(hgnc_match.group(1))
        elif line.startswith("DR   Ensembl;"):
            for ensembl_gene_id in ENSEMBL_GENE_PATTERN.findall(line):
                if ensembl_gene_id not in ensembl_gene_ids:
                    ensembl_gene_ids.append(ensembl_gene_id)
        elif line.startswith("RC   TISSUE="):
            for tissue in strip_evidence(line[len("RC   TISSUE="):]).rstrip(";,").split(","):
                if tissue.strip() and tissue.strip() not in tissues:
                    tissues.append(tissue.strip())
        elif line.startswith("CC   -!- FUNCTION:"):
            reading = "function"
            function.append(line[len("CC   -!- FUNCTION:"):].strip())
        elif line.startswith("CC   -!- SUBCELLULAR LOCATION:"):
            reading = "location"
            location.append(line[len("CC   -!- SUBCELLULAR LOCATION:"):].strip())
        elif line.startswith("CC   -!- COFACTOR:"):
            reading = "cofactor"
        elif line.startswith("CC   -!-") or line.startswith("CC   ---"):
            reading = None
        elif line.startswith("CC       ") and reading == "function":
            function.append(line.strip())
        elif line.startswith("CC       ") and reading == "location":
            location.append(line.strip())
        elif line.startswith("CC       ") and reading == "cofactor":
            for name in COFACTOR_NAME_PATTERN.findall(line):
                if name not in cofactors:
                    cofactors.append(name)
        elif line.startswith("SQ   SEQUENCE"):
            reading = "sequence"
            length_match = AA_LENGTH_PATTERN.match(line)
            aa_length = int(length_match.group(1)) if length_match else None
        elif line.startswith("     ") and reading == "sequence":
            sequence.append(line.replace(" ", ""))

    if not hgnc_ids or not accessions:
        return None
    return {"uniprot_id": accessions[0], "hgnc_ids": hgnc_ids, "full_name": full_name,
            "ec_number": ec_number, "canonical_aa_length": aa_length, "sequence": "".join(sequence) or None,
            "subcellular_location": strip_evidence(" ".join(location)) or None,
            "function": strip_evidence(" ".join(function)) or None,
            "ensembl_gene_ids": ensembl_gene_ids or None, "cofactors": cofactors or None,
            "tissues": tissues or None, "secondary_ids": accessions[1:] or None}


#########################################
def load(engine, path: str) -> None:
    stored = 0
    withdrawn = 0
    multi_gene = 0
    chunk: List[Dict[str, Any]] = []
    with Session(engine) as session:
        gene_id_of = {gene.hgnc_id: gene.id for gene in session.exec(select(Gene))}
        for block in entries(path):
            parsed = parse_entry(block)
            if parsed is None:
                continue
            hgnc_ids = parsed.pop("hgnc_ids")
            multi_gene += len(hgnc_ids) > 1
            gene_id = gene_id_of.get(hgnc_ids[0])
            if gene_id is None:  # an HGNC entry withdrawn since this flat file was cut
                withdrawn += 1
                continue
            parsed["gene_id"] = gene_id
            chunk.append(parsed)
            if len(chunk) < CHUNK:
                continue
            upsert_many(session, UniprotBasicInfo, chunk, conflict_on=["uniprot_id"])
            session.commit()
            stored += len(chunk)
            chunk = []
            progress = f"\t {stored} entries"
            print(progress)
        if chunk:
            upsert_many(session, UniprotBasicInfo, chunk, conflict_on=["uniprot_id"])
            session.commit()
            stored += len(chunk)

    report = f"{stored} uniprot entries stored"
    print(report)
    skipped = f"{withdrawn} entry(s) refer to an hgnc id that is not in 'genes' (withdrawn since)"
    print(skipped)
    shared = f"{multi_gene} entry(s) name more than one gene - attached to the first"
    print(shared)


#########################################
def main():
    if len(sys.argv) > 2:
        usage = f"usage: {sys.argv[0]} [uniprot_sprot.dat]"
        sys.exit(usage)
    path = sys.argv[1] if len(sys.argv) == 2 else UNIPROT_SPROT
    if not os.path.exists(path):
        errmsg = f"{path} not found - see the source url in the docstring of this script"
        sys.exit(errmsg)

    engine = db_engine()
    require_genes(engine)
    load(engine, path)


#########################################
if __name__ == '__main__':
    main()
