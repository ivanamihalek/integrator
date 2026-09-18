#! /usr/bin/python3
"""Fill the reactome tables: pathways, their hierarchy, and which gene takes part in which.

usage:  16_pathways_from_reactome.py [NCBI2Reactome_All_Levels.txt [ReactomePathwaysRelation.txt]]

Selection: human rows ('Homo sapiens', stable ids 'R-HSA-...') whose NCBI gene id is in 'genes'.
Reactome gives the gene as an Entrez id, so no symbol is involved. Ids not in 'genes' are
counted in the report.

The downloads, both from https://reactome.org/download/current/ :
  NCBI2Reactome_All_Levels.txt
  ReactomePathwaysRelation.txt
"""

import os
import sys
from typing import Any, Dict, Iterator, List, Set, Tuple

from sqlmodel import Session, select

from integrator_utils.python.db import db_engine, require_genes, upsert_many
from integrator_utils.python.models_core import Gene

from reactome_models import ReactomeGenePathway, ReactomePathway, ReactomePathwayRelation

REACTOME = "/storage/databases/reactome"
NCBI2REACTOME = os.path.join(REACTOME, "NCBI2Reactome_All_Levels.txt")
RELATIONS = os.path.join(REACTOME, "ReactomePathwaysRelation.txt")
HUMAN_PREFIX = "R-HSA-"
SPECIES = "Homo sapiens"


#########################################
def participations(path: str) -> Iterator[Tuple[int, str, str, str]]:
    """(entrez gene id, pathway stable id, pathway name, evidence code) for the human rows."""
    with open(path, encoding="utf-8", errors="replace") as inf:
        for line in inf:
            field = line.rstrip("\n").split("\t")
            if len(field) < 6 or field[5] != SPECIES:
                continue
            if not field[0].isdigit() or not field[1].startswith(HUMAN_PREFIX):
                continue
            yield int(field[0]), field[1], field[3].strip(), field[4].strip()


def relations(path: str) -> Iterator[Tuple[str, str]]:
    with open(path) as inf:
        for line in inf:
            field = line.rstrip("\n").split("\t")
            if len(field) == 2 and field[0].startswith(HUMAN_PREFIX) and field[1].startswith(HUMAN_PREFIX):
                yield field[0], field[1]


#########################################
def load(engine, participation_path: str, relation_path: str) -> None:
    pathway_name: Dict[str, str] = {}
    parsed: List[Tuple[int, str, str]] = []
    for entrez_id, stable_id, name, evidence in participations(participation_path):
        pathway_name.setdefault(stable_id, name)
        parsed.append((entrez_id, stable_id, evidence))
    pairs = list(relations(relation_path))
    # a pathway may appear in the hierarchy without taking part in the participation file
    for parent, child in pairs:
        pathway_name.setdefault(parent, "")
        pathway_name.setdefault(child, "")

    with Session(engine) as session:
        pathways = [{"stable_id": stable_id, "name": name} for stable_id, name in pathway_name.items()]
        upsert_many(session, ReactomePathway, pathways, conflict_on=["stable_id"])
        session.commit()

        pathway_id_of = {pathway.stable_id: pathway.id for pathway in session.exec(select(ReactomePathway))}
        gene_id_of = {gene.ncbi_gene_id: gene.id for gene in session.exec(select(Gene)) if gene.ncbi_gene_id}

        edges = [{"parent_id": pathway_id_of[parent], "child_id": pathway_id_of[child]} for parent, child in pairs]
        upsert_many(session, ReactomePathwayRelation, edges, conflict_on=["parent_id", "child_id"])
        session.commit()

        unknown_genes: Set[int] = set()
        rows: Dict[str, Dict[str, Any]] = {}
        for entrez_id, stable_id, evidence in parsed:
            gene_id = gene_id_of.get(entrez_id)
            if gene_id is None:
                unknown_genes.add(entrez_id)
                continue
            pathway_id = pathway_id_of[stable_id]
            # TAS - traceable author statement - is the curated one, and wins over an inference
            key = f"{gene_id}|{pathway_id}"
            if key not in rows or evidence == "TAS":
                rows[key] = {"gene_id": gene_id, "pathway_id": pathway_id, "evidence_code": evidence}
        upsert_many(session, ReactomeGenePathway, list(rows.values()), conflict_on=["gene_id", "pathway_id"])
        session.commit()

    report = f"{len(pathways)} reactome pathways, {len(edges)} hierarchy edges, {len(rows)} participations stored"
    print(report)
    skipped = f"{len(unknown_genes)} reactome gene id(s) are not in 'genes'"
    print(skipped)


#########################################
def main():
    if len(sys.argv) > 3:
        usage = f"usage: {sys.argv[0]} [NCBI2Reactome_All_Levels.txt [ReactomePathwaysRelation.txt]]"
        sys.exit(usage)
    participation_path = sys.argv[1] if len(sys.argv) > 1 else NCBI2REACTOME
    relation_path = sys.argv[2] if len(sys.argv) > 2 else RELATIONS
    for path in (participation_path, relation_path):
        if not os.path.exists(path):
            errmsg = f"{path} not found - see the download urls in the docstring of this script"
            sys.exit(errmsg)

    engine = db_engine()
    require_genes(engine)
    load(engine, participation_path, relation_path)


#########################################
if __name__ == '__main__':
    main()
