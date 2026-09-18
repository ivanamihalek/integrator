#! /usr/bin/python3
"""Fill 'kegg_pathways' and 'kegg_gene_pathways' from the KEGG pathway downloads.

usage:  18_kegg_id2pthwy.py [pathway_names.txt [pathway2gene.tsv]]

Selection: pathway memberships whose KEGG gene id is in 'genes'. KEGG writes those ids as
'hsa:10327', and the number is an Entrez gene id, so genes are attached by ncbi_gene_id rather
than by any symbol. Ids not in 'genes' - KEGG still carries genes HGNC has retired - are
counted in the report and dropped.

The legacy tables this replaces, 'kegg_human' and 'kegg_pathway_name', kept the pathway list as
a semicolon-joined text column on the gene row; the membership is a table of its own here.

The downloads:
  https://rest.kegg.jp/list/pathway/hsa     -> pathway_list_hsa.tsv
  https://rest.kegg.jp/link/hsa/pathway     -> pathway2gene_current.tsv

The older brite listing (pathway_names.txt) is read too, and is the only one of the two that
carries the Metabolism / Carbohydrate metabolism hierarchy - but it has to be current, since a
pathway in the link file with no name stops this script.
"""

import os
import re
import sys
from typing import Any, Dict, List, Optional, Tuple

from sqlmodel import Session, delete, select

from integrator_utils.python.db import db_engine, require_genes, upsert_many
from integrator_utils.python.models_core import Gene

from kegg_models import KeggGenePathway, KeggPathway

PATHWAY_NAMES = "/storage/databases/kegg/pathway_list_hsa.tsv"
PATHWAY2GENE = "/storage/databases/kegg/pathway2gene_current.tsv"

# the rest api writes "hsa00010\tGlycolysis / Gluconeogenesis - Homo sapiens (human)"
FLAT_PATTERN = re.compile(r"^hsa(\d{5})\t(.+?)\s*$")
# the brite listing writes "     00010  Glycolysis / Gluconeogenesis", indentation as hierarchy
PATHWAY_PATTERN = re.compile(r"^(\s+)(\d{5})\s+(.+?)\s*$")
HEADING_PATTERN = re.compile(r"^(\s+)(\D.*?)\s*$")
SPECIES_SUFFIX = " - Homo sapiens (human)"


#########################################
def pathways(path: str) -> List[Dict[str, Any]]:
    """Either format the pathway list comes in: the flat rest listing, or the brite hierarchy."""
    category: Optional[str] = None
    subcategory: Optional[str] = None
    rows: List[Dict[str, Any]] = []
    with open(path) as inf:
        for line in inf:
            flat_match = FLAT_PATTERN.match(line)
            if flat_match:
                name = flat_match.group(2).removesuffix(SPECIES_SUFFIX)
                rows.append({"kegg_pathway_id": f"hsa{flat_match.group(1)}", "name": name,
                             "category": None, "subcategory": None})
                continue
            pathway_match = PATHWAY_PATTERN.match(line)
            if pathway_match:
                rows.append({"kegg_pathway_id": f"hsa{pathway_match.group(2)}", "name": pathway_match.group(3),
                             "category": category, "subcategory": subcategory})
                continue
            heading_match = HEADING_PATTERN.match(line)
            if not heading_match:
                continue
            if len(heading_match.group(1)) <= 1:
                category = heading_match.group(2)
                subcategory = None
            else:
                subcategory = heading_match.group(2)
    return rows


def memberships(path: str) -> List[Tuple[str, int]]:
    """(kegg pathway id, entrez gene id) for every line of the link file."""
    rows: List[Tuple[str, int]] = []
    with open(path) as inf:
        for line in inf:
            field = line.split()
            if len(field) != 2:
                continue
            pathway = field[0].replace("path:", "")
            gene = field[1].replace("hsa:", "")
            if gene.isdigit():
                rows.append((pathway, int(gene)))
    return rows


#########################################
def load(engine, names_path: str, links_path: str) -> None:
    pathway_rows = pathways(names_path)
    links = memberships(links_path)
    with Session(engine) as session:
        upsert_many(session, KeggPathway, pathway_rows, conflict_on=["kegg_pathway_id"])
        session.commit()

        # the file is the authority on which pathways exist, so the list is replaced rather than
        # merged: a pathway KEGG has retired since the last run goes, and its memberships with it
        current = [row["kegg_pathway_id"] for row in pathway_rows]
        retired = session.exec(select(KeggPathway).where(KeggPathway.kegg_pathway_id.not_in(current))).all()
        if retired:
            retired_ids = [pathway.id for pathway in retired]
            session.exec(delete(KeggGenePathway).where(KeggGenePathway.pathway_id.in_(retired_ids)))
            session.exec(delete(KeggPathway).where(KeggPathway.id.in_(retired_ids)))
            session.commit()
            dropped = f"{len(retired)} pathway(s) no longer in the source - dropped"
            print(dropped)

        pathway_id_of = {pathway.kegg_pathway_id: pathway.id for pathway in session.exec(select(KeggPathway))}
        gene_id_of = {gene.ncbi_gene_id: gene.id for gene in session.exec(select(Gene)) if gene.ncbi_gene_id}

        unknown_genes = set()
        unknown_pathways = set()
        rows: Dict[str, Dict[str, int]] = {}
        for kegg_pathway_id, entrez_id in links:
            pathway_id = pathway_id_of.get(kegg_pathway_id)
            gene_id = gene_id_of.get(entrez_id)
            if pathway_id is None:
                unknown_pathways.add(kegg_pathway_id)
            elif gene_id is None:
                unknown_genes.add(entrez_id)
            else:
                rows[f"{gene_id}|{pathway_id}"] = {"gene_id": gene_id, "pathway_id": pathway_id}
        # the link file is a complete snapshot, so the memberships of the pathways it names are
        # replaced: a gene KEGG has taken out of a pathway has to come out of the table too
        session.exec(delete(KeggGenePathway).where(KeggGenePathway.pathway_id.in_(list(pathway_id_of.values()))))
        upsert_many(session, KeggGenePathway, list(rows.values()), conflict_on=["gene_id", "pathway_id"])
        session.commit()

    report = f"{len(pathway_rows)} kegg pathways and {len(rows)} memberships stored, from {len(links)} links"
    print(report)
    skipped = f"{len(unknown_genes)} kegg gene id(s) are not in 'genes'"
    print(skipped)
    if unknown_pathways:
        # the link file is the authority on which pathways exist; a name missing for one of them
        # means the two downloads are out of step, and the pathway list has to be refreshed
        missing = ", ".join(sorted(unknown_pathways))
        errmsg = f"{len(unknown_pathways)} pathway(s) in the link file have no name: {missing}"
        sys.exit(errmsg)


#########################################
def main():
    if len(sys.argv) > 3:
        usage = f"usage: {sys.argv[0]} [pathway_names.txt [pathway2gene.tsv]]"
        sys.exit(usage)
    names_path = sys.argv[1] if len(sys.argv) > 1 else PATHWAY_NAMES
    links_path = sys.argv[2] if len(sys.argv) > 2 else PATHWAY2GENE
    for path in (names_path, links_path):
        if not os.path.exists(path):
            errmsg = f"{path} not found - see the download urls in the docstring of this script"
            sys.exit(errmsg)

    engine = db_engine()
    require_genes(engine)
    load(engine, names_path, links_path)


#########################################
if __name__ == '__main__':
    main()
