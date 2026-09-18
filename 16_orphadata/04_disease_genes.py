#! /usr/bin/python3
"""Fill 'orpha_disease_genes' from the orphadata disease-gene association product.

usage:  04_disease_genes.py [en_product6.xml]

Selection: associations whose gene is an actual gene - orphanet's GeneType says 'gene with
protein product' or 'Non-coding RNA' - carrying an HGNC cross-reference. Orphanet writes that
reference as a bare number ('30497'), so it is read as 'HGNC:30497'.

The GeneType test is not cosmetic: orphanet also records 36 associations against a
'Disorder-associated locus' (USH1E, DYT13, SPG14 and the like), which have an HGNC id but are
phenotype loci, not approved genes, and so have no row in 'genes' and never will. With those
excluded, and ids HGNC has withdrawn since orphanet's cut (two snoRNA clusters at the time of
writing) left out, any HGNC id that still fails to resolve stops the run.

02_diseases_and_xrefs.py has to have run first: the diseases named here are looked up, not made.

The download:
  wget -O /storage/databases/orphadata/en_product6.xml https://www.orphadata.com/data/xml/en_product6.xml
"""

import os
import sys
import xml.etree.ElementTree as ElementTree
from typing import Any, Dict, List, Optional, Tuple

from sqlmodel import Session, select

from integrator_utils.python.db import db_engine, require_genes, upsert_many, withdrawn_hgnc_ids
from integrator_utils.python.models_core import Gene

from orphadata_models import OrphaDisease, OrphaDiseaseGene

PRODUCT6 = "/storage/databases/orphadata/en_product6.xml"
# what orphanet calls a gene; the fourth kind, "Disorder-associated locus", is not one
GENE_TYPES = ("gene with protein product", "Non-coding RNA")


#########################################
def text_of(element: Optional[ElementTree.Element], path: str) -> Optional[str]:
    if element is None:
        return None
    found = element.find(path)
    return found.text.strip() if found is not None and found.text else None


def hgnc_id_of(gene: ElementTree.Element) -> Optional[str]:
    for reference in gene.iterfind("ExternalReferenceList/ExternalReference"):
        if text_of(reference, "Source") == "HGNC":
            value = text_of(reference, "Reference")
            return f"HGNC:{value}" if value else None
    return None


def associations(path: str) -> Tuple[List[Dict[str, Any]], int]:
    """The selected associations, and the count of the ones dropped for naming a locus."""
    rows: List[Dict[str, Any]] = []
    loci = 0
    for event, disorder in ElementTree.iterparse(path, events=("end", )):
        if disorder.tag != "Disorder":
            continue
        orpha_code = text_of(disorder, "OrphaCode")
        for association in disorder.iterfind("DisorderGeneAssociationList/DisorderGeneAssociation"):
            gene = association.find("Gene")
            if gene is None or text_of(gene, "GeneType/Name") not in GENE_TYPES:
                loci += 1
                continue
            hgnc_id = hgnc_id_of(gene)
            if orpha_code is None or hgnc_id is None:
                continue
            rows.append({"orpha_code": int(orpha_code), "hgnc_id": hgnc_id,
                         "association_type": text_of(association, "DisorderGeneAssociationType/Name"),
                         "association_status": text_of(association, "DisorderGeneAssociationStatus/Name"),
                         "source_of_validation": text_of(association, "SourceOfValidation")})
        disorder.clear()
    return rows, loci


#########################################
def load(engine, path: str) -> None:
    parsed, loci = associations(path)
    with Session(engine) as session:
        gene_id_of = {gene.hgnc_id: gene.id for gene in session.exec(select(Gene))}
        disease_id_of = {disease.orpha_code: disease.id for disease in session.exec(select(OrphaDisease))}

        # an id HGNC has retired names no current gene: outside the selection, not an error
        retired = withdrawn_hgnc_ids(session)
        withdrawn = sorted({row["hgnc_id"] for row in parsed if row["hgnc_id"] in retired})
        parsed = [row for row in parsed if row["hgnc_id"] not in retired]

        unresolved = sorted({row["hgnc_id"] for row in parsed if row["hgnc_id"] not in gene_id_of})
        if unresolved:
            offenders = ", ".join(unresolved)
            errmsg = f"{len(unresolved)} hgnc id(s) from orphanet are not in 'genes': {offenders}"
            sys.exit(errmsg)
        missing = sorted({row["orpha_code"] for row in parsed if row["orpha_code"] not in disease_id_of})
        if missing:
            codes = ", ".join(str(code) for code in missing[:20])
            errmsg = f"{len(missing)} orpha code(s) are not in 'orpha_diseases' - run 02_diseases_and_xrefs.py: {codes}"
            sys.exit(errmsg)

        # keyed, so that a single statement never touches the same (disease, gene) pair twice
        rows: Dict[str, Dict[str, Any]] = {}
        for row in parsed:
            disease_id = disease_id_of[row["orpha_code"]]
            gene_id = gene_id_of[row["hgnc_id"]]
            association_type = row["association_type"] or "unspecified"
            rows[f"{disease_id}|{gene_id}|{association_type}"] = {
                "disease_id": disease_id, "gene_id": gene_id, "association_type": association_type,
                "association_status": row["association_status"],
                "source_of_validation": row["source_of_validation"]}
        conflict_on = ["disease_id", "gene_id", "association_type"]
        upsert_many(session, OrphaDiseaseGene, list(rows.values()), conflict_on=conflict_on)
        session.commit()

    report = f"{len(rows)} disease-gene associations stored, from {len(parsed)} selected"
    print(report)
    skipped = f"{loci} association(s) name a disorder-associated locus rather than a gene - not stored"
    print(skipped)
    retired_note = f"{len(withdrawn)} gene(s) named here have been withdrawn by HGNC - not stored: {', '.join(withdrawn)}"
    print(retired_note)


#########################################
def main():
    if len(sys.argv) > 2:
        usage = f"usage: {sys.argv[0]} [en_product6.xml]"
        sys.exit(usage)
    path = sys.argv[1] if len(sys.argv) == 2 else PRODUCT6
    if not os.path.exists(path):
        errmsg = f"{path} not found - see the download command in the docstring of this script"
        sys.exit(errmsg)

    engine = db_engine()
    require_genes(engine)
    load(engine, path)


#########################################
if __name__ == '__main__':
    main()
