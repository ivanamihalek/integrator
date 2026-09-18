#! /usr/bin/python3
"""Fill 'orpha_diseases' and 'orpha_disease_xrefs' from the orphadata cross-reference product.

usage:  02_diseases_and_xrefs.py [en_product1.xml]

No gene is involved here, so no resolution is needed; 04_disease_genes.py does that part.

The download:
  wget -O /storage/databases/orphadata/en_product1.xml https://www.orphadata.com/data/xml/en_product1.xml
"""

import os
import sys
import xml.etree.ElementTree as ElementTree
from typing import Any, Dict, Iterator, List, Optional, Tuple

from sqlmodel import Session, select

from integrator_utils.python.db import db_engine, upsert_many

from orphadata_models import OrphaDisease, OrphaDiseaseXref

PRODUCT1 = "/storage/databases/orphadata/en_product1.xml"


#########################################
def text_of(element: Optional[ElementTree.Element], path: str) -> Optional[str]:
    if element is None:
        return None
    found = element.find(path)
    return found.text.strip() if found is not None and found.text else None


def disorders(path: str) -> Iterator[ElementTree.Element]:
    """Stream the file: the xml is 54 MB, and only one disorder is needed at a time."""
    for event, element in ElementTree.iterparse(path, events=("end", )):
        if element.tag == "Disorder":
            yield element
            element.clear()


def parse(path: str) -> Tuple[List[Dict[str, Any]], Dict[int, List[Dict[str, str]]]]:
    diseases: List[Dict[str, Any]] = []
    xrefs: Dict[int, List[Dict[str, str]]] = {}
    for disorder in disorders(path):
        orpha_code = text_of(disorder, "OrphaCode")
        name = text_of(disorder, "Name")
        if not orpha_code or not name:
            continue
        disease = {"orpha_code": int(orpha_code), "name": name,
                   "disorder_type": text_of(disorder, "DisorderType/Name"),
                   "disorder_group": text_of(disorder, "DisorderGroup/Name")}
        diseases.append(disease)
        seen: Dict[str, Dict[str, str]] = {}
        for reference in disorder.iterfind("ExternalReferenceList/ExternalReference"):
            source = text_of(reference, "Source")
            value = text_of(reference, "Reference")
            if not source or not value:
                continue
            seen[f"{source}|{value}"] = {"source": source, "reference": value,
                                         "relation": text_of(reference, "DisorderMappingRelation/Name")}
        xrefs[int(orpha_code)] = list(seen.values())
    return diseases, xrefs


#########################################
def load(engine, path: str) -> None:
    diseases, xrefs = parse(path)
    with Session(engine) as session:
        upsert_many(session, OrphaDisease, diseases, conflict_on=["orpha_code"])
        session.commit()

        orpha_codes = [disease["orpha_code"] for disease in diseases]
        disease_id_of = {disease.orpha_code: disease.id
                         for disease in session.exec(select(OrphaDisease).where(OrphaDisease.orpha_code.in_(orpha_codes)))}
        rows: List[Dict[str, Any]] = []
        for orpha_code, references in xrefs.items():
            for reference in references:
                rows.append({"disease_id": disease_id_of[orpha_code], **reference})
        upsert_many(session, OrphaDiseaseXref, rows, conflict_on=["disease_id", "source", "reference"])
        session.commit()

    report = f"{len(diseases)} orphanet diseases and {len(rows)} cross-references stored"
    print(report)


#########################################
def main():
    if len(sys.argv) > 2:
        usage = f"usage: {sys.argv[0]} [en_product1.xml]"
        sys.exit(usage)
    path = sys.argv[1] if len(sys.argv) == 2 else PRODUCT1
    if not os.path.exists(path):
        errmsg = f"{path} not found - see the download command in the docstring of this script"
        sys.exit(errmsg)

    load(db_engine(), path)


#########################################
if __name__ == '__main__':
    main()
