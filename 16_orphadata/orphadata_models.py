#! /usr/bin/python3
"""Tables filled from the Orphanet rare disease nomenclature (orphadata xml products)."""

from typing import Optional

from sqlmodel import Field, SQLModel


#########################################
class OrphaDisease(SQLModel, table=True):
    __tablename__ = "orpha_diseases"

    id: Optional[int] = Field(default=None, primary_key=True)
    orpha_code: int = Field(index=True, unique=True)
    name: str = Field(index=True)
    disorder_type: Optional[str] = Field(default=None, index=True)   # Disease, Malformation syndrome, ...
    disorder_group: Optional[str] = None                             # Disorder, Group of disorders, Subtype


class OrphaDiseaseXref(SQLModel, table=True):
    """Orphanet's cross-references: OMIM, ICD-10/11, UMLS, MeSH, MONDO, GARD, MedDRA."""

    __tablename__ = "orpha_disease_xrefs"

    disease_id: int = Field(foreign_key="orpha_diseases.id", primary_key=True, index=True)
    source: str = Field(primary_key=True, index=True)
    reference: str = Field(primary_key=True, index=True)
    relation: Optional[str] = None                                   # "E" (exact), "NTBT", "BTNT", ...


class OrphaDiseaseGene(SQLModel, table=True):
    """Disease-gene associations, attached by the HGNC id orphanet carries for each gene."""

    __tablename__ = "orpha_disease_genes"

    disease_id: int = Field(foreign_key="orpha_diseases.id", primary_key=True, index=True)
    gene_id: int = Field(foreign_key="genes.id", primary_key=True, index=True)
    # part of the key: orphanet records the same disease and gene under several association
    # types - a germline mutation in one entry, a modifying mutation in another
    association_type: str = Field(primary_key=True, index=True)
    association_status: Optional[str] = None
    source_of_validation: Optional[str] = None                       # pubmed ids, as orphanet writes them
