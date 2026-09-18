#! /usr/bin/python3
"""Tables filled from the Reactome pathway downloads."""

from typing import Optional

from sqlmodel import Field, SQLModel


#########################################
class ReactomePathway(SQLModel, table=True):
    __tablename__ = "reactome_pathways"

    id: Optional[int] = Field(default=None, primary_key=True)
    stable_id: str = Field(index=True, unique=True)            # "R-HSA-109582"
    name: str = Field(index=True)                              # "Hemostasis"


class ReactomePathwayRelation(SQLModel, table=True):
    """The pathway hierarchy, as parent-child pairs.

    The legacy script built the tree in memory and stored it flattened; kept as edges here, so
    that a recursive query can walk it in either direction.
    """

    __tablename__ = "reactome_pathway_relations"

    parent_id: int = Field(foreign_key="reactome_pathways.id", primary_key=True, index=True)
    child_id: int = Field(foreign_key="reactome_pathways.id", primary_key=True, index=True)


class ReactomeGenePathway(SQLModel, table=True):
    """Which gene takes part in which pathway, attached by the NCBI gene id reactome provides."""

    __tablename__ = "reactome_gene_pathways"

    gene_id: int = Field(foreign_key="genes.id", primary_key=True, index=True)
    pathway_id: int = Field(foreign_key="reactome_pathways.id", primary_key=True, index=True)
    evidence_code: Optional[str] = Field(default=None, index=True)   # TAS (asserted), IEA (inferred)
