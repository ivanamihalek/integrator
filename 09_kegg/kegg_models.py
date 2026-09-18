#! /usr/bin/python3
"""Tables filled from the KEGG pathway downloads."""

from typing import Optional

from sqlmodel import Field, SQLModel


#########################################
class KeggPathway(SQLModel, table=True):
    __tablename__ = "kegg_pathways"

    id: Optional[int] = Field(default=None, primary_key=True)
    kegg_pathway_id: str = Field(index=True, unique=True)      # "hsa00010"
    name: str = Field(index=True)                              # "Glycolysis / Gluconeogenesis"
    category: Optional[str] = Field(default=None, index=True)  # "Metabolism"
    subcategory: Optional[str] = None                          # "Carbohydrate metabolism"


class KeggGenePathway(SQLModel, table=True):
    """Which gene is in which pathway. KEGG's 'hsa:' gene ids are Entrez ids, so the genes are
    attached by 'genes.ncbi_gene_id' and no symbol is involved."""

    __tablename__ = "kegg_gene_pathways"

    gene_id: int = Field(foreign_key="genes.id", primary_key=True, index=True)
    pathway_id: int = Field(foreign_key="kegg_pathways.id", primary_key=True, index=True)
