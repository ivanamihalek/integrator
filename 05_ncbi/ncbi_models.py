#! /usr/bin/python3
"""Tables filled from the NCBI RefSeq annotation of the human genome.

Intervals are int4range with a GiST index rather than a start/end pair: every consumer of these
tables is asking whether something overlaps something else, and '&&' answers that directly.
"""

from typing import Any, Optional

from sqlalchemy.dialects.postgresql import INT4RANGE
from sqlalchemy import UniqueConstraint
from sqlmodel import Field, Index, SQLModel

ASSEMBLY = "GRCh37.p13"


#########################################
class RegulatoryRegion(SQLModel, table=True):
    """RefSeqFE regulatory features: enhancers, silencers, promoters and the like.

    Not linked to genes: the Dbxref of a regulatory feature is the feature's own GeneID, not
    that of the gene it regulates, and the gff says nothing about which gene that would be.
    """

    __tablename__ = "ncbi_regulatory_regions"
    __table_args__ = (Index("ix_ncbi_regulatory_regions_region", "region", postgresql_using="gist"), )

    id: Optional[int] = Field(default=None, primary_key=True)
    chrom: str = Field(index=True)
    region: Any = Field(sa_type=INT4RANGE)          # zero based, half open, like bed
    strand: Optional[str] = None
    regulatory_class: str = Field(index=True)
    feature_id: str = Field(index=True, unique=True)  # the gff ID attribute, "id-GeneID:127898558"
    note: Optional[str] = None


class TranscriptionStartSite(SQLModel, table=True):
    """One row per transcript: the gff has no transcription_start_site features, so the site is
    taken from the transcript itself - its start on the plus strand, its end on the minus one."""

    __tablename__ = "ncbi_tss"
    # a transcript can be placed twice on the primary assembly - the pseudoautosomal regions of
    # X and Y - so the natural key the loader upserts on is the placement, not the transcript
    __table_args__ = (UniqueConstraint("transcript_id", "chrom", "position", name="uq_ncbi_tss_placement"), )

    id: Optional[int] = Field(default=None, primary_key=True)
    gene_id: int = Field(foreign_key="genes.id", index=True)
    transcript_id: str = Field(index=True)          # RefSeq accession, "NM_130786.3"
    transcript_type: str = Field(index=True)        # mRNA, ncRNA, transcript, ...
    chrom: str = Field(index=True)
    position: int                                   # zero based, like bed
    strand: str
