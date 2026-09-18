#! /usr/bin/python3
"""Tables filled from the UCSC refGene annotation."""

from typing import Any, Optional

from sqlalchemy import UniqueConstraint
from sqlalchemy.dialects.postgresql import INT4RANGE
from sqlmodel import Field, Index, SQLModel


#########################################
class RefgeneRegion(SQLModel, table=True):
    """One row per exon and per intron of every refGene transcript.

    The interval is an int4range with a GiST index instead of the start/end pair the legacy
    schema kept: everything that reads this table - the capture regions of 14_capture_regions,
    the variant mapping of 12_gnomad - is asking what overlaps what, and '&&' answers that.
    """

    __tablename__ = "ucsc_refgene_regions"
    __table_args__ = (
        # the region is part of the key: refGene places some transcripts more than once on the
        # same chromosome, and those placements are different rows, not a duplicate
        UniqueConstraint("assembly", "transcript_id", "chrom", "feature", "region",
                         name="uq_ucsc_refgene_regions_feature"),
        Index("ix_ucsc_refgene_regions_region", "region", postgresql_using="gist"),
    )

    id: Optional[int] = Field(default=None, primary_key=True)
    gene_id: int = Field(foreign_key="genes.id", index=True)
    assembly: str = Field(index=True)                 # hg19
    transcript_id: str = Field(index=True)            # RefSeq accession, unversioned: NM_130786
    chrom: str = Field(index=True)
    region: Any = Field(sa_type=INT4RANGE)            # zero based, half open, as in the bed files
    strand: str
    feature: str = Field(index=True)                  # exon_1, intron_1, ...
