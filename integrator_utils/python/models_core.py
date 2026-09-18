#! /usr/bin/python3
"""The tables shared across directories: the gene identity hub, and the publications every
annotation source cites.

'genes' and 'gene_aliases' are written by 01_hgnc alone. Everything else resolves symbols
against them through integrator_utils.python.db.gene_id().

Symbol columns are citext: MySQL compared them case-insensitively, postgres does not, and every
source spells gene names in its own case.
"""

from typing import List, Optional

from sqlalchemy import Text
from sqlalchemy.dialects.postgresql import ARRAY, CITEXT
from sqlmodel import Field, SQLModel

ALIAS_TYPES = ("alias", "previous")


#########################################
class Gene(SQLModel, table=True):
    __tablename__ = "genes"

    id: Optional[int] = Field(default=None, primary_key=True)
    hgnc_id: str = Field(index=True, unique=True)                       # "HGNC:5"
    symbol: str = Field(sa_type=CITEXT, index=True, unique=True)        # HGNC approved symbol
    name: Optional[str] = None                                          # HGNC approved name
    locus_group: Optional[str] = Field(default=None, index=True)        # "protein-coding gene", ...
    locus_type: Optional[str] = None
    location: Optional[str] = None                                      # cytogenetic band, "19q13.43"
    # not unique: HGNC maps a handful of distinct symbols onto one ncbi/ensembl id
    ncbi_gene_id: Optional[int] = Field(default=None, index=True)
    ensembl_gene_id: Optional[str] = Field(default=None, index=True)
    ucsc_id: Optional[str] = Field(default=None, index=True)
    # pipe-joined in the HGNC download, kept as arrays rather than strings to be split at every use
    uniprot_ids: Optional[List[str]] = Field(default=None, sa_type=ARRAY(Text))
    refseq_ids: Optional[List[str]] = Field(default=None, sa_type=ARRAY(Text))
    ccds_ids: Optional[List[str]] = Field(default=None, sa_type=ARRAY(Text))
    mgd_ids: Optional[List[str]] = Field(default=None, sa_type=ARRAY(Text))


class GeneAlias(SQLModel, table=True):
    """Aliases and previous symbols, one row each - the legacy comma-joined blobs unpacked.

    An alias is not unique: 1489 of them are shared by several genes, and 700 are somebody else's
    approved symbol. Resolution therefore tries 'genes' first and refuses an ambiguous alias.
    """

    __tablename__ = "gene_aliases"

    gene_id: int = Field(foreign_key="genes.id", primary_key=True, index=True)
    alias: str = Field(sa_type=CITEXT, primary_key=True, index=True)
    alias_type: str = Field(primary_key=True)  # one of ALIAS_TYPES


#########################################
class UniprotBasicInfo(SQLModel, table=True):
    """One row per reviewed (SwissProt) human protein, from the uniprot flat file.

    The canonical sequence lives here rather than in a separate table: the flat file carries it,
    so there is no reason to go back to a blast database for it. 'uniprot_seqs' is left to the
    transcript-to-genome coordinate mapping, which needs UCSC and is filled later.
    """

    __tablename__ = "uniprot_basic_infos"

    id: Optional[int] = Field(default=None, primary_key=True)
    uniprot_id: str = Field(index=True, unique=True)
    gene_id: int = Field(foreign_key="genes.id", index=True)
    full_name: Optional[str] = None
    ec_number: Optional[str] = Field(default=None, index=True)
    canonical_aa_length: Optional[int] = None
    sequence: Optional[str] = None
    subcellular_location: Optional[str] = None
    function: Optional[str] = None
    ensembl_gene_ids: Optional[List[str]] = Field(default=None, sa_type=ARRAY(Text))
    cofactors: Optional[List[str]] = Field(default=None, sa_type=ARRAY(Text))
    tissues: Optional[List[str]] = Field(default=None, sa_type=ARRAY(Text))
    # the accessions this entry has been known by before - 'old_ids' in the legacy schema
    secondary_ids: Optional[List[str]] = Field(default=None, sa_type=ARRAY(Text))


#########################################
class Publication(SQLModel, table=True):
    __tablename__ = "publications"

    id: Optional[int] = Field(default=None, primary_key=True)
    reference: str                                                      # "Harville (2010) J Med Genet 47, 262"
    pubmed_id: Optional[int] = Field(default=None, index=True, unique=True)


#########################################
class OmimGene(SQLModel, table=True):
    """The OMIM mim2gene mapping: which MIM number is which gene, and what kind of entry it is.

    The phenotype side of OMIM - genemap2.txt, and with it 'omim_genemaps' - needs a licensed
    download that is not on this machine; 20_omim/10_mark_inherited_metabolic_disorders.py stays
    on the old stack until it is.
    """

    __tablename__ = "omim_genes"

    id: Optional[int] = Field(default=None, primary_key=True)
    mim_number: int = Field(index=True, unique=True)
    entry_type: str = Field(index=True)          # gene, phenotype, gene/phenotype, moved/removed
    gene_id: Optional[int] = Field(default=None, foreign_key="genes.id", index=True)
    ncbi_gene_id: Optional[int] = Field(default=None, index=True)
    ensembl_gene_id: Optional[str] = Field(default=None, index=True)


class HgncWithdrawn(SQLModel, table=True):
    """HGNC ids that no longer name an approved gene, from the HGNC withdrawn list.

    Loaders need the distinction between 'this id was withdrawn' - a fact about HGNC, and a
    legitimate reason for a source record to fall outside the selection - and 'this id is
    unknown', which means the pipeline is wrong somewhere and has to stop.
    """

    __tablename__ = "hgnc_withdrawn"

    hgnc_id: str = Field(primary_key=True)
    status: str                                  # "Entry Withdrawn" or "Merged/Split"
    withdrawn_symbol: Optional[str] = Field(default=None, sa_type=CITEXT, index=True)
    merged_into: Optional[List[str]] = Field(default=None, sa_type=ARRAY(Text))
