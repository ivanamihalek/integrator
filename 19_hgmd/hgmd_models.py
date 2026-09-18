#! /usr/bin/python3
"""SQLModel description of the tables filled from HGMD per-gene mutation pages."""

import re
from typing import Any, Dict, List, Optional, Type

from sqlalchemy.dialects.postgresql import JSONB
from sqlmodel import Field, SQLModel

from integrator_utils.python.models_core import Publication  # noqa: F401  - the target of the links below


#########################################
# mutation tables - the columns common to all of them
class MutationBase(SQLModel):
    id: Optional[int] = Field(default=None, primary_key=True)
    gene_id: int = Field(foreign_key="genes.id", index=True)
    accession: str = Field(index=True, unique=True)
    variant_class: Optional[str] = None
    phenotype: Optional[str] = None
    # sa_type (rather than sa_column) so that each inheriting table gets its own Column instance
    xrefs: Optional[Dict[str, Any]] = Field(default=None, sa_type=JSONB)


class MissenseNonsense(MutationBase, table=True):
    __tablename__ = "hgmd_missense_nonsense"

    codon_change: Optional[str] = None
    aa_change: Optional[str] = None
    hgvs_nucleotide: Optional[str] = None
    hgvs_protein: Optional[str] = None
    vcf: Optional[str] = None
    # four genes (COL11A1, COL2A1, RBP4, TIMP3) carry a second, historical numbering
    legacy_change: Optional[str] = None


class Splicing(MutationBase, table=True):
    __tablename__ = "hgmd_splicing"

    splicing_mutation: Optional[str] = None
    hgvs_nucleotide: Optional[str] = None
    vcf: Optional[str] = None


class SmallIndel(MutationBase, table=True):
    __tablename__ = "hgmd_small_indels"

    variant_type: str = Field(index=True)  # "deletion", "insertion" or "indel"
    description: Optional[str] = None
    hgvs_nucleotide: Optional[str] = None
    hgvs_protein: Optional[str] = None
    vcf: Optional[str] = None


class GrossIndel(MutationBase, table=True):
    __tablename__ = "hgmd_gross_indels"

    # "deletion" or "insertion" from the section title; the gross insertion tables carry an
    # 'Insertion/ duplication' column of their own, and that one wins where it is present
    variant_type: str = Field(index=True)
    dna_level: Optional[str] = None
    description: Optional[str] = None
    hgvs_nucleotide: Optional[str] = None
    hgvs_protein: Optional[str] = None
    vcf: Optional[str] = None


class Regulatory(MutationBase, table=True):
    __tablename__ = "hgmd_regulatory"

    regulatory_sequence: Optional[str] = None
    hgvs_nucleotide: Optional[str] = None
    vcf: Optional[str] = None


class ComplexRearrangement(MutationBase, table=True):
    __tablename__ = "hgmd_complex_rearrangements"

    description: Optional[str] = None


class RepeatVariation(MutationBase, table=True):
    __tablename__ = "hgmd_repeat_variations"

    amplified_sequence: Optional[str] = None
    location: Optional[str] = None
    normal_range: Optional[str] = None
    pathological_range: Optional[str] = None


#########################################
# through tables, mapping mutation entries to publications
class MutationPublicationBase(SQLModel):
    publication_id: int = Field(foreign_key="publications.id", primary_key=True)
    is_primary: bool = False
    annotation: Optional[str] = None


class MissenseNonsensePublication(MutationPublicationBase, table=True):
    __tablename__ = "hgmd_missense_nonsense_publication"

    entry_id: int = Field(foreign_key="hgmd_missense_nonsense.id", primary_key=True)


class SplicingPublication(MutationPublicationBase, table=True):
    __tablename__ = "hgmd_splicing_publication"

    entry_id: int = Field(foreign_key="hgmd_splicing.id", primary_key=True)


class SmallIndelPublication(MutationPublicationBase, table=True):
    __tablename__ = "hgmd_small_indels_publication"

    entry_id: int = Field(foreign_key="hgmd_small_indels.id", primary_key=True)


class GrossIndelPublication(MutationPublicationBase, table=True):
    __tablename__ = "hgmd_gross_indels_publication"

    entry_id: int = Field(foreign_key="hgmd_gross_indels.id", primary_key=True)


class RegulatoryPublication(MutationPublicationBase, table=True):
    __tablename__ = "hgmd_regulatory_publication"

    entry_id: int = Field(foreign_key="hgmd_regulatory.id", primary_key=True)


class ComplexRearrangementPublication(MutationPublicationBase, table=True):
    __tablename__ = "hgmd_complex_rearrangements_publication"

    entry_id: int = Field(foreign_key="hgmd_complex_rearrangements.id", primary_key=True)


class RepeatVariationPublication(MutationPublicationBase, table=True):
    __tablename__ = "hgmd_repeat_variations_publication"

    entry_id: int = Field(foreign_key="hgmd_repeat_variations.id", primary_key=True)


#########################################
# the HGMD page layout: which section goes into which table, and how its
# column headers (normalized to lowercase, whitespace collapsed) map to fields
HEADER_TO_FIELD: Dict[str, str] = {
    "hgmd accession": "accession",
    "hgmd codon change": "codon_change",
    "hgmd amino acid change": "aa_change",
    "hgmd splicing mutation": "splicing_mutation",
    "hgmd deletion": "description",
    "hgmd insertion": "description",
    "hgmd deletion/insertion": "description",
    "insertion/ duplication": "variant_type",
    "regulatory sequence": "regulatory_sequence",
    "amplified sequence": "amplified_sequence",
    "location": "location",
    "normal range": "normal_range",
    "pathological range": "pathological_range",
    "description": "description",
    "dna level": "dna_level",
    "hgvs (nucleotide)": "hgvs_nucleotide",
    "hgvs (protein)": "hgvs_protein",
    "vcf chrom pos ref/alt": "vcf",
    "variant class": "variant_class",
    "reported phenotype": "phenotype",
    "reference": "reference",
    "extra information": "xrefs",
}

# a handful of genes number their variants twice, and say so in the column header:
# 'HGMD deletion (^legacy ATG=-528)' is the ordinary deletion column, 'Legacy change legacy
# ATG=-528' is the historical numbering of the variant
LEGACY_SUFFIX_PATTERN = re.compile(r"\s*\(\^legacy[^)]*\)\s*$")
LEGACY_CHANGE_PATTERN = re.compile(r"^legacy change\b")


def field_for(header: str) -> Optional[str]:
    """Model field for a normalized (lowercased, whitespace collapsed) column header."""
    header = LEGACY_SUFFIX_PATTERN.sub("", header)
    if LEGACY_CHANGE_PATTERN.match(header):
        return "legacy_change"
    return HEADER_TO_FIELD.get(header)


class Section:
    """One HGMD page section: its h3 title, the model it feeds, and the variant type it implies."""

    def __init__(self, title: str, model: Type[MutationBase], link_model: Type[MutationPublicationBase],
                 variant_type: Optional[str] = None) -> None:
        self.title: str = title
        self.model: Type[MutationBase] = model
        self.link_model: Type[MutationPublicationBase] = link_model
        self.variant_type: Optional[str] = variant_type


SECTIONS: List[Section] = [
    Section("Missense/nonsense", MissenseNonsense, MissenseNonsensePublication),
    Section("Splicing", Splicing, SplicingPublication),
    Section("Regulatory", Regulatory, RegulatoryPublication),
    Section("Small deletions", SmallIndel, SmallIndelPublication, "deletion"),
    Section("Small insertions", SmallIndel, SmallIndelPublication, "insertion"),
    Section("Small indels", SmallIndel, SmallIndelPublication, "indel"),
    Section("Gross deletions", GrossIndel, GrossIndelPublication, "deletion"),
    Section("Gross insertions", GrossIndel, GrossIndelPublication, "insertion"),
    Section("Complex rearrangements", ComplexRearrangement, ComplexRearrangementPublication),
    Section("Repeat variations", RepeatVariation, RepeatVariationPublication),
]
