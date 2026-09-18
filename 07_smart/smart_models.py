#! /usr/bin/python3
"""Tables filled from the SMART domain annotation service (https://smart.embl.de).

SMART annotates a protein sequence, not a gene, so the analysis is keyed by the uniprot accession
whose sequence was submitted; 'genes' is reached through it. The submitted sequence is recorded
as an md5, which is the only way to notice that uniprot has revised the sequence since the last
run and that the stored domains describe a protein that no longer exists.
"""

from datetime import datetime
from typing import List, Optional

from sqlalchemy import Text, UniqueConstraint
from sqlalchemy.dialects.postgresql import ARRAY
from sqlmodel import Field, SQLModel

# the TYPE= line of the text output, and with it the source of the annotation
FEATURE_TYPES = ("SMART", "Pfam", "unknown")


#########################################
class SmartAnalysis(SQLModel, table=True):
    """One row per protein submitted - the job, rather than what it found.

    Kept separate from the domains so that 'submitted, nothing found' is a fact the database can
    state: a protein with no row here has never been through SMART, one with a row and no domains
    has been through it and came back empty.
    """

    __tablename__ = "smart_analyses"

    id: Optional[int] = Field(default=None, primary_key=True)
    uniprot_id: str = Field(index=True, unique=True)
    gene_id: int = Field(foreign_key="genes.id", index=True)
    sequence_md5: str                                          # of the sequence actually submitted
    aa_length: int
    crc_passed: Optional[bool] = None                          # CRC_PASSED= of the text output
    feature_count: int = 0                                     # NUMBER_OF_FEATURES_FOUND=
    # the analysis options the submission carried: pfam, signalp, repeats, schnipsel
    options: Optional[List[str]] = Field(default=None, sa_type=ARRAY(Text))
    retrieved_at: Optional[datetime] = None


class SmartDomain(SQLModel, table=True):
    """One row per feature SMART reports, in the residue numbering of the submitted sequence.

    Hidden features are stored too. SMART hides everything below its threshold, but the e-value
    is right there in the row, so the caller can pick a cutoff of their own; throwing them away
    here would mean re-running the whole set to change one's mind.
    """

    __tablename__ = "smart_domains"
    __table_args__ = (
        UniqueConstraint("analysis_id", "name", "start", "end", name="uq_smart_domains_feature"),
    )

    id: Optional[int] = Field(default=None, primary_key=True)
    analysis_id: int = Field(foreign_key="smart_analyses.id", index=True)
    name: str = Field(index=True)                              # DOMAIN=, "TUDOR", "transmembrane_region"
    start: int                                                 # one based, inclusive, as SMART writes it
    end: int
    evalue: Optional[float] = None                             # None where SMART writes NaN
    feature_type: str = Field(index=True)                      # one of FEATURE_TYPES
    status: str                                                # "visible|OK", "hidden|threshold", ...
    visible: bool = Field(index=True)                          # the first field of status
