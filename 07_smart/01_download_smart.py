#! /usr/bin/python3
"""Submit protein sequences to SMART, keep the plain text results, and fill the smart tables.

usage:  01_download_smart.py [--genes ABCA4,MYO7A,USH2A] [--output-directory DIR]
                             [--fasta FILE] [--include-pfam] [--include-signalp]
                             [--include-repeats] [--include-schnipsel] [--no-store]

A python translation of SMART_batch.pl (Ivica Letunic, https://smart.embl.de), with the FASTA
file replaced by a selection out of 'uniprot_basic_infos': the sequence SMART is asked about is
the one already stored for the gene, so the domain boundaries it returns are in the residue
numbering the rest of the database uses. --fasta restores the original behaviour, for a sequence
that is not in the database; those results are written to disk but not stored, since there is no
gene to attach them to.

Selection: the reviewed uniprot entries of the named genes. A symbol that does not resolve stops
the run before anything is submitted, and a gene with no uniprot entry is an error too - SMART
cannot be asked about a protein we do not have.

The submission protocol, unchanged from the perl: one sequence at a time, always waiting for the
result before the next is sent, polling the queue every ten seconds. A results file that already
exists is not fetched again, so an interrupted run resumes where it stopped.
"""

import argparse
import hashlib
import math
import os
import re
import sys
import time
from datetime import datetime, timezone
from typing import Any, Dict, Iterator, List, Optional, Tuple

import requests
from sqlmodel import Session, delete, select

from integrator_utils.python.db import db_engine, gene_ids, require_genes, upsert, upsert_many
from integrator_utils.python.models_core import UniprotBasicInfo

from smart_models import SmartAnalysis, SmartDomain

SUBMIT_URL = "http://ismart.embl.de/smart/show_motifs.pl"
JOB_STATUS_URL = "http://ismart.embl.de/results.cgi"
USER_AGENT = "SMARTbatch1.0"
OUTPUT_DIRECTORY = "/storage/databases/smart"
DEFAULT_GENES = "ABCA4,MYO7A,USH2A"

RESULT_HEADER = "-- SMART RESULT"
JOB_ID_PATTERN = re.compile(r"results\.cgi\?id=(\d+)")
FIELD_PATTERN = re.compile(r"^([A-Z_]+)=(.*)$")

TIMEOUT = 120          # seconds a single request may take
FIRST_POLL = 5         # the perl waits this long before asking about a fresh job
POLL_INTERVAL = 10     # and this long between polls
BE_NICE = 5            # ... and this long between submissions

# the option name the form expects, per --include-<x> flag
FORM_FLAGS = {"pfam": "DO_PFAM", "signalp": "INCLUDE_SIGNALP", "repeats": "DO_PROSPERO",
              "schnipsel": "INCLUDE_BLAST"}


#########################################
# the sequences to submit
def fasta(path: str) -> Iterator[Tuple[str, str]]:
    """(display id, sequence) per entry - the Bio::SeqIO of the perl, in ten lines."""
    name: Optional[str] = None
    residues: List[str] = []
    with open(path) as inf:
        for line in inf:
            if not line.startswith(">"):
                residues.append(line.strip())
                continue
            if name is not None:
                yield name, "".join(residues)
            name = line[1:].split()[0]
            residues = []
    if name is not None:
        yield name, "".join(residues)


def sequences_from_db(engine, symbols: List[str]) -> List[Tuple[str, str, int]]:
    """(uniprot accession, sequence, gene id) for every reviewed entry of the named genes."""
    proteins: List[Tuple[str, str, int]] = []
    missing: List[str] = []
    with Session(engine) as session:
        resolved = gene_ids(session, symbols)   # a symbol that does not resolve exits here
        for symbol, gene_id in resolved.items():
            statement = select(UniprotBasicInfo).where(UniprotBasicInfo.gene_id == gene_id)
            entries = [entry for entry in session.exec(statement) if entry.sequence]
            if not entries:
                missing.append(symbol)
            for entry in sorted(entries, key=lambda entry: entry.uniprot_id):
                proteins.append((entry.uniprot_id, entry.sequence, gene_id))
    if missing:
        offenders = ", ".join(missing)
        errmsg = f"no uniprot sequence in 'uniprot_basic_infos' for: {offenders} - run 03_uniprot first"
        sys.exit(errmsg)
    return proteins


#########################################
# the submission, as in SMART_batch.pl
def is_result(text: str) -> bool:
    """SMART sometimes writes a blank line before the header, which is why two lines are looked at."""
    return any(line.startswith(RESULT_HEADER) for line in text.split("\n", 2)[:2])


def submit(session: requests.Session, sequence: str, options: List[str]) -> requests.Response:
    fields = {"SEQUENCE": sequence, "TEXTONLY": "1"}
    for option in options:
        fields[FORM_FLAGS[option]] = "1"
    # the form is multipart, and requests writes multipart only for 'files'; (None, value) is a
    # plain field rather than an attachment
    parts = {name: (None, value) for name, value in fields.items()}
    return session.post(SUBMIT_URL, files=parts, allow_redirects=False, timeout=TIMEOUT)


def wait_for_job(session: requests.Session, job_id: str) -> str:
    """Poll the queue until the job comes back with a result."""
    queued = f"Job entered the queue with ID {job_id}. Waiting for results."
    print(queued)
    time.sleep(FIRST_POLL)
    while True:
        response = session.get(JOB_STATUS_URL, params={"id": job_id}, timeout=TIMEOUT)
        if not response.ok:
            errmsg = f"SMART returned a web server error: {response.status_code} {response.reason}"
            sys.exit(errmsg)
        if is_result(response.text):
            return response.text
        time.sleep(POLL_INTERVAL)


def fetch(session: requests.Session, sequence: str, options: List[str], error_path: str) -> str:
    """The result text for one sequence, whether it comes from the queue or precomputed."""
    response = submit(session, sequence, options)
    if response.is_redirect or response.is_permanent_redirect:
        location = response.headers.get("Location", "")
        match = JOB_ID_PATTERN.search(location)
        if not match:
            errmsg = f"could not get the job ID from the redirect header ({location}). Aborting further submissions."
            sys.exit(errmsg)
        return wait_for_job(session, match.group(1))

    if not response.ok:
        errmsg = f"SMART returned a web server error: {response.status_code} {response.reason}"
        sys.exit(errmsg)
    if is_result(response.text):
        return response.text   # precomputed, SMART had seen this sequence before

    with open(error_path, "w") as outf:
        outf.write(response.text)
    errmsg = (f"SMART returned an error page, which was saved into '{error_path}'.\n"
              f"Please check the file for details. Aborting further submissions.")
    sys.exit(errmsg)


#########################################
# the text format: a header block, then one block per feature, then '-- FINISHED --'
def parse_result(text: str) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    header: Dict[str, str] = {}
    features: List[Dict[str, str]] = []
    for line in text.split("\n"):
        match = FIELD_PATTERN.match(line.strip())
        if not match:
            continue
        key, value = match.group(1), match.group(2).strip()
        if key == "DOMAIN":
            features.append({"DOMAIN": value})
        elif features:
            features[-1][key] = value
        else:
            header[key] = value

    domains: List[Dict[str, Any]] = []
    for feature in features:
        status = feature.get("STATUS", "")
        domains.append({"name": feature["DOMAIN"], "start": int(feature["START"]), "end": int(feature["END"]),
                        "evalue": as_float(feature.get("EVALUE")), "feature_type": feature.get("TYPE", "unknown"),
                        "status": status, "visible": status.split("|")[0] == "visible"})
    return header, domains


def as_float(value: Optional[str]) -> Optional[float]:
    """SMART writes NaN for the features it has no e-value for - a transmembrane region, say."""
    if not value:
        return None
    try:
        number = float(value)
    except ValueError:
        return None
    return None if math.isnan(number) else number


#########################################
def store(engine, uniprot_id: str, gene_id: int, sequence: str, options: List[str],
          header: Dict[str, Any], domains: List[Dict[str, Any]]) -> None:
    crc = header.get("CRC_PASSED")
    analysis = {"uniprot_id": uniprot_id, "gene_id": gene_id,
                "sequence_md5": hashlib.md5(sequence.encode()).hexdigest(), "aa_length": len(sequence),
                "crc_passed": None if crc is None else crc == "1", "feature_count": len(domains),
                "options": options or None, "retrieved_at": datetime.now(timezone.utc)}
    with Session(engine) as session:
        analysis_id = upsert(session, SmartAnalysis, analysis, conflict_on=["uniprot_id"])
        # the stored set is replaced rather than merged: re-run with other options, or against a
        # revised sequence, and the features that are gone have to go from the table too
        session.exec(delete(SmartDomain).where(SmartDomain.analysis_id == analysis_id))
        rows = [dict(domain, analysis_id=analysis_id) for domain in domains]
        upsert_many(session, SmartDomain, rows, conflict_on=["analysis_id", "name", "start", "end"])
        session.commit()


#########################################
def run(engine, proteins: List[Tuple[str, str, Optional[int]]], directory: str, options: List[str],
        no_store: bool) -> None:
    http = requests.Session()
    http.headers.update({"User-Agent": USER_AGENT})
    for uniprot_id, sequence, gene_id in proteins:
        output_path = os.path.join(directory, f"{uniprot_id}_SMART_results.txt")
        error_path = os.path.join(directory, f"{uniprot_id}_SMART_error.html")
        text = cached(output_path)
        if text is None:
            submitting = f"Submitting sequence {uniprot_id} ({len(sequence)} aa)..."
            print(submitting)
            text = fetch(http, sequence, options, error_path)
            with open(output_path, "w") as outf:
                outf.write(text)
            saved = f"Results saved to '{output_path}'"
            print(saved)
            time.sleep(BE_NICE)   # be nice to other users

        header, domains = parse_result(text)
        if no_store or gene_id is None:
            report = f"{uniprot_id}: {len(domains)} feature(s), not stored"
            print(report)
            continue
        store(engine, uniprot_id, gene_id, sequence, options, header, domains)
        visible = len([domain for domain in domains if domain["visible"]])
        report = f"{uniprot_id}: {len(domains)} feature(s) stored, {visible} of them visible"
        print(report)


def cached(output_path: str) -> Optional[str]:
    """The perl's skip rule: an existing results file is kept, an empty one is thrown away."""
    if not os.path.exists(output_path):
        return None
    if os.path.getsize(output_path) == 0:
        removing = f"Removing empty results file {output_path}."
        print(removing)
        os.unlink(output_path)
        return None
    skipping = f"Using the results file {output_path}, which already exists."
    print(skipping)
    with open(output_path) as inf:
        return inf.read()


#########################################
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="submit sequences to SMART and store the domains found")
    parser.add_argument("--genes", default=DEFAULT_GENES, help="comma separated gene symbols")
    parser.add_argument("--fasta", help="submit the sequences in this file instead (not stored)")
    parser.add_argument("--output-directory", default=OUTPUT_DIRECTORY, help="where the result files go")
    parser.add_argument("--include-pfam", action="store_true", help="include Pfam domains in the search")
    parser.add_argument("--include-signalp", action="store_true", help="include signal peptide predictions")
    parser.add_argument("--include-repeats", action="store_true", help="include internal repeat predictions")
    parser.add_argument("--include-schnipsel", action="store_true",
                        help="include outlier homologues and homologues of known structures")
    parser.add_argument("--no-store", action="store_true", help="download only, leave the database alone")
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_directory, exist_ok=True)
    options = [name for name in FORM_FLAGS if getattr(args, f"include_{name}")]

    engine = db_engine()
    require_genes(engine)
    if args.fasta:
        if not os.path.exists(args.fasta):
            errmsg = f"{args.fasta} does not exist"
            sys.exit(errmsg)
        proteins = [(name, sequence, None) for name, sequence in fasta(args.fasta)]
    else:
        proteins = sequences_from_db(engine, [symbol.strip() for symbol in args.genes.split(",") if symbol.strip()])

    announce = f"SMART batch analysis\n======================\n{len(proteins)} sequence(s) to do"
    print(announce)
    run(engine, proteins, args.output_directory, options, args.no_store)


#########################################
if __name__ == '__main__':
    main()
