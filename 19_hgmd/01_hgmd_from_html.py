#! /usr/bin/python3
"""Store the mutation tables from HGMD gene pages (saved as html) in the integrator database.

usage:  01_hgmd_from_html.py <hgmd_gene_page.html | directory of them>

The gene must already be in 'genes' - 01_hgnc fills that table, this script only resolves
against it. A directory argument is resolved in full before anything is stored, so a symbol
HGMD spells differently from HGNC stops the run before the first insert.
"""

import os
import re
import sys
from typing import Any, Dict, List, Optional, Tuple, Type

from bs4 import BeautifulSoup, Tag
from sqlalchemy.dialects.postgresql import insert
from sqlmodel import Session, select

from integrator_utils.python.db import db_engine, gene_ids, require_genes
from integrator_utils.python.models_core import Publication

from hgmd_models import SECTIONS, MutationBase, MutationPublicationBase, Section, field_for

NOT_AVAILABLE = "not yet available"
PUBMED_ID_PATTERN = re.compile(r"pubmed\.ncbi\.nlm\.nih\.gov/(\d+)")
GENE_TITLE_PATTERN = re.compile(r"^All\s+\d+\s+mutations", re.IGNORECASE)
ANNOTATION_PATTERN = re.compile(r"\[([^\[\]]+)\]\s*$")


#########################################
def normalize(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def cell_text(cell: Tag) -> Optional[str]:
    text = normalize(cell.get_text(" ", strip=True))
    if not text or text.lower() == NOT_AVAILABLE:
        return None
    return text


def accession_from_cell(cell: Tag) -> Optional[str]:
    hidden_input = cell.find("input", attrs={"name": "acc"})
    if hidden_input and hidden_input.get("value"):
        return normalize(hidden_input["value"])
    return cell_text(cell)


def xrefs_from_cell(cell: Tag) -> Optional[Dict[str, Dict[str, Optional[str]]]]:
    """Keep the linked entries only (dbSNP, gnomAD, ...); the hg19/hg38/COM/CpG spans carry no href."""
    xrefs: Dict[str, Dict[str, Optional[str]]] = {}
    for anchor in cell.find_all("a", href=True):
        label = normalize(anchor.get_text(" ", strip=True))
        if not label:
            continue
        inner_span = anchor.find("span", title=True)
        title = inner_span["title"] if inner_span else anchor.get("title")
        xrefs[label] = {"url": anchor["href"], "title": normalize(title) if title else None}
    return xrefs or None


def references_from_cell(cell: Tag) -> List[Tuple[str, Optional[int], bool, Optional[str]]]:
    """(citation, pubmed id, is_primary, annotation) for each linked publication, in document order."""
    references: List[Tuple[str, Optional[int], bool, Optional[str]]] = []
    for position, anchor in enumerate(cell.find_all("a", href=True)):
        citation = normalize(anchor.get_text(" ", strip=True))
        if not citation:
            continue
        pubmed_match = PUBMED_ID_PATTERN.search(anchor["href"])
        pubmed_id = int(pubmed_match.group(1)) if pubmed_match else None
        # the additional reports come wrapped in a span, annotated as e.g. "[Additional report]"
        wrapper = anchor.find_parent("span")
        annotation_match = ANNOTATION_PATTERN.search(normalize(wrapper.get_text(" ", strip=True))) if wrapper else None
        annotation = annotation_match.group(1) if annotation_match else None
        references.append((citation, pubmed_id, position == 0, annotation))
    return references


#########################################
def gene_name(soup: BeautifulSoup) -> Optional[str]:
    for header in soup.find_all("h3"):
        title = normalize(header.get_text(" ", strip=True))
        if not GENE_TITLE_PATTERN.match(title):
            continue
        anchor = header.find("a")
        if anchor:
            return normalize(anchor.get_text(" ", strip=True))
        return title.split()[-1]
    return None  # not an HGMD mutation page - see main()


def section_table(soup: BeautifulSoup, section_title: str) -> Optional[Tag]:
    for header in soup.find_all("h3"):
        title = normalize(header.get_text(" ", strip=True))
        if title.lower().startswith(section_title.lower()):
            return header.find_next("table")
    return None


def field_names(header_row: Tag, section: Section) -> List[Optional[str]]:
    """Column index to model field name; None for the columns we do not store."""
    fields: List[Optional[str]] = []
    for header_cell in header_row.find_all("th"):
        header = normalize(header_cell.get_text(" ", strip=True)).lower()
        field = field_for(header)
        if field is None:
            warning = f"\t warning: unrecognized column '{header}' in the '{section.title}' table - skipping it"
            print(warning)
        fields.append(field)
    return fields


def parse_section(table: Tag, section: Section) -> List[Dict[str, Any]]:
    """One dict per mutation row: model fields, plus the 'references' list to be stored separately."""
    rows = table.find_all("tr")
    header_rows = [row for row in rows if row.find("th")]
    if not header_rows:
        warning = f"\t warning: no header row in the '{section.title}' table - skipping the section"
        print(warning)
        return []
    header_row = header_rows[0]
    fields = field_names(header_row, section)

    entries: List[Dict[str, Any]] = []
    for row in rows[rows.index(header_row) + 1:]:
        cells = row.find_all("td")
        if len(cells) != len(fields):  # the column-toggle form row, or some other page furniture
            continue
        entry: Dict[str, Any] = {"references": []}
        if section.variant_type:
            entry["variant_type"] = section.variant_type
        for field, cell in zip(fields, cells):
            if field is None:
                continue
            elif field == "accession":
                entry["accession"] = accession_from_cell(cell)
            elif field == "xrefs":
                entry["xrefs"] = xrefs_from_cell(cell)
            elif field == "reference":
                entry["references"] = references_from_cell(cell)
            elif field == "variant_type":
                # the 'Insertion/ duplication' column of the gross insertion tables, which is
                # more specific than the "insertion" implied by the section title
                entry["variant_type"] = (cell_text(cell) or section.variant_type).lower()
            else:
                entry[field] = cell_text(cell)
        if entry.get("accession"):
            entries.append(entry)
    return entries


def page_soup(html_path: str) -> BeautifulSoup:
    """Some pages were saved from the browser's source viewer rather than from the page itself:
    the real html sits escaped in the viewer's line-content cells and is put back together here."""
    with open(html_path, encoding="utf-8", errors="replace") as inf:
        soup = BeautifulSoup(inf.read(), "lxml")
    if soup.find("h3"):
        return soup
    source_lines = soup.find_all("td", class_="line-content")
    if not source_lines:
        return soup
    return BeautifulSoup("\n".join(line.get_text() for line in source_lines), "lxml")


def parse_html(html_path: str) -> Tuple[Optional[str], List[Tuple[Section, List[Dict[str, Any]]]]]:
    soup = page_soup(html_path)
    name = gene_name(soup)
    if name is None:
        return None, []
    parsed: List[Tuple[Section, List[Dict[str, Any]]]] = []
    for section in SECTIONS:
        table = section_table(soup, section.title)
        if table is None:
            continue
        parsed.append((section, parse_section(table, section)))
    return name, parsed


#########################################
def store_publications(session: Session, parsed: List[Tuple[Section, List[Dict[str, Any]]]]) -> Dict[str, int]:
    """Store all publications found in the page; return a (citation or pmid) keyed map of their ids."""
    with_pubmed_id: Dict[int, str] = {}
    without_pubmed_id: Dict[str, None] = {}
    for _, entries in parsed:
        for entry in entries:
            for citation, pubmed_id, _, _ in entry["references"]:
                if pubmed_id is None:
                    without_pubmed_id[citation] = None
                elif pubmed_id not in with_pubmed_id:  # the first spelling of the citation wins
                    with_pubmed_id[pubmed_id] = citation

    publication_id: Dict[str, int] = {}
    if with_pubmed_id:
        rows = [{"pubmed_id": pubmed_id, "reference": citation} for pubmed_id, citation in with_pubmed_id.items()]
        statement = insert(Publication).values(rows).on_conflict_do_nothing(index_elements=["pubmed_id"])
        session.exec(statement)
        session.commit()
        stored = session.exec(select(Publication).where(Publication.pubmed_id.in_(with_pubmed_id.keys()))).all()
        for publication in stored:
            publication_id[str(publication.pubmed_id)] = publication.id

    for citation in without_pubmed_id:
        publication = session.exec(select(Publication).where(Publication.reference == citation)).first()
        if publication is None:
            publication = Publication(reference=citation)
            session.add(publication)
            session.commit()
            session.refresh(publication)
        publication_id[citation] = publication.id

    return publication_id


def store_entry(session: Session, model: Type[MutationBase], entry: Dict[str, Any]) -> int:
    values = {field: value for field, value in entry.items() if field != "references"}
    updatable = {field: value for field, value in values.items() if field != "accession"}
    statement = insert(model).values(values).on_conflict_do_update(index_elements=["accession"], set_=updatable)
    return session.exec(statement.returning(model.__table__.c.id)).one()[0]


def store_links(session: Session, link_model: Type[MutationPublicationBase], entry_id: int,
                references: List[Tuple[str, Optional[int], bool, Optional[str]]],
                publication_id: Dict[str, int]) -> None:
    # keyed by publication id: a single statement must not touch the same row twice
    rows: Dict[int, Dict[str, Any]] = {}
    for citation, pubmed_id, is_primary, annotation in references:
        key = str(pubmed_id) if pubmed_id is not None else citation
        if key not in publication_id or publication_id[key] in rows:
            continue
        link = {"entry_id": entry_id, "publication_id": publication_id[key]}
        link.update({"is_primary": is_primary, "annotation": annotation})
        rows[publication_id[key]] = link
    if not rows:
        return
    index_elements = ["entry_id", "publication_id"]
    statement = insert(link_model).values(list(rows.values()))
    updatable = {"is_primary": statement.excluded.is_primary, "annotation": statement.excluded.annotation}
    statement = statement.on_conflict_do_update(index_elements=index_elements, set_=updatable)
    session.exec(statement)


def store(session: Session, gene_id: int, parsed: List[Tuple[Section, List[Dict[str, Any]]]]) -> None:
    publication_id = store_publications(session, parsed)
    for section, entries in parsed:
        for entry in entries:
            entry["gene_id"] = gene_id
            entry_id = store_entry(session, section.model, entry)
            store_links(session, section.link_model, entry_id, entry["references"], publication_id)
        session.commit()
        report = f"\t {section.title}: {len(entries)} entries stored in {section.model.__tablename__}"
        print(report)
    summary = f"\t publications: {len(publication_id)} referenced"
    print(summary)


#########################################
def html_pages(path: str) -> List[str]:
    if os.path.isfile(path):
        return [path]
    pages = sorted(entry.path for entry in os.scandir(path) if entry.name.endswith(".html"))
    if not pages:
        errmsg = f"no html pages found in {path}"
        sys.exit(errmsg)
    return pages


def main():
    if len(sys.argv) != 2:
        usage = f"usage: {sys.argv[0]} <hgmd_gene_page.html | directory of them>"
        sys.exit(usage)
    path = sys.argv[1]
    if not os.path.exists(path):
        errmsg = f"{path} not found"
        sys.exit(errmsg)

    engine = db_engine()
    require_genes(engine)

    pages = html_pages(path)
    parsed_pages = [(page, ) + parse_html(page) for page in pages]
    # a page with no 'All <n> mutations in <gene>' heading is not an HGMD mutation page: the save
    # is a portal shell or otherwise unusable, and no amount of parsing will get data out of it
    unusable = [page for page, name, _ in parsed_pages if name is None]
    usable = [(page, name, parsed) for page, name, parsed in parsed_pages if name is not None]
    report = f"{len(usable)} of {len(pages)} page(s) parsed"
    print(report)

    with Session(engine) as session:
        # the whole symbol set first: an unresolvable one must stop the run before the first insert
        resolved = gene_ids(session, [name for _, name, _ in usable])
        for page, name, parsed in usable:
            print(f"gene: {name}  ({os.path.basename(page)})")
            store(session, resolved[name], parsed)

    if unusable:
        # non-zero exit, but only after the good pages are in: the loads are idempotent, so
        # re-running once these are saved again costs nothing
        names = "\n\t".join(os.path.basename(page) for page in unusable)
        errmsg = f"{len(unusable)} page(s) hold no mutation table and have to be saved again:\n\t{names}"
        sys.exit(errmsg)


#########################################
if __name__ == '__main__':
    main()
