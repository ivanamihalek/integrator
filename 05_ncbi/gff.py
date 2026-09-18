#! /usr/bin/python3
"""Reading the RefSeq gff and its assembly report - shared by the loaders in this directory."""

import gzip
from typing import Dict, Iterator, List
from urllib.parse import unquote


#########################################
def chromosome_names(path: str) -> Dict[str, str]:
    """RefSeq accession -> ucsc style chromosome name, for the assembled molecules only.

    Restricting to assembled molecules is the selection step of both loaders here: the alt loci
    and patch scaffolds repeat the same features in coordinates nothing downstream uses.
    """
    names: Dict[str, str] = {}
    with open(path) as inf:
        for line in inf:
            if line.startswith("#"):
                continue
            field = line.rstrip("\n").split("\t")
            if len(field) < 10 or field[1] != "assembled-molecule":
                continue
            names[field[6]] = field[9]
    return names


def attributes(column: str) -> Dict[str, str]:
    """The ninth gff column, percent-decoded, as a dictionary."""
    parsed: Dict[str, str] = {}
    for item in column.split(";"):
        key, _, value = item.partition("=")
        if key:
            parsed[key] = unquote(value)
    return parsed


def lines(path: str) -> Iterator[List[str]]:
    with gzip.open(path, "rt") as inf:
        for line in inf:
            if line.startswith("#"):
                continue
            field = line.rstrip("\n").split("\t")
            if len(field) == 9:
                yield field
