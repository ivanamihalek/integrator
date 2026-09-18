# Dockerized `integrator` Postgres + HGMD HTML loader

Files involved:
@docker-compose.yml (new — Postgres service)
@pyproject.toml (dependencies to fill in)
@19_hgmd/hgmd_models.py (new — SQLModel table definitions)
@19_hgmd/01_hgmd_from_html.py (new — parser/loader script)
@.env (POSTGRES_PASSWD / POSTGRES_HOST / POSTGRES_PORT, read via dotenv)
@.claude/CLAUDE.md (Python formatting rules to follow)

## Context

The `integrator` repo collects biomedical data into local databases; so far everything has
targeted MariaDB (`integrator_utils/python/mysql.py`), and there is no container definition in
the repo. We want a Postgres instance with a database named `integrator`, running in Docker with
its data on the portable drive, plus a first consumer: a loader that scrapes HGMD per-gene
mutation pages (HTML saved locally, e.g. `/home/ivana/scratch/ttc8.html`) into four tables via
the SQLModel ORM.

Findings from the example file (`/home/ivana/scratch/ttc8.html`, gene TTC8, 42 mutations):

- Gene name lives in `<h3>All 42 mutations in <a href="gene.php?gene=TTC8">TTC8</a></h3>`.
- Each section is an `<h3>` like `Splicing : 8 mutations [back to top]` followed by exactly one
  `<table>`. Sections present: Missense/nonsense, Splicing, Small deletions, Small insertions,
  Gross deletions (Gross insertions appears in other genes, same layout as gross deletions).
- Each table = one junk `<tr>` (column-toggle form/script) + one `<th>` header row + one row per
  mutation. Row counts match the section headline exactly.
- Accession is inside `<input type="hidden" name="acc" value="CM102429">`.
- Column sets differ per section (see mapping below); the VCF column is present but
  `style="display:none;"` — still parseable.
- "Reference" cells hold a primary publication anchor
  (`<a href="https://pubmed.ncbi.nlm.nih.gov/19797195/">Harville (2010) J Med Genet 47, 262</a>`)
  optionally followed by `<span class="td">`-wrapped extra anchors annotated in brackets
  (`[Additional report]`, `[Bardet-Biedl syndrome]`). One row can cite up to 4 publications here,
  and the same publication (e.g. PMID 25525159) recurs across many rows and sections — hence the
  separate `publications` table and the through tables.
- "Extra information" cells hold `<span class="gen" title="Chr14:...">hg38</span>` style spans
  (no link → skipped) and real links: `dbSNP` (16 rows) and `gnomAD` (7 rows) in this file, with
  the useful text in `title=` (on the `<a>` for dbSNP, on the inner `<span>` for gnomAD).
- Missing values show up as `Not yet available` → stored as NULL.

Decisions already taken: upsert on HGMD accession (re-runnable), dependencies declared in
`pyproject.toml` and installed with `uv`, database created by compose (not by the script).

All Python written here follows `.claude/CLAUDE.md`: horizontally aligned argument lists within
120 columns, message strings and dict literals defined on their own lines before the call,
f-strings only, `typing` hints throughout, no dangling lone brackets.

## 1. `docker-compose.yml` (repo root)

```yaml
services:
  integrator-db:
    image: postgres:latest
    container_name: integrator-postgres
    restart: unless-stopped
    environment:
      POSTGRES_DB: integrator
      POSTGRES_USER: postgres
      POSTGRES_PASSWORD: ${POSTGRES_PASSWD}
      PGDATA: /var/lib/postgresql/data/pgdata
    ports:
      - "${POSTGRES_PORT:-5432}:5432"
    volumes:
      - /media/ivana/portable/postgres:/var/lib/postgresql/data
    shm_size: 256mb
    healthcheck:
      test: ["CMD-SHELL", "pg_isready -U postgres -d integrator"]
      interval: 10s
      timeout: 5s
      retries: 5
```

Notes:
- Compose auto-reads `./.env` for `${POSTGRES_PASSWD}` / `${POSTGRES_PORT}` interpolation, so no
  `env_file:` is needed and no secret is hardcoded.
- `PGDATA` is pinned to a subdirectory of the mount: required for bind mounts (init refuses a
  non-empty dir such as `lost+found`) and keeps the file working with `postgres:18`, whose image
  moved its default data dir.
- `/media/ivana/portable` is currently not readable (I/O error) — the drive must be mounted, and
  formatted as a POSIX filesystem (ext4/xfs); exFAT/NTFS will not work for a PGDATA.

## 2. Dependencies (`pyproject.toml`)

Fill the empty `dependencies = []` with: `sqlmodel`, `psycopg[binary]`, `python-dotenv`,
`beautifulsoup4`, `lxml`. Install into the existing `.venv` via `uv sync`.

## 3. `19_hgmd/hgmd_models.py` — SQLModel definitions

```python
class Gene(SQLModel, table=True):          # __tablename__ = "gene"
    id: Optional[int] = Field(default=None, primary_key=True)
    name: str = Field(index=True, unique=True)


class Publication(SQLModel, table=True):   # __tablename__ = "publications"
    id: Optional[int] = Field(default=None, primary_key=True)
    reference: str                          # "Harville (2010) J Med Genet 47, 262"
    pubmed_id: Optional[int] = Field(default=None, index=True, unique=True)
```

A non-table mixin carries the columns shared by all four mutation tables:
`gene_id: int = Field(foreign_key="gene.id", index=True)`, `accession: str` (unique, indexed),
`variant_class`, `phenotype`, and
`xrefs: Optional[dict] = Field(default=None, sa_column=Column(JSONB))`.

| table | extra columns |
|---|---|
| `hgmd_missense_nonsense` | `codon_change`, `aa_change`, `hgvs_nucleotide`, `hgvs_protein`, `vcf` |
| `hgmd_splicing` | `splicing_mutation`, `hgvs_nucleotide`, `vcf` |
| `hgmd_small_indels` | `variant_type` ("deletion"/"insertion"), `description`, `hgvs_nucleotide`, `hgvs_protein`, `vcf` |
| `hgmd_gross_indels` | `variant_type`, `dna_level`, `description`, `hgvs_nucleotide`, `hgvs_protein`, `vcf` |

`variant_type` is what distinguishes the merged deletion/insertion sections. All text columns are
`Optional[str]`. No `reference` text column on the mutation tables — it is fully normalized into
`publications` plus the link tables below.

Four through tables (one per mutation table, so the FKs stay real rather than polymorphic), all
sharing a mixin and named `<mutation table>_publication`, e.g.

```python
class MissenseNonsensePublication(SQLModel, table=True):
    __tablename__ = "hgmd_missense_nonsense_publication"
    entry_id: int = Field(foreign_key="hgmd_missense_nonsense.id", primary_key=True)
    publication_id: int = Field(foreign_key="publications.id", primary_key=True)
    is_primary: bool = False            # first anchor in the cell
    annotation: Optional[str] = None    # bracketed note, e.g. "Additional report"
```

plus `hgmd_splicing_publication`, `hgmd_small_indels_publication`,
`hgmd_gross_indels_publication`. Composite PK `(entry_id, publication_id)` makes relinking on
re-run idempotent.

The module also exposes the header→field mapping used by the parser, keyed by normalized header
text, e.g. `"hgmd accession" → accession`, `"hgmd codon change" → codon_change`,
`"hgmd amino acid change" → aa_change`, `"hgvs (nucleotide)" → hgvs_nucleotide`,
`"hgvs (protein)" → hgvs_protein`, `"vcf chrom pos ref/alt" → vcf`,
`"variant class" → variant_class`, `"reported phenotype" → phenotype`,
`"reference" → publications (handled separately)`, `"extra information" → xrefs`,
`"hgmd splicing mutation" → splicing_mutation`,
`"hgmd deletion"/"hgmd insertion"/"description" → description`, `"dna level" → dna_level`.
Unrecognized headers are warned about and skipped, so HGMD layout drift does not crash the load.

## 4. `19_hgmd/01_hgmd_from_html.py` — loader

Follows the repo script style (`#! /usr/bin/python3`, `main()`, `if __name__ == '__main__'`),
takes the HTML path as the single CLI argument (usage/exit if missing or unreadable).

1. `load_dotenv()` from the repo root; build
   `postgresql+psycopg://postgres:{POSTGRES_PASSWD}@{POSTGRES_HOST}:{POSTGRES_PORT}/integrator`,
   `create_engine(...)`, then `SQLModel.metadata.create_all(engine)` — creates any missing table
   (the DB itself is created by compose).
2. Parse with BeautifulSoup/lxml:
   - gene name from the `<h3>` matching `All \d+ mutations` (text of its `<a>`, fallback to the
     trailing token);
   - for each section title in {Missense/nonsense, Splicing, Small deletions, Small insertions,
     Gross deletions, Gross insertions} find the matching `<h3>` and its `find_next("table")`;
     missing sections are simply skipped;
   - header row = first `<tr>` with `<th>`; build a column-index → field list from the mapping;
   - data rows = following `<tr>`s containing `<td>`, skipping the toggle-form row;
   - cell → text via `get_text(" ", strip=True)` with whitespace collapsed; `""` and
     `"Not yet available"` → `None`; accession read from the hidden `acc` input;
   - `xrefs`: for each `<a href>` in the extra-information cell, `{label: {"url": href,
     "title": inner-span title or the anchor's title}}`; anchors are the only entries kept, so
     the link-less `hg38`/`hg19`/`COM`/`CpG` spans drop out. Empty dict → `None`. (Repeated
     labels in one cell: last one wins.)
   - references: from the Reference cell collect every `<a href>` in document order as
     `(reference_text, pubmed_id, is_primary, annotation)`; `pubmed_id` from
     `pubmed.ncbi.nlm.nih.gov/(\d+)` (`None` if the anchor is not a PubMed link), `is_primary`
     for the first anchor, `annotation` from a trailing `[...]` in the enclosing span.
3. Gene lookup/insert: `select(Gene).where(Gene.name == name)`, create if absent, reuse `gene.id`
   as the `gene_id` foreign key on every row.
4. Publication upsert: dedup by `pubmed_id` within the file, then
   `insert(Publication).on_conflict_do_nothing(index_elements=["pubmed_id"])` followed by a
   select to resolve ids (first-seen `reference` text wins; HGMD reformats the same citation
   slightly between rows). Anchors without a PubMed id get a plain insert with `pubmed_id=None`.
5. Write each mutation row with `sqlalchemy.dialects.postgresql.insert(...).on_conflict_do_update(
   index_elements=["accession"], set_=<all non-PK columns>).returning(id)`, so a re-run of the
   same or an updated dump updates in place instead of duplicating and still yields the row id.
6. Link rows: `insert(<through table>).on_conflict_do_update(
   index_elements=["entry_id", "publication_id"], set_={is_primary, annotation})`.
7. Print a per-table summary (`hgmd_missense_nonsense: 18 rows`, …, `publications: N rows`) at
   the end.

## Verification

```bash
cd /home/ivana/academia/projects/pypeworks/integrator
uv sync
docker compose up -d && docker compose ps          # healthcheck must go healthy
.venv/bin/python 19_hgmd/01_hgmd_from_html.py /home/ivana/scratch/ttc8.html
```

Expected counts for TTC8: `gene` 1 row (TTC8), `hgmd_missense_nonsense` 18, `hgmd_splicing` 8,
`hgmd_small_indels` 13 (11 deletions + 2 insertions), `hgmd_gross_indels` 3; every mutation row
has at least one publication link, and `publications` holds strictly fewer rows than the total
number of links (PMID 25525159 alone is cited by several rows).

```bash
docker exec -it integrator-postgres psql -U postgres -d integrator \
  -c "select count(*) from hgmd_missense_nonsense;" \
  -c "select accession, xrefs from hgmd_missense_nonsense where xrefs is not null limit 3;" \
  -c "select m.accession, p.pubmed_id, p.reference, l.is_primary, l.annotation
        from hgmd_missense_nonsense m
        join hgmd_missense_nonsense_publication l on l.entry_id = m.id
        join publications p on p.id = l.publication_id
       where m.accession = 'CM102429';"
```

Spot checks: accession `CM102429` has `codon_change=TGG-TAG`, `aa_change=Trp41Term`,
`hgvs_nucleotide=c.122G>A`, `vcf=Chr14:88839459G/A`, `variant_class=DM`, and exactly two
publications — PMID 19797195 (`is_primary=true`) and 25525159 (`annotation='Additional report'`);
`CS111278` (splicing) has a `dbSNP` xref pointing at `rs139234943`; gross deletion `CG2217184` has
`dna_level=gDNA`, `description=ex, 7-8`, and NULL HGVS columns. Then re-run the loader and confirm
all counts — mutation rows, `publications`, and link rows — are unchanged (upsert works).
