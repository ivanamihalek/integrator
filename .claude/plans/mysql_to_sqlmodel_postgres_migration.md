# Rebuilding the `integrator` database on SQLModel + PostgreSQL

## Status, 2026-09-18

Executed through stage 4, with the parts of later stages whose sources were reachable.

| stage | state | rows loaded |
|---|---|---|
| 0 foundation + `01_hgnc` | **done** | 45,054 genes, 60,673 aliases, 5,291 withdrawn ids |
| 1 `05_ncbi` | **done** | 144,141 regulatory regions, 87,580 tss, 585 ncbi ids filled |
| 1 `03_uniprot` basic info | **done** | 20,133 entries (`uniprot_seqs` blocked, see below) |
| 2 `19_hgmd` | **done** | 38,711 mutations over 210 genes, 7,819 publications |
| 2 `20_omim` | **partial** | 26,109 mim2gene entries; genemap2 needs a licence |
| 2 `16_orphadata` | **done** | 11,645 diseases, 52,456 xrefs, 8,460 disease-gene links |
| 3 `06_ucsc` | **partial** | 1,531,511 refgene regions; `14_capture_regions` not started |
| 4 `09_kegg` | **done** | 372 pathways, 39,298 memberships |
| 4 `08_reactome` | **done** | 2,883 pathways, 2,899 edges, 139,614 participations |
| 4 `10_metacyc` | **blocked** | flat files are behind a licence |
| 5 `24_pdb` | not started | sources on disk |
| 6 `28_conservation` | **blocked** | exolocator answers 503 - the service is down |
| 7 `12_gnomad` | not started | download not attempted |
| 8 retire the old stack | not started | |

Twelve alembic revisions, `alembic check` clean, no orphan `gene_id` in any table, every loader
re-run to confirm it is idempotent.

Blocked, with the reason:

- **`uniprot_seqs`** wants exon coordinates from UCSC's `ensGene`, which is neither on disk nor
  reachable: `~/.ucsc_mysql_conf` does not exist, so the remote client cannot connect. The
  canonical sequence, which the legacy pipeline fetched with `blastdbcmd`, is in
  `uniprot_basic_infos` already - it is in the flat file.
- **`omim_genemaps`** and with it `10_mark_inherited_metabolic_disorders.py`: `genemap2.txt` is a
  licensed download. The hand-curated `replace`/`fix` dictionaries in that script are the
  curation that would be lost - leave the script in place until the file is available.
- **`10_metacyc`**, **`28_conservation`**: no source.
- **4 of the 214 HGMD pages** (`bbip1`, `bbs4`, `bbs5`, `cep78`) were saved as the portal's
  javascript shell and hold no mutation table; they have to be saved again. Seven others were
  saved from the browser's source viewer and are recovered by the loader.

## Premise

The `integrator` Postgres database is populated **from scratch, from primary sources**. The old
MariaDB instance and its contents are immaterial: nothing is dumped, copied or compared against
it. The legacy scripts matter only as *reference implementations* of the parsing and derivation
logic — how a MetaCyc dat file becomes an edge list, how a UCSC refGene row becomes an exon
region. Their table definitions are a starting suggestion, not a constraint.

Decisions taken:

1. **Prefixed tables in the `public` schema** — the six legacy MySQL databases become table-name
   prefixes, as `19_hgmd` already does.
2. **Alembic from the start** — the initial revision is generated before the first load.
3. **Re-derive everything** — every table is filled by running a loader against a primary source.
4. **`01_hgnc` first, and everything else fails hard without it** — the directories are renumbered
   so the run order is the numbering. Any loader outside `01_hgnc` that finds no `genes` table, or
   a gene it cannot resolve, prints the reason and exits.

This changes the shape of the work from "port 167 `search_db` calls" to "write ~15 loaders and
their models". Read-only analysis scripts still need their cursors replaced, but the write path —
where the `%`-formatted SQL and the MyISAM assumptions live — gets rewritten rather than
translated.

## Target architecture

One database, `integrator`, served by `@docker-compose.yml` (container `integrator-postgres`,
data in `/media/ivana/portable/postgres/integrator-pgdata`).

| legacy database | prefix | example |
|---|---|---|
| `blimps_*`, `monogenic_development` | none (core tables) | `genes`, `uniprot_basic_infos` |
| `identifier_maps` | `id_` | `id_translation`, `id_kegg` |
| `ucsc` | `ucsc_` | `ucsc_refgene_regions` |
| `gnomad` | `gnomad_` | `gnomad_freqs` |
| `reactome` | `reactome_` | `reactome_pathways` |
| (metacyc) | `metacyc_` | `metacyc_edges` |

Code layout mirrors `19_hgmd`:

```
integrator_utils/python/
    db.py             # engine, session, upsert helpers
    models_core.py    # tables shared across directories: gene, publications, uniprot_*, omim_*
    mysql.py          # only for remote read-only clients (UCSC public server); renamed in stage 8
alembic/              # versions/ holds one revision per stage
<nn>_<source>/
    <source>_models.py   # tables owned by that directory
    <nn>_*.py            # loaders
```

Directory-local models import the shared ones, so every `Field(foreign_key="genes.id")` resolves
inside a single `SQLModel.metadata`.

`integrator_utils/python/db.py` — the only new code every stage depends on:

```python
def db_engine(echo: bool = False) -> Engine              # dotenv from repo root, POSTGRES_* vars
def db_session() -> Iterator[Session]                    # contextmanager, commit/rollback
def upsert(session, model, values, conflict_on) -> int   # returns the row id
def upsert_many(session, model, rows, conflict_on) -> None
def insert_ignore(session, model, rows, conflict_on) -> None
def copy_from_tsv(engine, model, path, columns) -> int    # psycopg COPY, for the bulk loads
def require_genes(engine) -> None                         # genes table missing -> sys.exit
def gene_id(session, symbol) -> int                       # symbol unknown -> sys.exit
def gene_ids(session, symbols) -> Dict[str, int]          # any symbol unknown -> sys.exit
```

`upsert()` is `19_hgmd/01_hgmd_from_html.py:store_entry` lifted out and parameterized on the
conflict columns; `insert_ignore()` is the bulk branch of `store_publications`. `copy_from_tsv()`
is new and replaces every `mysqlimport` invocation in the comments of the old `.sql` files —
`COPY … FROM STDIN` is the only sane way to load gnomad-scale data. `require_genes()` and
`gene_id()` are the enforcement points for the rule in the next section.

## The gene hub — settle this before anything loads

A from-scratch build gets one thing the legacy schema never had: a single decision about gene
identity, made before the first row exists. `19_hgmd` created `gene` as `(id, name)` because that
was all it needed; every directory below joins on genes, so that table becomes the hub:

```python
class Gene(SQLModel, table=True):
    __tablename__ = "genes"
    id: Optional[int] = Field(default=None, primary_key=True)
    symbol: str = Field(index=True, unique=True)          # HGNC approved symbol, citext
    hgnc_id: Optional[str] = Field(default=None, index=True, unique=True)
    ncbi_gene_id: Optional[int] = Field(default=None, index=True, unique=True)
    ensembl_gene_id: Optional[str] = Field(default=None, index=True, unique=True)
    name: Optional[str] = None                             # approved full name
```

with aliases and previous symbols in a companion `gene_aliases` table
(`gene_id`, `alias`, `alias_type`) rather than the comma-joined `blob` columns of
`@01_hgnc/09_hgnc_table.sql` — every legacy consumer splits those strings back apart anyway.

**HGNC is the authority for gene identity.** The directory numbering now says so: `01_hgnc` runs
first and is the only writer of `genes`. Every other directory resolves symbols to `gene_id` and
*never* inserts into `genes`.

**The rule, enforced in `db.py` and applied in every directory except `01_hgnc`:**

1. **No `genes` table → error and exit.** Loaders do not call `SQLModel.metadata.create_all()`;
   schema creation belongs to Alembic alone. Left in, `create_all` would quietly conjure an empty
   `genes` and every subsequent lookup would fail one row at a time instead of once, up front.
   `19_hgmd/01_hgmd_from_html.py:274` has that call today and loses it in stage 0.
2. **Symbol not in `genes` → error and exit.** Not a warning, not a skipped row, not an
   autovivified hub entry. An unknown symbol means either HGNC has not been loaded, or the source
   uses a name HGNC does not approve — both are facts about the pipeline that have to be fixed in
   the loader, and both are invisible if the row is merely skipped.

```python
def require_genes(engine: Engine) -> None:
    if not inspect(engine).has_table("genes"):
        errmsg = "the 'genes' table does not exist - run 01_hgnc first (alembic upgrade head)"
        sys.exit(errmsg)


def gene_id(session: Session, symbol: str) -> int:
    gene = session.exec(select(Gene).where(Gene.symbol == symbol)).first()
    if gene is None:
        errmsg = f"gene symbol '{symbol}' not found in 'genes' - not an HGNC approved symbol?"
        sys.exit(errmsg)
    return gene.id
```

"Present in `genes`" means: approved symbol first, then `gene_aliases` (previous symbols and
synonyms, which is what that table is for), and only when both miss does the loader exit. An
alias hit is resolution, not a fallback — HGNC published that alias precisely so a source using
the old name can be matched.

`require_genes()` is called once, right after `db_engine()`, before any parsing — a loader that
cannot resolve genes should not spend forty minutes reading a dat file first. Bulk loaders call
`gene_ids()` on the whole symbol set up front and exit on the first unknown one, listing all of
them, so a 214-page or 24-chromosome run fails in seconds rather than partway through.

This is what keeps a from-scratch rebuild from re-acquiring the identifier mess the old
`id_translation`/`genes`/`uniprot_basic_infos` triangle had.

One caveat to state now rather than discover at 3 a.m. during the NCBI load: some sources
legitimately carry records outside the HGNC set — `gene_info` includes withdrawn, pseudo- and
non-coding entries, UCSC carries every transcript it knows, MetaCyc is not human-only. "Unknown
symbol → exit" applied to the raw file would make those loaders unable to run at all. The rule
holds; what those loaders need is an explicit, stated *selection* step first ("human, current,
protein-coding"), after which any remaining unresolvable symbol is a genuine error and exits. The
selection criterion goes in the loader's docstring, so a skipped record is a documented decision
rather than a silent swallow. `05_ncbi`, `06_ucsc`, `10_metacyc` and `12_gnomad` are the four that
need one.

`05_ncbi/01` and the other identifier loaders *update* `genes` (filling `ncbi_gene_id`,
`ensembl_gene_id`) rather than insert — an `UPDATE … WHERE symbol = …` that touches zero rows is
the same error as an unknown symbol and exits the same way.

Consequence for `19_hgmd`: its `Gene(id, name)` is replaced by the above, its `store_gene()`
becomes `gene_id()`, and its `create_all` goes. Its 42 TTC8 rows are re-loaded after the change —
and TTC8 has to be in `genes` first, which means `01_hgnc` now runs before the stage 0 round trip
can be verified at all.

## Source inventory

Checked on this machine, since "re-derive" makes source availability the actual schedule driver:

| directory | primary source | status |
|---|---|---|
| `01_hgnc` | `/storage/databases/hgnc/` | on disk, 8 MB (`hgnc2refseq.tsv`, `hgnc_name_res.tsv`) |
| `05_ncbi` | `/storage/databases/ncbi/`, `gene_info.gz` | on disk, 4.9 GB |
| `03_uniprot` | `/storage/databases/uniprot/` | on disk, 3.5 GB (`uniprot_basic_info.tsv`, `uniprot_seq.tab`, `uniprot_hgnc.tsv`) |
| `19_hgmd` | `/storage/databases/hgmd-ird-related/` | on disk, 80 MB — **214 html pages** |
| `20_omim` | `/storage/databases/omim/` | on disk, 1.4 MB (`mim2gene.tsv`, `mim2hgnc.tsv`) |
| `06_ucsc` | `/storage/databases/ucsc/`, UCSC public MySQL | on disk, 45 GB + remote |
| `09_kegg` | `/storage/databases/kegg/`, REST | on disk, 1.9 MB |
| `24_pdb` | `/storage/databases/pdb/`, `swissmodel/` | on disk, 939 MB + 22 GB |
| `08_reactome` | reactome download | **missing — download first** |
| `10_metacyc` | `/storage/databases/meta_cyc/` | **missing — download first** |
| `16_orphadata` | orphadata xml | **missing — download first** |
| `12_gnomad` | gnomad vcfs | **missing — large download** |
| `14_capture_regions` | CCDS, Agilent bed | **missing — download first** |
| `28_conservation` | exolocator REST | remote, availability unverified |

Directories whose source is on disk are scheduled first; the missing ones get a download step as
their first task, and `12_gnomad` is scheduled last of the real work because the download alone
is a multi-hour affair.

## Stage 0 — foundation and the gene hub

```txt
@integrator_utils/python/mysql.py      (the API being replaced)
@19_hgmd/hgmd_models.py                (the model style to follow)
@19_hgmd/01_hgmd_from_html.py          (db_engine + upsert code to lift)
@01_hgnc/09_hgnc_table.sql
@01_hgnc/06_hgnc_tables_in_identifier_maps.py.sql
@01_hgnc/10_add_hgnc_synonyms_to_genes.py
@docker-compose.yml
@pyproject.toml
```

`01_hgnc` belongs here rather than in stage 1: under the fail-hard rule nothing else can be run,
let alone verified, until `genes` is populated.

1. `integrator_utils/python/db.py` with the API above.
2. `integrator_utils/python/models_core.py` with `Gene`, `GeneAlias`, `Publication` (both moved
   out of `19_hgmd/hgmd_models.py`), plus `UniprotBasicInfo`, `UniprotSeq`, `OmimGenemap`,
   `Disease` — columns taken from `@03_uniprot/02_uniprot_table.sql`,
   `@03_uniprot/04_uniprot_seq_table.sql`, `@20_omim/02_omim_genemap_table.sql`, with `blob` →
   `text`, comma-joined lists → `ARRAY(Text)` or companion tables, symbols → `citext`.
3. `uv add alembic`, `alembic init`, configure `env.py` to take the URL from `db_engine()` and
   `target_metadata = SQLModel.metadata`; autogenerate revision 1 from `models_core` + `19_hgmd`.
4. `01_hgnc` loader: fill `genes` and `gene_aliases` from the HGNC download — the first and only
   writer of the hub. `10_add_hgnc_synonyms_to_genes.py` (which patched synonyms into an existing
   `genes` table) is folded in and deleted. This loader is the one place that creates gene rows,
   so it is also the only one that does not call `require_genes()`.
5. Point `19_hgmd` at `db.py`/`models_core.py`, drop its local `Gene`/`Publication`/`db_engine`
   and its `create_all`, add `require_genes()` + `gene_id()`. Then: drop the database in the
   container, `alembic upgrade head`, load HGNC, re-load TTC8, confirm 18 / 8 / 13 / 3 mutations,
   33 publications, 54 links. Also confirm the two failure modes on purpose — run the HGMD loader
   against a database with no `genes` table, and against one where TTC8 is absent; both must exit
   non-zero with the message naming the cause, having written nothing. That round trip — empty
   database to verified counts and verified failures through Alembic — is the template every
   stage below repeats.

**Every stage from here on ends with `alembic revision --autogenerate`, review of the generated
diff, and `alembic upgrade head`.** Autogenerate does not notice column renames or type changes
it cannot infer; the review step is not optional.

## Stage 1 — identifiers on top of the hub: `05_ncbi`, `03_uniprot`

```txt
@05_ncbi/01_id_gene_name_map.py
@05_ncbi/03_regulatory_regions_from_gff.py
@05_ncbi/05_tss_from_gff.py
@03_uniprot/02_uniprot_table.sql
@03_uniprot/04_uniprot_seq_table.sql
@03_uniprot/11_add_uniprot_alt_names_to_id_table.py
@03_uniprot/12_cleanup_alias_names_in_uniprot.py
@03_uniprot/14_gene_coordinates_patch.py
@integrator_utils/python/db.py
@integrator_utils/python/models_core.py
```

- `05_ncbi/01` currently writes `/storage/databases/ncbi/gene_id_gene_name.tsv`; it becomes the
  step that fills `genes.ncbi_gene_id`. `03_regulatory_regions_from_gff.py` and `05_tss_from_gff.py`
  emit BED today — give them `ncbi_regulatory_regions` and `ncbi_tss` tables, keeping BED output
  as an option, since the downstream consumers are bedtools pipelines.
- `03_uniprot` fills `uniprot_basic_infos` and `uniprot_seqs`, resolving `gene_name` → `gene_id`.
  The three `*_patch*.py` / `*cleanup*.py` scripts are hand-fix scripts against the old database:
  read them for the corrections they encode, fold those into the loader as data rules, and drop
  the scripts. **This is the one place where a from-scratch rebuild can silently lose work** —
  hand-curated fixes are not in any primary source.
- `08_cDNA_seq_to_uniprot.py`, `09_check_uniprot_ens_mismatch_cases.py`,
  `10_align_uniprot_and_ensembl_pepseqs.py` are read-only QC; convert their cursors last.

## Stage 2 — `19_hgmd` at full scale, `20_omim`, `16_orphadata`

```txt
@19_hgmd/01_hgmd_from_html.py
@19_hgmd/hgmd_models.py
@20_omim/02_omim_genemap_table.sql
@20_omim/10_mark_inherited_metabolic_disorders.py
@16_orphadata/02_diseases_and_xrefs.py
@26_rare_diseases/01_screen_proficiency_table.sql
```

`/storage/databases/hgmd-ird-related/` holds 214 saved gene pages — add a batch driver that walks
a directory, and this stage populates the mutation tables for the whole IRD gene set in one run.
It is also the first real test of the gene hub, and the first place the fail-hard rule bites at
scale: the batch driver reads the gene name out of all 214 pages and resolves the whole set
through `gene_ids()` *before* storing anything, so a symbol HGMD spells differently from HGNC
stops the run in seconds with the full list of offenders — rather than 180 pages in, with a
half-loaded database to clean up. Expect this list to be non-empty on the first run; HGMD is not
scrupulous about tracking symbol changes. Each entry is then either an HGNC alias (the loader
resolves through `gene_aliases`) or a genuine rename to record.

`20_omim` fills `omim_genemaps` and `diseases`; `16_orphadata/02_diseases_and_xrefs.py` parses
disease xrefs and stores nothing today — give it `diseases` + a `disease_xrefs` through table
(same shape as HGMD's publication links). Orphadata xml needs downloading first.

## Stage 3 — `06_ucsc`, `14_capture_regions`

```txt
@06_ucsc/01_create_ucsc_db.sql
@06_ucsc/11_refgene_regions_store.py
@06_ucsc/13_get_gene_regions_from_ucsc.py
@06_ucsc/15_ucsc_ensgene_to_regions.py
@14_capture_regions/50_ucsc_ensgene_to_bed.py
@14_capture_regions/54_ccds_to_bed.py
```

**These are MySQL clients of UCSC's public server** (`genome-mysql.cse.ucsc.edu`, via
`/home/ivana/.ucsc_mysql_conf`, which no longer exists on this machine and has to be recreated).
The remote half stays on MySQL — UCSC serves MySQL and that is not ours to change. Only the local
write half migrates: `11_refgene_regions_store.py` opens two connections and becomes "read UCSC
over MySQL, write Postgres through SQLModel". `13_get_gene_regions_from_ucsc.py` is remote-read
plus file output, so it keeps importing `mysql.py` with a comment saying why.

Region tables (`ucsc_refgene_regions`, `ucsc_ensembl_gene_regions`) should use `int4range`/
`int8range` with GiST indices instead of paired `start`/`end` integers: every consumer in
`14_capture_regions` and `12_gnomad/64_map_variants_to_regions.py` is doing interval overlap by
hand today, and `&&` with a GiST index replaces that logic outright.

## Stage 4 — pathways: `09_kegg`, `08_reactome`, `10_metacyc`

```txt
@09_kegg/15_kegg_tables_in_identifier_maps.sql
@09_kegg/18_kegg_id2pthwy.py
@08_reactome/16_pathways_from_reactome.py
@10_metacyc/23_metacyc_edges.sql
@10_metacyc/22_metacyc_protein_id_mapper.py
@10_metacyc/24_human_enzyme_edges.py
@10_metacyc/25_metacyc_enzrxn_pathways.py
@10_metacyc/27_metacyc_pwy_hierarchy.py
@10_metacyc/31_metacyc_enzrxn_substrate_edges.py
@10_metacyc/33_metacyc_enrxn_regulators.py
```

KEGG's source is on disk; Reactome and MetaCyc have to be downloaded first. One
`metacyc_models.py` covers the seven tables the nine scripts share (`metacyc_genes`,
`metacyc_proteins`, `metacyc_enzrxns`, `metacyc_reactions`, `metacyc_pathways`,
`metacyc_compounds`, `metacyc_edges`). The `*_inspector.py` scripts (`26`, `30`, `32`) are
read-only reporting — plain `select()` conversions with no write path to test.

## Stage 5 — `24_pdb`

```txt
@24_pdb/14_pdb_ligands_table.sql
@24_pdb/15_pdb_uniprot_ec_maps_table.sql
@24_pdb/16_model_elements_table.sql
@24_pdb/10_find_pdb_models_for_iem_proteins.py
@24_pdb/12_compile_model_from_swissmodel.py
@24_pdb/18_compile_model_elements.py
@24_pdb/20_compile_model.py
@24_pdb/08_build_on_backbone.py
@24_pdb/compile_model_test.py
```

~2,150 lines, the largest code volume in the repo, but the database use is shallow: three tables
and lookups keyed on uniprot id; most of the length is structure manipulation that never touches
the database. Sources (`/storage/databases/pdb`, `/storage/databases/swissmodel`) are on disk.
`compile_model_test.py` is the repo's only test — keep it green throughout.

## Stage 6 — `28_conservation`

```txt
@28_conservation/01_exolocator_pep_alignments.py
@28_conservation/02_alignment_cleanup.py
@28_conservation/04_conservation_vert_orthologues.py
@28_conservation/06_conservation_distant_homologues.py
```

Read-only consumers of `uniprot_seqs` and `genes`, plus the exolocator REST service
(`exolocator.bii.a-star.edu.sg`) — check that it still answers before budgeting this stage; if it
is gone, the alignments have to come from elsewhere and this stage becomes a rewrite rather than
a conversion.

## Stage 7 — `12_gnomad`

```txt
@12_gnomad/41_gnomad_mysql_tables.py
@12_gnomad/42_gnomad_freq_to_tsv.py
@12_gnomad/45_gnomad_freq_table_indices.py
@12_gnomad/46_find_hotspot_regions.py
@12_gnomad/48_map_gnomad_freq_to_hotspots_table.py
@12_gnomad/64_map_variants_to_regions.py
```

Last because the download is the long pole and nothing else depends on it. The only stage where
the schema genuinely changes shape:

- `41_gnomad_mysql_tables.py` creates **24 identical `freqs_chr_*` tables in a loop**. Replace
  with **one partitioned table** `gnomad_freqs`, `partition by list (chrom)`, 24 partitions: same
  per-chromosome I/O locality, one model, one set of index declarations, and cross-chromosome
  queries that work. SQLModel cannot express partitioning, so the parent DDL and the partition
  loop live in the models file and go into an Alembic revision as raw `op.execute()`.
- `42_gnomad_freq_to_tsv.py` already writes TSV for bulk import — point it at `copy_from_tsv()`.
- `45_gnomad_freq_table_indices.py` mostly dissolves into `Field(index=True)`, except that
  indices must be created *after* the bulk load at this scale; keep that as an explicit script
  step, not an Alembic side effect.
- `12_gnomad/original_exac_processing/` and `12_gnomad/tests/` are historical: leave them alone
  and say so in the README.

## Stage 8 — retire the old stack

1. Delete the `.sql` files whose tables now live in models: `@03_uniprot/02_uniprot_table.sql`,
   `@03_uniprot/04_uniprot_seq_table.sql`, `@01_hgnc/09_hgnc_table.sql`,
   `@01_hgnc/06_hgnc_tables_in_identifier_maps.py.sql`, `@09_kegg/15_kegg_tables_in_identifier_maps.sql`,
   `@10_metacyc/23_metacyc_edges.sql`, `@20_omim/02_omim_genemap_table.sql`,
   `@24_pdb/14_pdb_ligands_table.sql`, `@24_pdb/15_pdb_uniprot_ec_maps_table.sql`,
   `@24_pdb/16_model_elements_table.sql`, `@26_rare_diseases/01_screen_proficiency_table.sql`,
   `@06_ucsc/01_create_ucsc_db.sql`. Git history keeps them.
2. Reduce `integrator_utils/python/mysql.py` to `connect_to_mysql` + `search_db` for the UCSC
   client and rename it `remote_mysql.py`.
3. Drop `COOKIEMONSTER_PASSWORD`, `BLIMPS_DATABASE_PASSWORD` and the
   hostname-based development/production switch in `connect()`. `.env` + `POSTGRES_*` is the only
   configuration left; `.env` is now gitignored.
4. Rewrite `@README.md`: the stack, `docker compose up -d`, `alembic upgrade head`, and the
   directory run order — which is exactly the stage order above, since a from-scratch rebuild has
   to be reproducible by someone who has only the repo and the source downloads.
5. `22_ensembl/01_ensembl_phenotypes_table.sql`, `02_biogrid/06_homo_sapiens.sql`,
   `40_expression/01_bgee.sql`, `18_mesh`, `30_incidence` have only perl/bash behind them. Mark
   them legacy in the README; convert only if those pipelines are to be revived.

## What a from-scratch rebuild changes, beyond the driver

- **Hand patches are lost unless folded in.** `03_uniprot/07_uniprot_seq_patch_by_hand.py`,
  `12_cleanup_alias_names_in_uniprot.py`, `14_gene_coordinates_patch.py` and
  `01_hgnc/10_add_hgnc_synonyms_to_genes.py` exist to fix data *after* loading. Read each one for
  the corrections it encodes and move those into the loader before deleting it.
- **Case sensitivity.** MySQL's default collation is case-insensitive; Postgres is not. Every
  gene-symbol lookup changes behaviour. `citext` on `genes.symbol`, `gene_aliases.alias` and the
  uniprot/omim symbol columns, decided at model-definition time — retrofitting it later means a
  migration plus a reload.
- **No MyISAM.** The old loaders relied on non-transactional inserts surviving an interrupt.
  Loaders now need explicit per-chunk commits, or a resumable "already loaded" check — for
  gnomad-scale loads, per-chromosome transactions.
- **Comma-joined `blob` columns become arrays or child tables.** `previous_symbols`, `synonyms`,
  `pubmed_ids`, `refseq_ids`, `uniprot_ids` in `@01_hgnc/09_hgnc_table.sql` are all
  split-on-comma at every use site.
- **Interval columns become ranges.** See stage 3; this removes hand-written overlap logic from
  three directories.
- **Every loader is idempotent**, upserting on a natural key the way `19_hgmd` upserts on HGMD
  accession. Re-running a loader must not change row counts — that is what makes a from-scratch
  rebuild repeatable rather than a one-shot.
- **Loaders fail, they do not cope.** Missing `genes` table or unresolvable symbol → exit. This
  trades convenience for a database in which every `gene_id` means what it says, and it is only
  affordable because the loaders are idempotent: the cost of a hard exit is re-running a loader
  that has already been made safe to re-run.

## Verification

Per stage, in order:

1. `uv run python -c "import <dir>.<source>_models"` plus `alembic upgrade head` against a scratch
   database — catches duplicate `__tablename__`, unresolved foreign keys and bad autogenerated
   diffs before any data moves.
2. Run the loader on the real source; check row counts against the source itself (file line
   counts, the `All n mutations` headers HGMD prints, the record counts the dat files declare).
   There is no old database to compare with — the source is the only reference.
3. Re-run the loader; counts must not change.
4. Foreign-key sanity: no orphan `gene_id` — which under the fail-hard rule should be structurally
   impossible, so a violation here means a loader bypassed `gene_id()` and that is the bug to fix.
5. The two failure modes, deliberately provoked once per directory: run the loader against a
   database with no `genes` table, and against one missing a symbol the source needs. Both exit
   non-zero, name the cause, and leave the database untouched.
6. For the read-only analysis scripts, run old and new against the same input and diff the output
   files — they write TSV/BED, so the diff is exact.

## Effort and sequencing

| stage | directories | scripts | LOC | source status |
|---|---|---|---|---|
| 0 | `integrator_utils`, `01_hgnc`, `19_hgmd` | 6 | ~450 new | hgnc on disk |
| 1 | `05_ncbi`, `03_uniprot` | 9 | ~800 | on disk |
| 2 | `19_hgmd`, `20_omim`, `16_orphadata` | 4 | ~300 | on disk / orphadata to download |
| 3 | `06_ucsc`, `14_capture_regions` | 5 | ~370 | on disk + UCSC remote; CCDS to download |
| 4 | `09_kegg`, `08_reactome`, `10_metacyc` | 11 | ~1,170 | kegg on disk; two to download |
| 5 | `24_pdb` | 8 | ~2,150 | on disk |
| 6 | `28_conservation` | 4 | ~415 | remote, verify first |
| 7 | `12_gnomad` | 7 | ~1,700 | large download |
| 8 | cleanup, README | — | — | — |

Stage 0 blocks everything, now literally and not just logically: without `01_hgnc` loaded, every
other loader exits on its first line. Within 1-7 the order is set by source availability rather
than by dependency, so the download-first directories can be started in parallel with the ones
whose data is already here — the only hard sequencing left in the plan is "HGNC first".
