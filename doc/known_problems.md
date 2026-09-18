# Known problems

The state of the postgres/SQLModel rebuild as of 2026-09-18: what is blocked, what is loaded but
imperfect, and what has not been started. Row counts and stage numbering refer to
`.claude/plans/mysql_to_sqlmodel_postgres_migration.md`.


## Blocked on a source we cannot reach

### `uniprot_seqs` - no UCSC login
The table maps uniprot residue positions onto genomic coordinates, and needs the exon boundaries
of the UCSC `ensGene` track. `/home/ivana/.ucsc_mysql_conf` - the credentials file
`06_ucsc/13_get_gene_regions_from_ucsc.py` and `06_ucsc/11_refgene_regions_store.py` read - does
not exist on this machine, so the remote client cannot connect.

The canonical *sequence* is not affected: `uniprot_basic_infos.sequence` is taken straight from
the SwissProt flat file, so nothing here waits on `blastdbcmd` any more. Only the
transcript-to-genome mapping does.

### `omim_genemaps` and the inherited-metabolic-disorder marking - licensed download
`genemap2.txt` is behind an OMIM licence and is not on disk; `mim2gene.txt` (which is public) is
loaded, so `omim_genes` has 26,109 rows but the phenotype side has none.

`20_omim/10_mark_inherited_metabolic_disorders.py` therefore stays on the old MySQL stack. **Do
not delete it during the stage-8 cleanup**: it carries hand-curated `replace` and `fix`
dictionaries - disease-name corrections that exist nowhere else - and they have to be moved into
the new loader before the script goes.

### `10_metacyc` - licensed download
Same situation, no licensed copy on disk. Nothing loaded.

### `28_conservation` - exolocator is down
`exolocator.bioinfo.hr` answers HTTP 503. Stage 6 was always going to be a rewrite rather than a
conversion, but it cannot even be tested until the service returns.


## Loaded, with a caveat

### HGMD: 4 of 214 pages carry no data
`bbip1`, `bbs4`, `bbs5` and `cep78` were saved from the HGMD portal as javascript shells - the
mutation table is fetched by script, so the saved HTML holds nothing to parse. The loader reports
them and exits non-zero *after* loading the other 210 pages. **They need re-saving from the
portal.**

Seven further pages had been saved through the browser's *source viewer* rather than as HTML;
`19_hgmd/01_hgmd_from_html.py` reconstructs those itself, so they need no action.

The 210 usable pages give 38,711 mutations, which is exactly the total HGMD's own
"All N mutations" headlines declare.

### UCSC: 1,414 transcripts do not resolve to a gene
refGene names them `LOC…`, clone ids, or fusion constructs, and HGNC has no such symbol. They are
written to `/tmp/ucsc_unresolved_transcripts.tsv` and the loader exits non-zero, but only *after*
the 1,531,511 rows that do resolve are committed - the alternative would be to hold a whole
usable track hostage to its junk.

This is the one place where the fail-hard rule is relaxed to "report and exit non-zero" instead
of "exit before writing anything".

### OMIM is attached by Entrez id, not by symbol
353 of the 2,019 symbols in the OMIM download no longer resolve against HGNC (`QARS`→`QARS1`,
`SEPT2`→`SEPTIN2`, and so on). Resolution is by `ncbi_gene_id`, which survives a rename. Any new
loader for a source of that vintage should do the same.

### Orphanet excludes 32 records by selection
30 are entries Orphanet itself types as "Disorder-associated locus" rather than a gene - they are
outside the selection, which is `GENE_TYPES = ("gene with protein product", "Non-coding RNA")`.
The other 2 are snoRNA clusters whose HGNC ids have since been withdrawn.

That second case is why `hgnc_withdrawn` exists: a loader has to be able to tell "HGNC retired
this id" - a fact about HGNC, and a legitimate reason to skip a record - from "this id is
unknown", which means the pipeline is wrong somewhere and has to stop.

### KEGG: the copy on disk was five years stale
`/storage/databases/kegg/pathway_names.txt` dates from 2020 and names 89 pathways that the
current link file no longer has. The loader now reads the current REST listing
(`https://rest.kegg.jp/list/pathway/hsa`) and treats the downloads as a snapshot: pathways KEGG
has retired are deleted, and memberships are replaced rather than merged. The brite listing is
still the only source of the Metabolism / Carbohydrate-metabolism hierarchy, so if that hierarchy
is wanted, the brite file has to be re-downloaded too - a stale one stops the run.

### SMART holds three proteins, not the proteome
`07_smart/01_download_smart.py` has been run for ABCA4, MYO7A and USH2A only (92 features in
`smart_domains`). SMART is a public queue and the protocol is one sequence at a time, waiting for
each result, so the whole of `uniprot_basic_infos` would take days and should not be attempted
without asking EMBL first. Existing result files are never re-fetched, so widening the selection
later costs only the sequences that are new.

### `genes.ncbi_gene_id` and `genes.ensembl_gene_id` are indexed but not unique
HGNC maps 2 distinct symbols onto one Entrez id and 3 onto one Ensembl id. Any lookup by those
columns can therefore return more than one gene; the loaders that resolve this way build a dict
and accept last-one-wins.


## Not started

| | |
|---|---|
| `24_pdb` (stage 5) | sources are on disk - `/storage/databases/pdb` 939 MB, swissmodel 22 GB - but it is ~2,150 lines and three tables |
| `12_gnomad` (stage 7) | nothing downloaded; the partitioned `gnomad_freqs` table and its COPY path are designed but unbuilt. A multi-hour download |
| `14_capture_regions` | needs the CCDS and Agilent downloads |
| stage 8 cleanup | the superseded `.sql` files still exist; `integrator_utils/python/mysql.py` still holds both clients and should be cut down to the remote-UCSC one and renamed `remote_mysql.py`; the legacy password environment variables are still read; the README still describes the MySQL pipeline |
