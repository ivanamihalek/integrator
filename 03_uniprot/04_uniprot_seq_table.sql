-- to create the table:
-- bash> mysql -u root blimps_environment < 04_uniprot_seq_table.sql
-- to load dat, start mysql in local infile mode:
-- mysql --local-infile -uroot
-- then load a cleaned biogrid table from mysql shell
-- mysql>  load data local infile 'uniprot.csv' into table uniprot_basic_infos ignore 1 lines;
-- the horrible table name is to comply with rails pluralization convention
-- simpler -- not the matching file and table name
-- mysqlimport --local -u root -p blimps_development  uniprot_basic_infos.csv

CREATE TABLE `uniprot_seqs` (
`uniprot_id` varchar(40)   NOT NULL,
`ensembl_transcript_id`  varchar(40)  DEFAULT NULL,
`ensembl_protein_id`  varchar(40)  DEFAULT NULL,
`chrom` varchar (2)  DEFAULT NULL,
`strand` char (1)  DEFAULT NULL,
`hg19_exon_starts` text   DEFAULT NULL,
`hg19_exon_ends` text   DEFAULT NULL,
`hg19_cds_start` int(11)   DEFAULT NULL,
`hg19_cds_end` int(11)   DEFAULT NULL,
`sequence`  text   DEFAULT NULL,
  PRIMARY KEY (`uniprot_id`)
) ENGINE=MyISAM;
