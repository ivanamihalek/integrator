-- to create the table:
-- bash> mysql -u root blimps_environment < 02_uniprot_table.sql
-- to load data,
-- mysqlimport --local -u root -p blimps_development  uniprot_basic_infos.csv

CREATE TABLE `uniprot_basic_infos` (
`uniprot_id` varchar(40)   NOT NULL,
`gene_name` varchar(40)   DEFAULT NULL,
`ensembl_gene_id` varchar(40)   DEFAULT NULL,
`ec_number` varchar(40)   DEFAULT NULL,
`cofactors` varchar(250)   DEFAULT NULL,
`canonical_aa_length` int   DEFAULT NULL,
`full_name`  varchar(100)   DEFAULT NULL,
`tissue`  text   DEFAULT NULL,
`function` text  DEFAULT NULL,
`old_ids` text   DEFAULT NULL,
  PRIMARY KEY (`uniprot_id`),
  KEY `gene_name_idx` (`gene_name`)
) ENGINE=MyISAM;
