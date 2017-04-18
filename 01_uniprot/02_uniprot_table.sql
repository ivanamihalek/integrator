-- to create the table:
-- bash> mysql -u root blimps_environment < 02_uniprot_table.sql
-- to load dat, start mysql in local infile mode:
-- mysql --local-infile -uroot
-- then load a cleaned biogrid table from mysql shell
-- mysql>  load data local infile 'uniprot.csv' into table uniprot_basic_infos ignore 1 lines;
-- the horrible table name is to comply with rails pluralization convention

CREATE TABLE `uniprot_basic_infos` (
`uniprot_id` varchar(40)   NOT NULL,
`gene_name` varchar(40)   DEFAULT NULL,
`full_name`  varchar(100)   DEFAULT NULL,
`tissue`  blob   DEFAULT NULL,
`function` blob  DEFAULT NULL,
  PRIMARY KEY (`uniprot_id`),
  KEY `gene_name_idx` (`gene_name`)
) ENGINE=MyISAM;
