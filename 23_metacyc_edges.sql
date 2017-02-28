-- to create the table:
-- bash> mysql -u root blimps_environment < 02_uniprot_table.sql
-- to load dat, start mysql in local infile mode:
-- mysql --local-infile -uroot
-- then load a cleaned biogrid table from mysql shell
-- mysql>  load data local infile 'uniprot.csv' into table uniprot_basic_infos ignore 1 lines;
-- the horrible table name is to comply with rails pluralization convention

CREATE TABLE `metacyc_edges` (
`id` int(11)   NOT NULL AUTO_INCREMENT,
`from_field` varchar(40)   NOT NULL,
`from_id` varchar(40)   NOT NULL,
`to_field` varchar(40)   NOT NULL,
`to_id` varchar(40)   NOT NULL,
`via`  varchar(40)   NOT NULL,
  PRIMARY KEY (`id`),
  KEY `from_field_idx` (`from_field`)
) ENGINE=MyISAM;
