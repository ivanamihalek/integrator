-- to create the table:
-- bash> mysql -u root blimps_environment < 06_homo_sapiens.sql 
-- to load dat, start mysql in local infile mode:
-- mysql --local-infile -uroot
-- then load a cleaned biogrid table from mysql shell
-- mysql>  load data local infile '/home/ivana/databases/biogrid/BIOGRID-ORGANISM-Homo_sapiens-3.4.140.tab2.clean.txt' into table biogrid_human_interactions ignore 1 lines;

-- if done as a preparation for daddy, follow with
-- rake db:populate:mark_problematic_biogrid
-- rake db:output:biogrid_edges

CREATE TABLE `biogrid_human_interactions` (
`biogrid_id` int(10) unsigned NOT NULL,
`entrez_gene_A` int(10) unsigned DEFAULT NULL,
`entrez_gene_B` int(10) unsigned DEFAULT NULL,
`biogrid_id_A` int(10) unsigned DEFAULT NULL,
`biogrid_id_B` int(10) unsigned DEFAULT NULL,
`systematic_name_A` varchar(40)   DEFAULT NULL,
`systematic_name_B` varchar(40)   DEFAULT NULL,
`official_symbol_A` varchar(40)   DEFAULT NULL,
`official_symbol_B` varchar(40)   DEFAULT NULL,
`synonyms_A`  varchar(300)   DEFAULT NULL,
`synonyms_B`  varchar(300)   DEFAULT NULL,
`experimental_system`  varchar(100)   DEFAULT NULL,
`experimental_system_type`  varchar(100)   DEFAULT NULL,
`author`  blob  DEFAULT NULL,
`pubmed_id` int(10) unsigned DEFAULT NULL,
`organism_A` int(10) unsigned DEFAULT NULL,
`organism_B` int(10) unsigned DEFAULT NULL,
`throughput`  varchar(40)   DEFAULT NULL,
`score`  varchar(40)    DEFAULT NULL,
`modification` varchar(100)   DEFAULT NULL,
`phenotypes`  varchar(100)   DEFAULT NULL,
`qualifications` blob   DEFAULT NULL,
`tags`  varchar(10)   DEFAULT NULL,
`source_database` varchar(40)   DEFAULT NULL,
`problematic` tinyint(1)  DEFAULT 1,
  PRIMARY KEY (`biogrid_id`),
  KEY `gene_a_idx` (`entrez_gene_A`),
  KEY `gene_b_idx` (`entrez_gene_B`)
) ENGINE=MyISAM;
