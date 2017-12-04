-- to create the table:
-- bash> mysql -u root blimps_environment < 01_bgee.sql
-- Gene ID	"Gene name"	Anatomical entity ID	"Anatomical entity name"	Expression	Call quality	Expression rank
CREATE TABLE `expression_bgee` (
     `ensembl_gene_id` varchar(40)   NOT NULL,
    `gene_name` varchar(40)   DEFAULT NULL,
    `anatomical_entity_id` varchar(40)   DEFAULT NULL,
    `anatomical_entity_name`  text   DEFAULT NULL,
    `expression` varchar(40)   DEFAULT NULL,
    `call_quality` varchar(40)   DEFAULT NULL,
    `expression_rank` varchar(40)   DEFAULT NULL,
    `id` int(11) NOT NULL AUTO_INCREMENT,
     PRIMARY KEY (`id`),
      KEY `ensembl_gene_id` (`ensembl_gene_id`),
      KEY `gene_name` (`gene_name`)
) ENGINE=MyISAM;

-- it is important for the autoincremented id to be the last,
-- because it is not  preslent in the original table
-- delete quotes from the input file:
-- sed s'/\"//g' Homo_sapiens_expr_simple.tsv -i
-- to load data, start mysql in local infile mode:
-- mysql --local-infile -uroot
-- then load bgee table from mysql shell
-- mysql>  load data local infile '/databases/bgee/Homo_sapiens_expr_simple.tsv' into table expression_bgee ignore 1 lines;
