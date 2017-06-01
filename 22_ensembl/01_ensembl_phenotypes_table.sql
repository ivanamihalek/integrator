-- bash> mysql -u root blimps_environment < 02_ensembl_phenotype_table.sql
-- ensembl_phenotypes.csv made  by awking from  Homo_sapiens_phenotype_associated.vcf
-- available from ftp://ftp.ensembl.org/pub/grch37/update/variation/vcf/homo_sapiens/
-- mysqlimport --local -u root -p blimps_development  ensembl_phenotypes.csv

CREATE TABLE `ensembl_phenotypes` (
`chrom`  char(2)   NOT NULL,
`position` int(11)   NOT NULL,
`ref`  varchar(100)   DEFAULT NULL,
`alt`  varchar(100)   DEFAULT NULL,
`phenotype`  text   DEFAULT NULL,
  KEY `chrom_pos_idx` (`chrom`,`position`)
) ENGINE=MyISAM;
