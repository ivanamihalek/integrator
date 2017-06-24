-- bash> mysql -u root monogenic_development -p < 18_model_elements_table.sql
-- ENSG00000137992
-- 123456789012345
CREATE TABLE `model_elements` (
    `id`            int(11)   NOT NULL AUTO_INCREMENT,
    `gene_symbol`       varchar(10)   NOT NULL,
    `ensembl_gene_id`   varchar(20)   NOT NULL,
    `uniprot_id`        varchar(10) DEFAULT NULL,
    `ec_number`         varchar(20) DEFAULT NULL,
    `main_model`       varchar(100) DEFAULT NULL,
    `pdb_chain`         varchar(10)   NOT  NULL,
    `pdb_pct_identical_to_uniprot` TINYINT DEFAULT NULL,
    `pdb_ligand`        varchar(3)  NOT NULL,
    `metacyc_ligand`   varchar(100) DEFAULT NULL,
    `ligand_tanimoto`  float DEFAULT NULL,
      PRIMARY KEY (`id`),
      KEY `gene_symbol_idx` (`gene_symbol`)
) ENGINE=MyISAM;
