-- bash> mysql -u root monogenic_developmen -p < 14_pdb_ligands_table.sql
CREATE TABLE `pdb_ligands` (
    `id`            int(11)   NOT NULL AUTO_INCREMENT,
    `pdb_id`        char(4)   NOT NULL,
    `chemical_id`   char(3)   NOT NULL,
    `ligand_type`    varchar(40) DEFAULT NULL,
    `chemical_name` varchar(200) DEFAULT NULL,
    `smiles` text DEFAULT NULL,
      PRIMARY KEY (`id`),
      KEY `pdbid_idx` (`pdb_id`)
) ENGINE=MyISAM;
