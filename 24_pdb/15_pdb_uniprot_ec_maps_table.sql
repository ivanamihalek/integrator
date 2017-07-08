-- bash> mysql -u root monogenic_developmen -p < 15_pdb_uniprot_ec_maps_table.sql
CREATE TABLE `pdb_uniprot_ec_maps` (
    `pdb_chain`        char(5)   NOT NULL,
    `uniprot_id`   varchar(10)   NOT NULL,
    `ec_numbers`    varchar(200) DEFAULT NULL,
    `cofactors`    varchar(250),
      PRIMARY KEY (`pdb_chain`)
) ENGINE=MyISAM;
