CREATE TABLE `omim_genemaps` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `mim_number` int(11) DEFAULT NULL,
  `gene_symbols` text,
  `gene_name` varchar(255) DEFAULT NULL,
  `approved_symbol` varchar(255) DEFAULT NULL,
  `entrez_gene_id` varchar(255) DEFAULT NULL,
  `ensembl_gene_id` varchar(255) DEFAULT NULL,
  `phenotypes` text,
  `mouse_gene_symbol` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=15447 DEFAULT CHARSET=latin1;
