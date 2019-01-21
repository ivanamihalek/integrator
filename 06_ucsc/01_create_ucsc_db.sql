
CREATE DATABASE ucsc;
USE ucsc;

--gene_name, assembly, chromosome, strand, txStart, txEnd, exonCount, exonStarts,  exonEnds
CREATE TABLE refgenes(
  id int  NOT NULL AUTO_INCREMENT PRIMARY KEY,
  name varchar(50) NOT NULL,
  assembly varchar(50) NOT NULL,
  chromosome varchar(50) NOT NULL,
  strand varchar(1) DEFAULT NULL,
  tx_start bigint  NOT NULL,
  tx_end bigint  NOT NULL,
  exon_count int NOT NULL,
  exon_starts longblob NOT NULL,
  exon_ends longblob NOT NULL
) ENGINE=InnoDB;

