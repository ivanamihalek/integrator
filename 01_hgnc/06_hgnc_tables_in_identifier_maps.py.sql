CREATE TABLE uniprot_hgnc(
  uniprot_id varchar(20)  NOT NULL  PRIMARY KEY,
  hgnc_approved varchar(20),
  uniprot_old varchar(150)
) ENGINE=InnoDB;



