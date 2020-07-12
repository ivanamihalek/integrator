CREATE TABLE meta (
  id int  NOT NULL AUTO_INCREMENT PRIMARY KEY,
  description varchar (150) NOT NULL,
  payload text NOT NULL
) ENGINE=InnoDB;


--
CREATE TABLE kegg_human(
  kegg_id int  NOT NULL  PRIMARY KEY,
  uniprot varchar(50),
  uniprot_old varchar(150),
  kegg_pathways text
) ENGINE=InnoDB;


CREATE TABLE kegg_pathway_name(
   kegg_pathway_id varchar(5)   NOT NULL  PRIMARY KEY,
   name varchar(150)
) ENGINE=InnoDB;