-- to create the table:
-- bash> mysql -u root blimps_environment < 23_metacyc_edges.py

CREATE TABLE `metacyc_edges` (
`id` int(11)   NOT NULL AUTO_INCREMENT,
`from_field` varchar(40)   NOT NULL,
`from_id` varchar(40)   NOT NULL,
`to_field` varchar(40)   NOT NULL,
`to_id` varchar(40)   NOT NULL,
`via`  varchar(40)   NOT NULL,
  PRIMARY KEY (`id`),
  KEY `from_field_idx` (`from_field`)
) ENGINE=MyISAM;
