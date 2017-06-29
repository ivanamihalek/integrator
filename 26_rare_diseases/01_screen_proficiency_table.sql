-- data from Newborn Screening Quality Assurance Program, annual report 2014
-- to create the table:
-- bash> mysql -u root blimps_environment < 02_screen_proficiency.sql
-- mysqlimport --local -u root -p blimps_development  screen_proficiencies.csv

CREATE TABLE `screen_proficiencies` (
`screen` varchar(255)   NOT NULL,
`number_of_positive_specimens` int   DEFAULT NULL,
`pct_false_negative_error` float  DEFAULT NULL,
`number_of_negative_specimens` int   DEFAULT NULL,
`pct_false_positive_error` float  DEFAULT NULL,
  PRIMARY KEY (`screen`)
) ENGINE=MyISAM;
