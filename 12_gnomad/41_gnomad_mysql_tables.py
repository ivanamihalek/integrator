#!/usr/bin/python3


from integrator_utils.python.mysql import *

##########################################
def main():

	mysql_conf_file = "/home/ivana/.tcga_conf"
	db     = connect_to_mysql(conf_file=mysql_conf_file)
	cursor = db.cursor()

	dbname = "gnomad"

	switch_to_db(cursor,dbname)

	check_and_drop_table(cursor, dbname, "meta")
	qry = "create table meta "
	qry += " (id int(11) NOT NULL AUTO_INCREMENT, "
	qry += " name varchar(255), value text, "
	qry += " PRIMARY KEY (id) )"
	search_db(cursor, qry, verbose=True)

	check_and_drop_table(cursor, dbname, "hotspots")
	qry = "create table hotspots "
	qry += " (id int(11) NOT NULL AUTO_INCREMENT, "
	qry += " chrom varchar(255), start int(11), end int(11), "
	qry += " PRIMARY KEY (id) )"
	search_db(cursor, qry, verbose=True)

	for chrom in [str(i) for i in range(1,23)] + ['X','Y']:
		table = "freqs_chr_" + chrom
		check_and_drop_table(cursor, dbname, table)
		print(table)
		#continue
		qry  = "CREATE TABLE %s " % table
		qry += " (id int(11) NOT NULL AUTO_INCREMENT, position int(11), reference varchar(400),  variant varchar(400), "
		qry += " variant_count int, total_count int, "
		qry += " afr_count int, afr_tot_count int, "
		qry += " amr_count int, amr_tot_count int, "
		qry += " asj_count int, asj_tot_count int, "
		qry += " eas_count int, eas_tot_count int, "
		qry += " fin_count int, fin_tot_count int, "
		qry += " nfe_count int, nfe_tot_count int, "
		qry += " oth_count int, oth_tot_count int, "
		qry += " sas_count int, sas_tot_count int, "
		qry += " homozygote_count int, "
		qry += " hotspot_id int, "
		qry += " PRIMARY KEY (id) )" # the default engine is InnoDB
		search_db(cursor, qry, verbose=True)
	cursor.close()
	db.close()
	return True

#########################################
if __name__ == '__main__':
	main()
