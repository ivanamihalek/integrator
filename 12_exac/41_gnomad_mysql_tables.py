#!/usr/bin/python


from integrator_utils.mysql import *

##########################################
def main():
	db, cursor = connect()
	for chrom in [str(i) for i in range(1,23)] + ['X','Y']:
		table = "gnomad_freqs_chr_" + chrom
		#qry = "DROP TABLE %s " % table
		#search_db(cursor, qry, verbose=True)
		print table
		qry  = "CREATE TABLE %s " % table
		qry += " (position int(11), reference text,  variants text, consequences text, "
		qry += " variant_counts text, total_count int, "
		qry += " afr_counts text, afr_tot_count int, "
		qry += " amr_counts text, amr_tot_count int, "
		qry += " eas_counts text, eas_tot_count int, "
		qry += " fin_counts text, fin_tot_count int, "
		qry += " nfe_counts text, nfe_tot_count int, "
		qry += " oth_counts text, oth_tot_count int, "
		qry += " sas_counts text, sas_tot_count int, "
		qry += " hotspot_id int, "
		qry += " PRIMARY KEY (position) )" # the default engine is InnoDB
		search_db(cursor, qry, verbose=True)
	cursor.close()
	db.close()
	return True

#########################################
if __name__ == '__main__':
	main()
