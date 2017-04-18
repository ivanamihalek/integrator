#!/usr/bin/python


from integrator_utils.mysql import *

##########################################
def main():
	db, cursor = connect()
	for chrom in [str(i) for i in range(1,23)] + ['X','Y']:
		table = "exac_freqs_chr_" + chrom
		#print table
		qry  = "CREATE TABLE %s " % table
		qry += " ( position int(11), reference text,  variants text, "
		qry += "  variant_counts text, total_count int,"
		qry += "  PRIMARY KEY (position) )" # the default engine is InnoDB
		search_db(cursor, qry, verbose=True)
	cursor.close()
	db.close()
	return True

#########################################
if __name__ == '__main__':
	main()
