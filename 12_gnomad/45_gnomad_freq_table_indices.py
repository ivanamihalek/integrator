#!/usr/bin/python3
import math

from integrator_utils.python.mysql import *

##########################################
def main():

	mysql_conf_file = "/home/ivana/.tcga_conf"
	for dep in [ mysql_conf_file]:
		if not os.path.exists(dep):
			print(dep, "not found")
	db     = connect_to_mysql(conf_file=mysql_conf_file)
	cursor = db.cursor()

	switch_to_db(cursor,"gnomad")

	for chrom in [str(i) for i in range(1,23)] + ['X']:
		table = "freqs_chr_" + chrom
		#qry = "create index variant_id on %s (position,reference,variant)" % table
		qry = "create index position_id on %s (position)" % table
		search_db(cursor, qry, verbose=True)

	cursor.close()
	db.close()
	return True

#########################################
if __name__ == '__main__':
	main()
