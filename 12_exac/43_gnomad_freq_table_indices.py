#!/usr/bin/python
import math

from integrator_utils.python.mysql import *

##########################################
def main():

	db     = connect_to_mysql()
	cursor = db.cursor()

	switch_to_db(cursor,"gnomad")

	# to put an index on a column it must have a fixed number of chars, it cannot be text
	#  find the longest reference and variant and dimension the fields accordingly

	# absmax = -1
	# for chrom in [str(i) for i in range(1,23)] + ['X']:
	# 	table = "gnomad_freqs_chr_" + chrom
	# 	qry = "select max(char_length(reference)) from %s" % table
	# 	maxr = search_db(cursor, qry, verbose=True)[0][0]
	# 	if maxr> absmax: absmax=maxr
	# 	qry = "select max(char_length(variant)) from %s" % table
	# 	maxv = search_db(cursor, qry, verbose=True)[0][0]
	# 	if maxv> absmax: absmax=maxv
	#
	# print absmax
	# # round to a 100
	# var_ref_size = int(100*math.ceil(float(absmax)/100))
	# print "new col size:", var_ref_size
	#
	# for chrom in [str(i) for i in range(1,23)] + ['X']:
	# 	table = "gnomad_freqs_chr_" + chrom
	# 	qry = "alter table  %s change reference reference varchar(%d)" % (table,var_ref_size)
	# 	search_db(cursor, qry, verbose=True)
	# 	qry = "alter table  %s change variant variant varchar(%d)" % (table,var_ref_size)
	# 	search_db(cursor, qry, verbose=True)

	for chrom in [str(i) for i in range(1,23)] + ['X']:
		table = "gnomad_freqs_chr_" + chrom
		qry = "create index variant_id on %s (position,reference,variant)" % table
		search_db(cursor, qry, verbose=True)



	cursor.close()
	db.close()
	return True

#########################################
if __name__ == '__main__':
	main()
