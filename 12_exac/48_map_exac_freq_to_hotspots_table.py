#!/usr/bin/python

from integrator_utils.mysql import *
from time import time

#########################################
def main():

	db, cursor = connect()
	chroms = [str(i) for i in range(1,23)] + ['X','Y']
	chroms.reverse()

	for chrom in chroms:
		t0 = time()
		qry = "select id, start, end from exac_hotspots where chrom='%s' " % chrom
		for hotspot in  search_db(cursor, qry):
			[hotspot_id, start, end] = hotspot
			# find positions in the exac_freq table falling in this region
			# update their hotspot_id to the current value
			table = "exac_freqs_chr_" + chrom
			qry = "update %s " % table
			qry += "set hotspot_id=%d " % hotspot_id
			qry += "where %d<=position  and  position<=%d" % (start,end)
			search_db(cursor, qry)

	cursor.close()
	db.close()

	return
#########################################
if __name__ == '__main__':
	main()

