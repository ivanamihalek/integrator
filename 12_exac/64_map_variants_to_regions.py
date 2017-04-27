#!/usr/bin/python

# this builds on monogenic db:populate pipeline
# this is even wore than ruby

from integrator_utils.mysql import *
from time import time
##########################################
def main():

	db, cursor = connect()
	switch_to_db(cursor, 'monogenic_development')
	chromosomes = [str(x) for x in range(1,23)] + ['X']
	count = 0
	for chrom in chromosomes:
		print chrom
		t0 = time()
		rows = search_db(cursor,"select * from ensembl_gene_regions where chrom='%s'" %chrom)
		for row in rows:
			count += 1
			[ens_region_id, ens_id, symbol, chrom, strand, start, end] = row
			table = "chr%s_variant_annotations" % chrom
			qry = "select * from %s where %d<=position and position<=%d" %(table, start, end)
			rows2  = search_db(cursor, qry)
			if not rows2: continue
			for row2 in rows2:
				var_id =  row2[0]
				qry = "update %s set ensembl_gene_region_id=%d "  %(table, ens_region_id)
				qry += "where id=%d" % var_id
				search_db(cursor, qry)
			if not count%2: print " %d  %.2f min" % (count, (time()-t0)/60.0)
		print chrom, " done"

	cursor.close()
	db.close()

	return
#########################################
if __name__ == '__main__':
	main()

