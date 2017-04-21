#!/usr/bin/python

from integrator_utils.mysql import *

##########################################
def main():
	db, cursor = connect() # this puts me to blimps_*  directly

	# find genes of interest (the ones that are mapped onto a mendelian disease)
	qry  = "select distinct(ensembl_gene_id) from omim_genemaps "
	qry += "where ensembl_gene_id is not null and ensembl_gene_id!='' "
	qry += "and phenotypes is not null and phenotypes!=''"
	rows = search_db(cursor, qry)
	switch_to_db(cursor, "monogenic_development")
	not_found = 0
	# annovar wants to  see 8 fields
	print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
	for row in rows:
		ensembl_gene_id = row[0]
		if len(ensembl_gene_id) == 0: continue
		# for each gene find the region
		qry = "select chrom, start, end from ensembl_gene_regions where ensembl_gene_id='%s'" % ensembl_gene_id
		rows2 = search_db(cursor, qry)
		if not rows2 or len(rows2)<1:
			not_found += 1
			continue
		if len(rows2)>1:
			print "two regions for %s (?)" % ensembl_gene_id
			exit()
		for row2 in rows2:
			[chrom, start, end] = row2
			# for each region find the variants
			qry  = "select * from blimps_development.exac_freqs_chr_%s " % chrom
			qry += "where %d<=position and position<= %d" % (int(start), int(end))
			rows3 = search_db(cursor, qry)
			if not rows3: continue
			for row3 in rows3:
				[pos, ref, variants, variant_counts, total_count] = row3
				vars = variants.split(",")
				for alt in vars:
					print "\t".join([chrom, str(pos), "-", ref, alt, "1000","PASS","-"])
					pass

	#print "total", len(rows), "not mapped:", not_found

	cursor.close()
	db.close()

	return
#########################################
if __name__ == '__main__':
	main()
