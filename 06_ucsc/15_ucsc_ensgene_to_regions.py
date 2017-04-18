#!/usr/bin/python

from integrator_utils.mysql import *


##########################################
def get_symbol(cursor,ensembl_gene_id):

	qry = "select symbol from genes where ensembl_gene_id='%s'" % ensembl_gene_id
	rows = search_db(cursor,qry)
	if not rows or rows[0][0]==None or rows[0][0]=='' or rows[0][0].lower()=='error':
		return None
	if len(rows) > 1:
		print "multiple gene reurns for", ensembl_gene_id
		exit()

	return rows[0][0]


##########################################
def main():
	db, cursor = connect()

	for chrom in [str(i) for i in range(1,23)] +["X","Y"]:
		switch_to_db(cursor, "ucsc")
		# name is Ensembl transcript ID, and name2 is ensembl gene id
		gene_ids = [] # keep to maintain the original order
		min_start = {}
		max_end   = {}
		chromosome  = {}
		strand = {}
		qry  = "select chrom, strand, exonStarts, exonEnds,  name2 "
		qry += "from ensGene where chrom='chr%s'  " % chrom
		rows = search_db(cursor, qry)
		for row in rows:
			[chrm, strnd,  exonStarts, exonEnds,  gene_id] = row
			start = exonStarts.split(",")[0]
			end = exonEnds.split(",")[-1]
			if end=="": end = exonEnds.split(",")[-2]
			if not gene_id in gene_ids:
				gene_ids.append(gene_id)
				min_start[gene_id] = start
				max_end[gene_id] =  end
				chromosome[gene_id] = chrm.upper().replace("CHR","")
				strand[gene_id] = strnd
			else:
				if min_start[gene_id]> start: min_start[gene_id] = start
				if max_end[gene_id]<end:  max_end[gene_id] = end

		switch_to_db(cursor, "monogenic_development")
		for gene_id in gene_ids:
			fixed_fields = {"ensembl_gene_id":gene_id}
			update_fields = {"start":min_start[gene_id], "end":max_end[gene_id],
							 "chrom":chromosome[gene_id], "strand": strand[gene_id]}
			qry  = " select  symbol from blimps_development.genes where ensembl_gene_id='%s'" % gene_id
			rows = search_db(cursor, qry)
			if rows:
				if len(rows) > 2:
					print "more than 2 entries for a single ENSG"
					exit()
				update_fields["symbol"] = rows[0][0]
			store_or_update(cursor, "ensembl_gene_regions", fixed_fields, update_fields)

	cursor.close()
	db.close()
	return True

#########################################
if __name__ == '__main__':
	main()
