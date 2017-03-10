#!/usr/bin/python
# this assumes we have read in the knownGene from ucsc
# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/
#
# mysql -u root ucsc < /databases/ucsc/knownGene.sql
#
# mysql -u root --local-infile=1
# use ucsc
# LOAD DATA LOCAL INFILE '/databases/ucsc/knownGene.txt.gz' INTO TABLE knownGene;
# GRANT ALL PRIVILEGES ON ucsc.* TO  cookiemonster@localhost;

# it quickly turned out that knowGene is using some antiquated uniprot (!) identifiers
# switching to

# now export exon and gene regions to bedfile, to be use by seqmule in reporting overall coverage



from integrator_utils.mysql import *

def  process_exons(exonStarts, exonEnds):
	exon_ranges = set()
	if exonStarts[-1] == ',':  exonStarts = exonStarts[0:-1]
	if exonEnds[-1] == ',':    exonEnds   = exonEnds[0:-1]
	starts = exonStarts.split(',')
	ends   = exonEnds.split(',')
	if len(starts) != len(ends):
		print "starts/ends length mismatch"
		exit(1)
	for i in range(len(starts)):
		exon_ranges.add(starts[i] + "_" + ends[i])
	return exon_ranges

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
def remove_duplicate_regions(exon_coordinates):

	processed_coords = []
	for coords in sorted(list(exon_coordinates), key=lambda x: int(x.split("_")[0])):
		[a_new, b_new] = coords.split("_")
		placed = False
		for i in range(len(processed_coords)):
			[a_old,b_old] = processed_coords[i]
			if b_new < a_old:
				# insert
				new_processed_coords = processed_coords[:i] + [[a_new, b_new]] +  processed_coords[i:]
				processed_coords = new_processed_coords
				placed = True
				break
			elif a_new > b_old:
				#move on
				pass
			else:
				a = min(a_old,a_new)
				b = max(b_old,b_new)
				processed_coords[i] = [a,b]
				placed = True
				break
		if not placed: processed_coords.append([a_new, b_new])

	return processed_coords

##########################################
def main():
	db, cursor = connect()
	switch_to_db(cursor, "ucsc")
	gene_ids = [] # keep to maintain the original order
	general_info = {}
	exon_coordinates = {}
	for chrom in [str(i) for i in range(1,23)] +["X","Y"]:
		# name is Ensembl transcript ID, and name2 is ensembl gene id
		qry  = "select chrom, strand, exonStarts, exonEnds,  name2 "
		qry += "from ensGene where chrom='chr%s' order by txStart " % chrom
		rows = search_db(cursor, qry)
		for row in rows:
			[chrom, strand,  exonStarts, exonEdns,  gene_id] = row
			if not gene_id in gene_ids:
				gene_ids.append(gene_id)
				general_info[gene_id] = [chrom, strand]
				exon_coordinates[gene_id] = set()
			exon_coordinates[gene_id] |= process_exons(exonStarts, exonEdns)

	switch_to_db(cursor, "blimps_development")

	for gene_id  in gene_ids:
		symbol = get_symbol(cursor,gene_id)
		if not symbol: continue
		chrom  = general_info[gene_id][0].replace('chr','')
		strand = general_info[gene_id][1]
		processed_coords = remove_duplicate_regions(exon_coordinates[gene_id])
		number_of_exons = len(processed_coords)
		exon_ct = 0
		for [start, stop]   in processed_coords:
			if strand=='-':
				name = "%s_exon_%d" %(symbol, number_of_exons-exon_ct)
			else:
				name = "%s_exon_%d" %(symbol, exon_ct + 1)
			print "\t".join([chrom, start, stop, name])
			exon_ct += 1


	cursor.close()
	db.close()
	return True

#########################################
if __name__ == '__main__':
	main()
