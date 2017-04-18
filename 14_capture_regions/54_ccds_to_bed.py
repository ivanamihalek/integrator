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

import re

##########################################
def remove_duplicate_regions(exon_coordinates):

	processed_coords = []
	for coords in sorted( list(exon_coordinates), key=lambda x: int(x.split("-")[0])):
		[a_new, b_new] = [int(x) for x in coords.split("-")]
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
# ccds current txt columns:
# 0               1          2       3        4       5               6          7         8       9              10
#chromosome	nc_accession	gene	gene_id	ccds_id	ccds_status	cds_strand	cds_from	cds_to	cds_locations	match_type
def main():

	infile = open("/databases/ccds/15/CCDS.current.txt","r")

	gene_ids = {} # keep to maintain the original order
	general_info = {}
	exon_coordinates = {}
	for line in infile:
		if line[0]=='#': continue
		[chrom, nc_accession, gene, gene_id, ccds_id, ccds_status,
		 cds_strand, cds_from, cds_to, cds_locations, match_type] = line.rstrip().split("\t")
		if ccds_status!='Public': continue
		if not gene_ids.has_key(chrom):
			gene_ids[chrom] = []
		if not gene in gene_ids[chrom]:
			gene_ids[chrom].append(gene)
			general_info[gene] = [chrom, cds_strand, int(cds_from)]
			exon_coordinates[gene] = set()
		exon_coordinates[gene] |= set(re.sub(r'[\[\]\s]','', cds_locations).rstrip(',').split(','))

	# sort by the starting position on the chromosome
	for chrom, chrom_gene_ids in gene_ids.iteritems():
		gene_ids[chrom] = sorted(chrom_gene_ids, key=lambda x: general_info[x][2])

	for chrom in [str(x) for x in range(1,23)] + ['X','Y']:
		for gene_id  in gene_ids[chrom]:
			chrom  = general_info[gene_id][0].replace('chr','')
			strand = general_info[gene_id][1]
			processed_coords = remove_duplicate_regions(exon_coordinates[gene_id])
			number_of_exons = len(processed_coords)
			exon_ct = 0
			for [start, stop] in processed_coords:
				if strand=='-':
					name = "%s_exon_%d" %(gene_id, number_of_exons-exon_ct)
				else:
					name = "%s_exon_%d" %(gene_id, exon_ct + 1)
				print "\t".join([chrom, str(start), str(stop), name])
				exon_ct += 1
	return True

#########################################
if __name__ == '__main__':
	main()
