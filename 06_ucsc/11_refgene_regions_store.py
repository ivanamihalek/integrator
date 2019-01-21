#!/usr/bin/python3

from  integrator_utils.python.mysql import *
import os

#########################################
def main():

	print ("drop the ucsc.refgenes table - this script does not check for duplicate entries")
	exit(1)

	chromosomes = ["chr"+str(x) for x in range(1,23)] + ["chrX", "chrY"]
	assembly = "hg19"

	local_conf_file = "/home/ivana/.mysql_conf"
	ucsc_conf_file  = "/home/ivana/.ucsc_mysql_conf"

	for dependency in [local_conf_file, ucsc_conf_file]:
		if not os.path.exists(dependency):
			print(dependency, "not found")
			exit()

	local_db = connect_to_mysql(conf_file=local_conf_file)
	local_cursor = local_db.cursor()
	# autocommit is on by default, except when it is not
	search_db(local_cursor,"set autocommit=1")
	switch_to_db(local_cursor,'ucsc')

	ucsc_db     = connect_to_mysql(conf_file=ucsc_conf_file)
	ucsc_cursor = ucsc_db.cursor()
	switch_to_db(ucsc_cursor, assembly)
	#
	for chrom in chromosomes:
		print("downloading data for", assembly, chrom)
		qry  = "select name,  name2, strand, txStart, txEnd, exonCount, exonStarts,  exonEnds "
		qry += "from refGene "
		qry += "where chrom='%s' " % chrom
		qry += "and name like 'NM_%'"   # refseq says: NM_	mRNA	Protein-coding transcripts (usually curated)
		rows = search_db(ucsc_cursor,qry)
		print("loading ...")
		for row in rows:
			[gene_name, strand, txStart, txEnd, exonCount, exonStarts,  exonEnds] = row[1:]
			exonStarts = exonStarts.decode("utf-8").rstrip().rstrip(',')
			exonEnds = exonEnds.decode("utf-8").rstrip().rstrip(',')
			#print(gene_name, strand, txStart, txEnd, exonCount, exonStarts,  exonEnds)
			# store region
			# name, assembly, chromosome, strand, tx_tart, tx_end, exon_count, exon_starts,  exon_ends
			fields  = {'name':gene_name, 'assembly':assembly,
						'chromosome':chrom, 'strand':strand,
						'tx_start':txStart, 'tx_end':txEnd, 'exon_count':exonCount,
						'exon_starts':exonStarts,  'exon_ends':exonEnds}
			region_id = store_without_checking(local_cursor, 'refgenes', fields)
			if region_id<0:
				print("insert failure for", gene_name)
				exit()
		print ("\t{} done".format(chrom))
	ucsc_cursor.close()
	ucsc_db.close()


	local_cursor.close()
	local_db.close()
	return True


#########################################
if __name__ == '__main__':
	main()


