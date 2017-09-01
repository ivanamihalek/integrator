#! /usr/bin/python

# in particular, find the cases where uniprot sequence is substantially longer than the ensembl translation
from integrator_utils.mysql import *
import os,  subprocess, requests, re


###################################################
def main():

	db, cursor = connect()
	# make sure that the uniprot_seqs table has cds columns
	for col in  ['ensembl_coding_sequence', 'ensembl_cds_translation']:
		if 'coding' in col:
			coltype= "mediumtext"
		else:
			coltype= "text"
		check_or_make_column (cursor, 'monogenic_development', 'uniprot_seqs', col, coltype)

	# limit ourselves to iems related genes
	qry = "select ensembl_gene_id from blimps_development.omim_genemaps where inborn_error_of_metabolism=1"
	ret = search_db(cursor, qry)
	total = 0
	in_agreement = 0
	for line in ret:
		ensembl_gene_id = line[0]
		uniprot_id = get_uniprot_id (cursor, ensembl_gene_id)
		qry = "select * from monogenic_development.uniprot_seqs where uniprot_id= '%s' " % uniprot_id
		[uniprot_id, ensembl_transcript_id, ensembl_protein_id, chrom, strand,\
		hg19_exon_starts, hg19_exon_ends, hg19_cds_start, hg19_cds_end, \
		sequence, cons_in_verts_map_to_uni, ensembl_coding_sequence, \
		ensembl_cds_translation, uni2ens_cigar, cons_in_verts_map_to_ens] = search_db(cursor,qry)[0]
		if not sequence:
			print "no uniprot seqeunce found for", uniprot_id
			exit()
		if not ensembl_cds_translation: continue
		length_ratio =   float(len(ensembl_cds_translation))/len(sequence)
		if length_ratio < 0.9:
			print "%s %s  %4d %4d   %8.3f" %(uniprot_id,  ensembl_gene_id, len(ensembl_cds_translation), len(sequence), length_ratio)

	return


#########################################
if __name__ == '__main__':
	main()
